#!/usr/bin/env python
"""Concatenate per-peak ``re_variants`` tables and pick which peaks to plot.

Reads every ``{indir}/{peak}.tsv`` listed in the run manifest, concatenates
them into one parquet (the unified deliverable) plus a gzipped TSV, and writes
three plot manifests that drive data-driven plotting:

  * ``plot_drivers.txt``             — peaks with >=1 selected driver variant.
  * ``plot_heritable_no_driver.txt`` — H2 identified (RE or region) but 0 drivers.
  * ``plot_random.txt``              — N seeded-random peaks (any category).

Peak-level fields (``n_variants_selected``, ``h2_re_ok``, ``h2_region_ok``) are
repeated on every row of a peak, so classification reads the first row per peak.
"""
import argparse
import gzip
import random
from pathlib import Path

import polars as pl


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True)
    ap.add_argument("--peaks", required=True, help="run manifest, one peak_id per line")
    ap.add_argument("--out-parquet", required=True)
    ap.add_argument("--out-tsv", required=True, help="gzipped TSV path")
    ap.add_argument("--manifest-dir", required=True)
    ap.add_argument("--random-n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    indir = Path(args.indir)
    peaks = [p.strip() for p in Path(args.peaks).read_text().splitlines() if p.strip()]

    # Pin EVERY column's type up front so per-file inference can't disagree
    # across peaks (concat how="vertical" requires identical schemas). Known
    # numeric/bool columns get their type; read_csv reads "" as null for those,
    # so a peak with an all-blank numeric column (e.g. attribution_mode=none →
    # blank h2_re_*, non-carrier rows → blank n_other_diffs) stays numeric.
    # EVERY OTHER column is forced to Utf8 — critically the comma-list columns
    # (act_values, samples, driver_ids/positions/betas/p), which infer Float64
    # in a peak whose rows each hold a single value ("0.55") and String in a peak
    # with multiple ("0.1,0.55"); that mismatch is what breaks the concat.
    numeric_bool = {
        "gc_frac": pl.Float64, "act_mean": pl.Float64, "act_sd": pl.Float64,
        "act_min": pl.Float64, "act_max": pl.Float64, "total_cov": pl.Float64,
        "driver_delta": pl.Float64,
        "h2_re_mean": pl.Float64, "h2_re_lower": pl.Float64, "h2_re_upper": pl.Float64,
        "h2_region_mean": pl.Float64, "h2_region_lower": pl.Float64,
        "h2_region_upper": pl.Float64,
        "re_len": pl.Int64, "n_haps": pl.Int64, "n_other_diffs": pl.Int64,
        "n_unique_seqs": pl.Int64, "n_variants_selected": pl.Int64,
        "has_driver": pl.Boolean, "clean_pair": pl.Boolean,
        "h2_re_ok": pl.Boolean, "h2_region_ok": pl.Boolean,
    }
    first = next(
        (indir / f"{p}.tsv" for p in peaks if (indir / f"{p}.tsv").exists()), None
    )
    header = (
        pl.read_csv(first, separator="\t", n_rows=0).columns if first is not None else []
    )
    schema = {c: numeric_bool.get(c, pl.Utf8) for c in header}

    frames, missing = [], 0
    for peak in peaks:
        f = indir / f"{peak}.tsv"
        if not f.exists():
            missing += 1
            continue
        df = pl.read_csv(f, separator="\t", schema_overrides=schema)
        if df.height:
            frames.append(df)
    if missing:
        print(f"[aggregate] WARNING: {missing}/{len(peaks)} peak TSVs missing "
              f"(failed/incomplete?) — omitted from the table", flush=True)

    table = pl.concat(frames, how="vertical") if frames else pl.DataFrame()

    Path(args.out_parquet).parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(args.out_parquet)
    with gzip.open(args.out_tsv, "wb") as fh:
        table.write_csv(fh, separator="\t")

    # ── classify peaks for plotting ───────────────────────────────────────────
    # One representative row per peak carries the peak-level fields.
    drivers, heritable = [], []
    present = []
    if table.height:
        per_peak = table.group_by("peak_id").first()
        for r in per_peak.iter_rows(named=True):
            pid = r["peak_id"]
            present.append(pid)
            n_sel = r.get("n_variants_selected") or 0
            if n_sel and int(n_sel) >= 1:
                drivers.append(pid)
            elif r.get("h2_re_ok") or r.get("h2_region_ok"):
                heritable.append(pid)

    # Truly-random sample across all peaks-with-data (may overlap drivers/heritable;
    # its own subdir, so a duplicate render is harmless).
    rng = random.Random(args.seed)
    pool = sorted(present)
    rand = sorted(rng.sample(pool, min(args.random_n, len(pool)))) if pool else []

    mdir = Path(args.manifest_dir)
    mdir.mkdir(parents=True, exist_ok=True)
    (mdir / "plot_drivers.txt").write_text("\n".join(sorted(drivers)) + "\n")
    (mdir / "plot_heritable_no_driver.txt").write_text("\n".join(sorted(heritable)) + "\n")
    (mdir / "plot_random.txt").write_text("\n".join(rand) + "\n")

    print(f"[aggregate] {table.height} rows over {len(present)} peaks → "
          f"{len(drivers)} drivers, {len(heritable)} heritable-no-driver, "
          f"{len(rand)} random")


if __name__ == "__main__":
    main()
