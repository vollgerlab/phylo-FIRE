#!/usr/bin/env python
"""Concatenate per-peak ``re_variants`` tables and pick which peaks to plot.

Reads every ``{indir}/{peak}.tsv`` for the run's peaks, concatenates
them into one parquet (the unified deliverable) plus a gzipped TSV, and writes
three plot manifests that drive data-driven plotting:

  * ``plot_drivers.txt``             — peaks with >=1 selected driver variant.
  * ``plot_heritable_no_driver.txt`` — H2 identified (RE or region) but 0 drivers.
  * ``plot_random.txt``              — N seeded-random peaks (any category).

Peak-level fields (``n_variants_selected``, ``h2_re_ok``, ``h2_region_ok``) are
repeated on every row of a peak, so classification reads the first row per peak.
Called in-process by re_variants_all.py.
"""
import gzip
import random
from pathlib import Path

import polars as pl

# Pin every column's type so per-file inference can't break the vertical concat.
# Numeric/bool columns stay numeric even when all blank (attribution_mode=none,
# non-carrier rows). The rest are Utf8: comma-list columns infer Float64 for a
# single value ("0.55") but String for a list ("0.1,0.55").
NUMERIC_BOOL = {
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


def read_tables(indir: Path, peaks: list[str]) -> pl.DataFrame:
    """Concatenate ``{indir}/{peak}.tsv`` in manifest order with pinned types."""
    first = next(
        (indir / f"{p}.tsv" for p in peaks if (indir / f"{p}.tsv").exists()), None
    )
    header = (
        pl.read_csv(first, separator="\t", n_rows=0).columns if first is not None else []
    )
    schema = {c: NUMERIC_BOOL.get(c, pl.Utf8) for c in header}

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
    return pl.concat(frames, how="vertical") if frames else pl.DataFrame()


def classify_peaks(table: pl.DataFrame, random_n: int, seed: int):
    """(drivers, heritable, random, present) peak lists for the plot manifests.

    One representative row per peak carries the peak-level fields. The random
    sample is a NEGATIVE control: peaks that are neither drivers nor
    heritable-no-driver (no selected driver, no identified H²; the epigenetic /
    not-sequence-driven majority), so it deliberately excludes them. Seeded."""
    drivers, heritable, present = [], [], []
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
    rng = random.Random(seed)
    exclude = set(drivers) | set(heritable)
    pool = sorted(p for p in present if p not in exclude)
    rand = sorted(rng.sample(pool, min(random_n, len(pool)))) if pool else []
    return drivers, heritable, rand, present


def aggregate(indir: Path, peaks: list[str], out_parquet: Path, out_tsv: Path,
              manifest_dir: Path, random_n: int = 100, seed: int = 42) -> None:
    """Write the unified parquet + gzipped TSV and the three plot manifests."""
    table = read_tables(Path(indir), peaks)
    Path(out_parquet).parent.mkdir(parents=True, exist_ok=True)
    table.write_parquet(out_parquet)
    with gzip.open(out_tsv, "wb") as fh:
        table.write_csv(fh, separator="\t")

    drivers, heritable, rand, present = classify_peaks(table, random_n, seed)
    mdir = Path(manifest_dir)
    mdir.mkdir(parents=True, exist_ok=True)
    (mdir / "plot_drivers.txt").write_text("\n".join(sorted(drivers)) + "\n")
    (mdir / "plot_heritable_no_driver.txt").write_text("\n".join(sorted(heritable)) + "\n")
    (mdir / "plot_random.txt").write_text("\n".join(rand) + "\n")

    print(f"[aggregate] {table.height} rows over {len(present)} peaks → "
          f"{len(drivers)} drivers, {len(heritable)} heritable-no-driver, "
          f"{len(rand)} random-background")
