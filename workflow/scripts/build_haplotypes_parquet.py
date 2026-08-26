#!/usr/bin/env python3
"""Backfill tables/haplotypes.parquet from an existing per-peak fastas/ dir.

For runs whose Stage 1 predates the two-stage split (so `fastas/*.fa` + `*.tsv`
are on disk but no consolidated parquet exists). Reads every peak's .fa (seq) +
.tsv (actuation/coverage/asm coords) and writes the one parquet Stage 2 reads.

Reads are parallelized (I/O-latency-bound on Lustre: ~2x46K tiny files) and
progress is printed, so a big run doesn't look hung.

    pixi run python workflow/scripts/build_haplotypes_parquet.py \
        --fastas-dir results/<run>/fastas \
        --out results/<run>/tables/haplotypes.parquet [--threads 32]
"""
import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

COLS = ["peak_id", "sample", "seq", "actuation", "coverage",
        "fire_coverage", "asm_chr", "asm_start", "asm_end"]


def read_fa(p: Path) -> dict[str, str]:
    d: dict[str, list[str]] = {}
    name = None
    for ln in p.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            d[name] = []
        elif name is not None:
            d[name].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


def read_peak(fa: Path):
    """Return list of per-hap tuples for one peak, or None if no .tsv."""
    tsv = fa.with_name(f"{fa.stem}.tsv")
    if not tsv.exists():
        return None
    seqs = read_fa(fa)
    lines = tsv.read_text().splitlines()
    ix = {h: i for i, h in enumerate(lines[0].split("\t"))}
    out = []
    for ln in lines[1:]:
        f = ln.split("\t")
        s = f[ix["sample"]]
        out.append((fa.stem, s, seqs.get(s, ""), f[ix["actuation"]],
                    f[ix["coverage"]], f[ix["fire_coverage"]],
                    f[ix["asm_chr"]], f[ix["asm_start"]], f[ix["asm_end"]]))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastas-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    fas = sorted(args.fastas_dir.glob("*.fa"))
    print(f"[haplotypes] reading {len(fas)} peaks from {args.fastas_dir} "
          f"with {args.threads} threads ...", file=sys.stderr)

    rows: dict[str, list] = {c: [] for c in COLS}
    done = 0
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        for res in ex.map(read_peak, fas):
            done += 1
            if res:
                for tup in res:
                    for k, v in zip(COLS, tup):
                        rows[k].append(v)
            if done % 1000 == 0:
                print(f"[haplotypes]   {done}/{len(fas)} peaks "
                      f"({len(rows['peak_id'])} rows)", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    (
        pl.DataFrame(rows, schema={c: pl.Utf8 for c in COLS})
        .with_columns(
            pl.col("actuation").cast(pl.Float64, strict=False),
            pl.col("coverage").cast(pl.Float64, strict=False),
            pl.col("fire_coverage").cast(pl.Float64, strict=False),
        )
        .write_parquet(args.out)
    )
    print(f"[haplotypes] wrote {len(rows['peak_id'])} rows over {len(fas)} "
          f"peaks -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
