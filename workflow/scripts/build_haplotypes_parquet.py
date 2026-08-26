#!/usr/bin/env python3
"""Backfill tables/haplotypes.parquet from an existing per-peak fastas/ dir.

For runs whose Stage 1 predates the two-stage split (so `fastas/*.fa` + `*.tsv`
are on disk but no consolidated parquet exists). Reads every peak's .fa (seq) +
.tsv (actuation/coverage/asm coords) and writes the one parquet Stage 2 reads.

    pixi run python workflow/scripts/build_haplotypes_parquet.py \
        --fastas-dir results/<run>/fastas \
        --out results/<run>/tables/haplotypes.parquet
"""
import argparse
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastas-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    rows: dict[str, list] = {c: [] for c in COLS}
    n_peaks = 0
    for fa in sorted(args.fastas_dir.glob("*.fa")):
        tsv = args.fastas_dir / f"{fa.stem}.tsv"
        if not tsv.exists():
            continue
        n_peaks += 1
        seqs = read_fa(fa)
        lines = tsv.read_text().splitlines()
        ix = {h: i for i, h in enumerate(lines[0].split("\t"))}
        for ln in lines[1:]:
            f = ln.split("\t")
            s = f[ix["sample"]]
            for k, v in zip(COLS, (fa.stem, s, seqs.get(s, ""), f[ix["actuation"]],
                                   f[ix["coverage"]], f[ix["fire_coverage"]],
                                   f[ix["asm_chr"]], f[ix["asm_start"]],
                                   f[ix["asm_end"]])):
                rows[k].append(v)

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
    print(f"wrote {len(rows['peak_id'])} haplotype rows over {n_peaks} peaks -> {args.out}")


if __name__ == "__main__":
    main()
