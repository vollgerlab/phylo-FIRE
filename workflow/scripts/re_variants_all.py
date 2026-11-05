#!/usr/bin/env python3
"""Single-pass re_variants for a WHOLE run (Stage 2, no snakemake DAG).

Loops every peak, builds its per-unique-RE-sequence table via
`build_re_variants.build` (MSA + phylo_mm from the Stage-1 summary tables,
actuation from haplotypes.parquet), writes a per-peak TSV to a scratch dir, then
runs `aggregate_re_variants.aggregate` on it to emit the unified parquet +
tsv.gz + the three plot manifests.

    pixi run python workflow/scripts/re_variants_all.py \
        --haplotypes results/<run>/tables/haplotypes.parquet \
        --tables-dir results/<run>/tables \
        --peaks results/<run>/fastas/peak_ids.txt \
        --out-parquet results/<run>/tables/re_variants.parquet \
        --out-tsv results/<run>/tables/re_variants.tsv.gz \
        --manifest-dir results/<run>/tables --random-n 100 --seed 42 --threads 16
"""
import argparse
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from aggregate_re_variants import aggregate  # noqa: E402
from build_re_variants import build, COLUMNS  # noqa: E402
from hap_io import load_parts, peak_actuation  # noqa: E402
from stage1_tables import Stage1Store  # noqa: E402


def write_peak_tables(peaks: list[str], parts: dict, store: Stage1Store,
                      outdir: Path, threads: int) -> list[str]:
    """Write ``{outdir}/{peak}.tsv`` for every peak; returns the peaks skipped
    for lack of an alignment or phylo_mm table."""
    skipped: list[str] = []

    def do_peak(peak: str) -> None:
        aligned = store.aligned(peak)
        pmm = store.phylo_mm(peak)
        if aligned is None or pmm is None:
            skipped.append(peak)
            return
        rows = build(str(aligned), peak_actuation(parts, peak), str(pmm))
        df = (
            pl.DataFrame(rows).select(COLUMNS)
            if rows
            else pl.DataFrame(schema={c: pl.Utf8 for c in COLUMNS})
        )
        df.write_csv(outdir / f"{peak}.tsv", separator="\t")

    with ThreadPoolExecutor(max_workers=threads) as ex:
        list(ex.map(do_peak, peaks))
    return skipped


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--haplotypes", required=True, type=Path)
    ap.add_argument("--tables-dir", required=True, type=Path,
                    help="results/<run>/tables (Stage-1 summary parquets)")
    ap.add_argument("--peaks", required=True, type=Path)
    ap.add_argument("--out-parquet", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--manifest-dir", required=True, type=Path)
    ap.add_argument("--random-n", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    parts = load_parts(args.haplotypes)
    store = Stage1Store(args.tables_dir)
    with tempfile.TemporaryDirectory(prefix="re_variants_all_") as tmpdir:
        tmp = Path(tmpdir)
        skipped = write_peak_tables(peaks, parts, store, tmp, args.threads)
        if skipped:
            print(f"[re_variants_all] {len(skipped)} peaks missing aligned/phylo_mm "
                  f"(skipped): {', '.join(skipped[:10])}", file=sys.stderr)
        aggregate(tmp, peaks, args.out_parquet, args.out_tsv, args.manifest_dir,
                  random_n=args.random_n, seed=args.seed)
    print(f"[re_variants_all] done -> {args.out_parquet}", file=sys.stderr)


if __name__ == "__main__":
    main()
