#!/usr/bin/env python3
"""Single-pass re_variants for a WHOLE run (Stage 2, no snakemake DAG).

Loops every peak, builds its per-unique-RE-sequence table via
`build_re_variants.build` (sequences from the persisted aligned/ MSA, actuation
from haplotypes.parquet, drivers from phylo_mm), writes a per-peak TSV to a
scratch dir, then reuses `aggregate_re_variants.py` to emit the unified parquet +
tsv.gz + the three plot manifests.

    pixi run python workflow/scripts/re_variants_all.py \
        --haplotypes results/<run>/tables/haplotypes.parquet \
        --aligned-dir results/<run>/aligned \
        --phylo-mm-dir results/<run>/phylo_mm \
        --peaks results/<run>/fastas/peak_ids.txt \
        --out-parquet results/<run>/tables/re_variants.parquet \
        --out-tsv results/<run>/tables/re_variants.tsv.gz \
        --manifest-dir results/<run>/tables --random-n 100 --seed 42 --threads 16
"""
import argparse
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_re_variants import build, COLUMNS  # noqa: E402
from hap_io import load_parts, peak_actuation  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--haplotypes", required=True, type=Path)
    ap.add_argument("--aligned-dir", required=True, type=Path)
    ap.add_argument("--phylo-mm-dir", required=True, type=Path)
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
    tmp = Path(tempfile.mkdtemp(prefix="re_variants_all_"))
    skipped: list[str] = []

    def do_peak(peak: str) -> None:
        aligned = args.aligned_dir / f"{peak}.fa"
        pmm = args.phylo_mm_dir / f"{peak}.tsv"
        if not (aligned.exists() and pmm.exists()):
            skipped.append(peak)
            return
        rows = build(str(aligned), peak_actuation(parts, peak), str(pmm))
        df = (
            pl.DataFrame(rows).select(COLUMNS)
            if rows
            else pl.DataFrame(schema={c: pl.Utf8 for c in COLUMNS})
        )
        df.write_csv(tmp / f"{peak}.tsv", separator="\t")

    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(do_peak, peaks))
    if skipped:
        print(f"[re_variants_all] {len(skipped)} peaks missing aligned/phylo_mm "
              f"(skipped): {', '.join(skipped[:10])}", file=sys.stderr)

    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    agg = Path(__file__).resolve().parent / "aggregate_re_variants.py"
    subprocess.run(
        ["python", str(agg), "--indir", str(tmp), "--peaks", str(args.peaks),
         "--out-parquet", str(args.out_parquet), "--out-tsv", str(args.out_tsv),
         "--manifest-dir", str(args.manifest_dir),
         "--random-n", str(args.random_n), "--seed", str(args.seed)],
        check=True,
    )
    print(f"[re_variants_all] done -> {args.out_parquet}", file=sys.stderr)


if __name__ == "__main__":
    main()
