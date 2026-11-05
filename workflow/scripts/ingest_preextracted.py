#!/usr/bin/env python3
"""Build tables/haplotypes.parquet from per-peak fastas/{peak}.fa + .tsv pairs.

The one haplotypes.parquet builder outside extract_fastas.py. Two uses:

- extract_fastas, preextracted_dir mode: --link-from symlinks a pre-extracted
  dir's pairs (e.g. PARE paralog groups) into --fastas-dir, writes its
  peak_ids.txt, then builds the parquet.
- Backfill for a run whose Stage 1 predates the two-stage split (fastas/ on
  disk, no parquet).

Reads are threaded (I/O-latency-bound on Lustre: ~2x46K tiny files).

    python workflow/scripts/ingest_preextracted.py \
        [--link-from <preextracted_dir>] \
        --fastas-dir results/<run>/fastas \
        --out results/<run>/tables/haplotypes.parquet [--threads 32]
"""

import argparse
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

COLS = [
    "peak_id",
    "sample",
    "seq",
    "actuation",
    "coverage",
    "fire_coverage",
    "asm_chr",
    "asm_start",
    "asm_end",
]


def sanitize(peak: str) -> str:
    """Mirror extract_fastas.sanitize: peak_ids become filenames, `#`->`_`."""
    return peak.replace("#", "_")


def link_preextracted(src_dir: Path, dst_dir: Path) -> int:
    """Symlink {peak}.fa/.tsv pairs from src_dir into dst_dir and write peak_ids.txt.

    Returns the number of pairs linked. A .fa without a .tsv is skipped.
    """
    src_dir = src_dir.resolve()
    dst_dir.mkdir(parents=True, exist_ok=True)
    linked: list[str] = []
    for fa in sorted(src_dir.glob("*.fa")):
        tsv = src_dir / f"{fa.stem}.tsv"
        if not tsv.exists():
            continue
        peak = sanitize(fa.stem)
        for src, dst in ((fa, dst_dir / f"{peak}.fa"), (tsv, dst_dir / f"{peak}.tsv")):
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src)
        linked.append(peak)
    (dst_dir / "peak_ids.txt").write_text("\n".join(linked) + "\n")
    return len(linked)


def read_fa(p: Path) -> dict[str, str]:
    """Read a FASTA into {name: seq}, name being the header up to the first space."""
    d: dict[str, list[str]] = {}
    name = None
    for ln in p.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            d[name] = []
        elif name is not None:
            d[name].append(ln.strip())
    return {k: "".join(v) for k, v in d.items()}


def read_peak(fa: Path) -> list[tuple] | None:
    """Return one tuple (in COLS order) per haplotype of a peak, or None if no .tsv."""
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
        out.append(
            (
                fa.stem,
                s,
                seqs.get(s, ""),
                f[ix["actuation"]],
                f[ix["coverage"]],
                f[ix["fire_coverage"]],
                f[ix["asm_chr"]],
                f[ix["asm_start"]],
                f[ix["asm_end"]],
            )
        )
    return out


def build_parquet(fastas_dir: Path, out: Path, threads: int) -> None:
    """Consolidate every fastas_dir/{peak}.fa + .tsv pair into one parquet at out."""
    fas = sorted(fastas_dir.glob("*.fa"))
    print(
        f"[haplotypes] reading {len(fas)} peaks from {fastas_dir} "
        f"with {threads} threads ...",
        file=sys.stderr,
    )
    rows: dict[str, list] = {c: [] for c in COLS}
    done = 0
    with ThreadPoolExecutor(max_workers=threads) as ex:
        for res in ex.map(read_peak, fas):
            done += 1
            if res:
                for tup in res:
                    for k, v in zip(COLS, tup):
                        rows[k].append(v)
            if done % 1000 == 0:
                print(
                    f"[haplotypes]   {done}/{len(fas)} peaks "
                    f"({len(rows['peak_id'])} rows)",
                    file=sys.stderr,
                )

    out.parent.mkdir(parents=True, exist_ok=True)
    (
        pl.DataFrame(rows, schema={c: pl.Utf8 for c in COLS})
        .with_columns(
            pl.col("actuation").cast(pl.Float64, strict=False),
            pl.col("coverage").cast(pl.Float64, strict=False),
            pl.col("fire_coverage").cast(pl.Float64, strict=False),
        )
        .write_parquet(out)
    )
    print(
        f"[haplotypes] wrote {len(rows['peak_id'])} rows over {len(fas)} "
        f"peaks -> {out}",
        file=sys.stderr,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--link-from",
        type=Path,
        help="pre-extracted dir of {peak}.fa/.tsv pairs to symlink into --fastas-dir",
    )
    ap.add_argument("--fastas-dir", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    if args.link_from:
        n = link_preextracted(args.link_from, args.fastas_dir)
        print(
            f"linked {n} pre-extracted FASTA/TSV pairs from {args.link_from}",
            file=sys.stderr,
        )
    build_parquet(args.fastas_dir, args.out, args.threads)


if __name__ == "__main__":
    main()
