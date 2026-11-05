#!/usr/bin/env python3
"""Stage-1 summary tables: the per-peak MSAs, rooted trees and phylo_mm stats
combined into three parquet files, so Stage 2 can run from a handful of files
instead of ~70K (e.g. locally, after downloading ``results/{run}/tables``).

  * alignments.parquet  peak_id, idx, name, seq   (one row per aligned haplotype;
                        ``idx`` keeps the FASTA order, ``seq`` keeps the gaps)
  * trees.parquet       peak_id, newick           (the ancestral-proxy-rooted tree)
  * phylo_mm.parquet    peak_id, text             (the full phylo_mm TSV, header
                        comments included, so it is lossless)

Build (end of Stage 1)::

    python stage1_tables.py --results results/<run> --peaks <peak_ids.txt>

Read (Stage 2): ``Stage1Store(tables_dir)`` hands back a per-peak file path
materialized from the parquet tables into a temp dir. Stage 2 reads only the
tables, never the per-peak directories.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

TABLES = {
    "aligned": "alignments.parquet",
    "trees": "trees.parquet",
    "phylo_mm": "phylo_mm.parquet",
}


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    """(name, seq) records in file order; name is the header up to whitespace."""
    recs, name, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(chunks)))
            name, chunks = line[1:].split()[0], []
        elif line:
            chunks.append(line.strip())
    if name is not None:
        recs.append((name, "".join(chunks)))
    return recs


def _read_peak(results: Path, peak: str) -> tuple[list, str | None, str | None]:
    """One peak's (aligned records, rooted newick, phylo_mm text); None if absent."""
    fa = results / "aligned" / f"{peak}.fa"
    nwk = results / "trees" / f"{peak}.rooted.nwk"
    tsv = results / "phylo_mm" / f"{peak}.tsv"
    return (
        _read_fasta(fa) if fa.exists() else [],
        nwk.read_text().strip() if nwk.exists() else None,
        tsv.read_text() if tsv.exists() else None,
    )


def build(results: Path, peaks: list[str], threads: int = 16) -> None:
    """Write the three Stage-1 summary tables to ``results/tables/``.

    Files are read in parallel: on a network filesystem (Lustre) the cost is
    per-file open latency, which threads hide. Peak order is preserved."""
    aln_rows: dict[str, list] = {"peak_id": [], "idx": [], "name": [], "seq": []}
    trees: dict[str, list] = {"peak_id": [], "newick": []}
    pmm: dict[str, list] = {"peak_id": [], "text": []}
    step = max(1, len(peaks) // 20)
    with ThreadPoolExecutor(max_workers=threads) as ex:
        results_iter = ex.map(lambda pk: _read_peak(results, pk), peaks)
        for n, (peak, (recs, newick, text)) in enumerate(zip(peaks, results_iter), 1):
            for i, (name, seq) in enumerate(recs):
                aln_rows["peak_id"].append(peak)
                aln_rows["idx"].append(i)
                aln_rows["name"].append(name)
                aln_rows["seq"].append(seq)
            if newick is not None:
                trees["peak_id"].append(peak)
                trees["newick"].append(newick)
            if text is not None:
                pmm["peak_id"].append(peak)
                pmm["text"].append(text)
            if n % step == 0 or n == len(peaks):
                print(f"[stage1_tables] read {n}/{len(peaks)} peaks", file=sys.stderr, flush=True)
    out = results / "tables"
    out.mkdir(parents=True, exist_ok=True)
    for key, data in (("aligned", aln_rows), ("trees", trees), ("phylo_mm", pmm)):
        pl.DataFrame(data).write_parquet(out / TABLES[key], compression="zstd")
    print(f"[stage1_tables] {len(peaks)} peaks -> {len(aln_rows['peak_id'])} aligned "
          f"haplotypes, {len(trees['peak_id'])} trees, {len(pmm['peak_id'])} phylo_mm "
          f"in {out}")


class Stage1Store:
    """Per-peak Stage-1 files served from the summary tables.

    ``aligned(peak)``, ``tree(peak)`` and ``phylo_mm(peak)`` return a Path to the
    peak's file (written to a temp dir on first use), or None if the peak is
    not in that table."""

    def __init__(self, tables_dir: Path):
        self.tables_dir = Path(tables_dir)
        self._tables: dict[str, dict[str, object]] = {}
        self._tmp = Path(tempfile.mkdtemp(prefix="stage1_store_"))
        self._lock = threading.Lock()

    def _table(self, key: str) -> dict[str, object]:
        """peak_id -> content for one summary table (loaded once, lazily)."""
        with self._lock:
            if key not in self._tables:
                path = self.tables_dir / TABLES[key]
                if key == "aligned":
                    df = pl.read_parquet(path).sort(["peak_id", "idx"])
                    self._tables[key] = {
                        p[0]: list(zip(g["name"].to_list(), g["seq"].to_list()))
                        for p, g in df.partition_by("peak_id", as_dict=True).items()
                    }
                else:
                    col = "newick" if key == "trees" else "text"
                    df = pl.read_parquet(path)
                    self._tables[key] = dict(zip(df["peak_id"].to_list(), df[col].to_list()))
            return self._tables[key]

    def _get(self, key: str, peak: str, name: str) -> Path | None:
        content = self._table(key).get(peak)
        if content is None:
            return None
        out = self._tmp / key / name
        if not out.exists():
            out.parent.mkdir(parents=True, exist_ok=True)
            if key == "aligned":
                content = "".join(f">{n}\n{s}\n" for n, s in content)
            elif key == "trees":
                content = content + "\n"
            out.write_text(content)
        return out

    def aligned(self, peak: str) -> Path | None:
        return self._get("aligned", peak, f"{peak}.fa")

    def tree(self, peak: str) -> Path | None:
        return self._get("trees", peak, f"{peak}.rooted.nwk")

    def phylo_mm(self, peak: str) -> Path | None:
        return self._get("phylo_mm", peak, f"{peak}.tsv")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", required=True, type=Path, help="results/<run>")
    ap.add_argument("--peaks", required=True, type=Path, help="peak id list")
    ap.add_argument("--threads", type=int, default=16, help="parallel file readers")
    args = ap.parse_args()
    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    build(args.results, peaks, args.threads)


if __name__ == "__main__":
    main()
