#!/usr/bin/env python3
"""Rebuild aligned/{peak}.fa for a whole run in ONE pass -- no snakemake DAG.

For a pre-split run whose aligned/ (mafft output) was temp() and got deleted.
Runs the SAME add_flanks + mafft steps the pipeline uses, with the SAME flags, so
in the same pixi env the rebuilt MSA is byte-identical to the original (mafft is
deterministic; verified). It reads the existing fastas/*.fa directly -- never the
extract_fastas checkpoint -- so it reuses the EXACT sequences the trees were built
from, and it loops peaks internally with a thread pool + progress, avoiding the
~46K-node DAG snakemake spends >10 min building for a 23K-peak run.

Resumable + safe to re-run: existing aligned/{peak}.fa are skipped, and each file
is written atomically (temp + rename) so a killed run leaves no partial output.

    pixi run python workflow/scripts/rebuild_aligned.py \
        --fastas-dir results/<run>/fastas --out-dir results/<run>/aligned \
        --agc resources/HPRC_r2_assemblies_0.6.1.agc --flank 250 --threads 16
"""
import argparse
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
ADD_FLANKS = HERE / "add_flanks.py"
# The exact rule-mafft invocation (workflow/Snakefile rule mafft).
MAFFT = ["mafft", "--maxiterate", "1000", "--reorder", "--auto", "--adjustdirection"]


def _strip_r(text: str) -> str:
    """Drop mafft's `_R_` reversed-prefix from headers, matching the pipeline."""
    out = []
    for l in text.splitlines():
        out.append(">" + l[4:] if l.startswith(">_R_") else l)
    return "\n".join(out) + "\n"


def _count_records(path: Path) -> int:
    """Number of `>` headers in a FASTA (0 if unreadable)."""
    try:
        return sum(l.startswith(">") for l in path.read_text().splitlines())
    except OSError:
        return -1


def _valid_aligned(out: Path, fa: Path) -> bool:
    """A skippable aligned/{peak}.fa must be a real alignment, not just present:
    >=2 records, all sequences the same non-zero length, and one record per
    haplotype in the RE-only input (so a truncated/partial file is rebuilt)."""
    try:
        lines = out.read_text().splitlines()
    except OSError:
        return False
    seqs, cur = [], None
    for l in lines:
        if l.startswith(">"):
            seqs.append([])
            cur = seqs[-1]
        elif cur is not None:
            cur.append(l.strip())
    joined = ["".join(s) for s in seqs]
    if len(joined) < 2:
        return False
    lens = {len(s) for s in joined}
    if lens == {0} or len(lens) != 1:
        return False
    return len(joined) == _count_records(fa)


def do_peak(peak, fastas_dir, out_dir, agc, flank, tool_threads):
    out = out_dir / f"{peak}.fa"
    fa, tsv = fastas_dir / f"{peak}.fa", fastas_dir / f"{peak}.tsv"
    if not (fa.exists() and tsv.exists()):
        return peak, "missing"
    if out.exists():
        if _valid_aligned(out, fa):
            return peak, "skip"
        print(f"[rebuild_aligned] {peak}: existing aligned invalid, rebuilding",
              file=sys.stderr)
    with tempfile.NamedTemporaryFile("w", suffix=".flk.fa", delete=False) as fh:
        flk = fh.name
    try:
        subprocess.run(
            ["python", str(ADD_FLANKS), "--fa", str(fa), "--tsv", str(tsv),
             "--agc", agc, "--flank", str(flank), "--threads", str(tool_threads),
             "--out", flk],
            check=True, capture_output=True, text=True,
        )
        res = subprocess.run(
            MAFFT + ["--thread", str(tool_threads), flk],
            check=True, capture_output=True, text=True,
        )
    except subprocess.CalledProcessError as e:
        return peak, f"FAIL:{(e.stderr or '').strip()[-200:]}"
    finally:
        Path(flk).unlink(missing_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".fa.tmp")
    tmp.write_text(_strip_r(res.stdout))
    tmp.rename(out)                       # atomic publish
    return peak, "ok"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastas-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--agc", required=True)
    ap.add_argument("--flank", type=int, default=250)
    ap.add_argument("--threads", type=int, default=16, help="peaks in parallel")
    ap.add_argument("--tool-threads", type=int, default=2,
                    help="threads per add_flanks/mafft call (output is thread-invariant)")
    ap.add_argument("--peaks", type=Path, default=None,
                    help="manifest; default = <fastas-dir>/peak_ids.txt, else glob *.fa")
    args = ap.parse_args()

    if args.peaks and args.peaks.exists():
        peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    else:
        man = args.fastas_dir / "peak_ids.txt"
        peaks = ([p.strip() for p in man.read_text().splitlines() if p.strip()]
                 if man.exists() else sorted(f.stem for f in args.fastas_dir.glob("*.fa")))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[rebuild_aligned] {len(peaks)} peaks, {args.threads} parallel, "
          f"flank={args.flank} -> {args.out_dir}", file=sys.stderr)
    done = ok = skip = missing = 0
    fails = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futs = [ex.submit(do_peak, p, args.fastas_dir, args.out_dir, args.agc,
                          args.flank, args.tool_threads) for p in peaks]
        for fut in as_completed(futs):
            peak, st = fut.result()
            done += 1
            if st == "ok":
                ok += 1
            elif st == "skip":
                skip += 1
            elif st == "missing":
                missing += 1
            else:
                fails.append((peak, st))
            if done % 500 == 0:
                print(f"[rebuild_aligned]   {done}/{len(peaks)} "
                      f"({ok} built, {skip} skipped, {len(fails)} failed)", file=sys.stderr)

    print(f"[rebuild_aligned] done: {ok} built, {skip} already present, "
          f"{missing} missing fastas, {len(fails)} failed", file=sys.stderr)
    for p, s in fails[:10]:
        print(f"  FAIL {p}: {s}", file=sys.stderr)
    if fails:
        sys.exit(1)


if __name__ == "__main__":
    main()
