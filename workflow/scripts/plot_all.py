#!/usr/bin/env python3
"""Single-pass plotting for a WHOLE run (Stage 2, no snakemake DAG).

Renders the categorized plots (drivers / heritable_no_driver / random) listed in
the aggregate manifests. For each selected peak it materializes a temp actuation
TSV from haplotypes.parquet and calls the existing per-peak plot script
(unchanged) -- so only the ~few-K plotted subset is rendered, in a worker pool.
For driver peaks it also builds the plasmid construct panel (via
build_plasmid_constructs) so the constructs draw under the MSA.

    pixi run python workflow/scripts/plot_all.py \
        --haplotypes results/<run>/tables/haplotypes.parquet \
        --aligned-dir results/<run>/aligned --trees-dir results/<run>/trees \
        --phylo-mm-dir results/<run>/phylo_mm --manifest-dir results/<run>/tables \
        --outdir results/<run>/plots --flank <F> [--pclai-bed BED] [--plasmid] [--html]
"""
import argparse
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hap_io import load_parts, peak_actuation  # noqa: E402
from build_plasmid_constructs import build_peak as plasmid_peak  # noqa: E402
from stage1_tables import Stage1Store  # noqa: E402

CATEGORIES = ["drivers", "heritable_no_driver", "random"]
PLOT_SCRIPT = Path(__file__).resolve().parent / "plot_tree_actuation_html.py"


def _write_actuation(parts, peak, path):
    act = peak_actuation(parts, peak)
    with open(path, "w") as fh:
        fh.write("sample\tactuation\tcoverage\n")
        for s, (a, c) in act.items():
            # coverage is an integer read count; the parquet carries it as f64,
            # but the plot's TSV reader does int(cov) -> emit it integral.
            fh.write(f"{s}\t{a}\t{int(round(c))}\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--haplotypes", required=True, type=Path)
    ap.add_argument("--aligned-dir", required=True, type=Path)
    ap.add_argument("--trees-dir", required=True, type=Path)
    ap.add_argument("--phylo-mm-dir", required=True, type=Path)
    ap.add_argument("--manifest-dir", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--pclai-bed", type=Path, default=None)
    ap.add_argument("--plasmid", action="store_true")
    ap.add_argument("--html", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    parts = load_parts(args.haplotypes)
    tmp = Path(tempfile.mkdtemp(prefix="plot_all_"))
    env = dict(os.environ)
    cp = env.get("CONDA_PREFIX", "")
    env["LD_LIBRARY_PATH"] = f"{cp}/lib:" + env.get("LD_LIBRARY_PATH", "")

    jobs = []
    for cat in CATEGORIES:
        manifest = args.manifest_dir / f"plot_{cat}.txt"
        if not manifest.exists():
            continue
        for peak in (p.strip() for p in manifest.read_text().splitlines() if p.strip()):
            jobs.append((cat, peak))

    failed: list[str] = []
    store = Stage1Store(args.aligned_dir.parent)  # per-peak dirs, else tables/*.parquet

    def render(job) -> None:
        cat, peak = job
        aligned = store.aligned(peak)
        tree = store.tree(peak)
        pmm = store.phylo_mm(peak)
        if aligned is None or tree is None or pmm is None:
            failed.append(peak)
            return
        act_tsv = tmp / f"{peak}.act.tsv"
        _write_actuation(parts, peak, act_tsv)
        outdir = args.outdir / cat
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "python", str(PLOT_SCRIPT),
            "--tree", str(tree), "--actuation", str(act_tsv),
            "--alignment", str(aligned), "--pdf", str(outdir / f"{peak}.pdf"),
            "--peak-id", peak, "--phylo-mm", str(pmm), "--flank", str(args.flank),
        ]
        if args.html:
            cmd += ["--html", str(outdir / f"{peak}.html")]
        if args.pclai_bed and args.pclai_bed.exists():
            cmd += ["--pclai-bed", str(args.pclai_bed)]
        if args.plasmid and cat == "drivers":
            _, panel = plasmid_peak(aligned, peak_actuation(parts, peak), pmm, peak)
            if panel:
                panel_fa = tmp / f"{peak}.panel.fa"
                with panel_fa.open("w") as fh:
                    for name, anchor, label, seq in panel:
                        fh.write(f">{name} anchor_leaf={anchor} type=plasmid "
                                 f"label={label}\n{seq}\n")
                cmd += ["--panel", str(panel_fa)]
        r = subprocess.run(cmd, env=env, capture_output=True, text=True)
        if r.returncode != 0:
            failed.append(peak)
            print(f"[plot_all] {peak} FAILED: {r.stderr.strip()[-300:]}", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(render, jobs))
    print(f"[plot_all] rendered {len(jobs) - len(failed)}/{len(jobs)} plots"
          + (f"; {len(failed)} failed" if failed else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
