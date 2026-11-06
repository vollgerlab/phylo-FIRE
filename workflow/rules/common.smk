"""Config shared by both entry points: Stage 1 (workflow/Snakefile) and
Stage 2 (workflow/analysis.smk). Include after configfile:."""

import os
from pathlib import Path

# Namespace outputs per run (results/{run}, temp/{run}) so runs never clobber each other.
# Write paths as "<results>/..." and "<temp>/..."; Snakemake expands them in
# input/output/log. Plain Python (and params) needs workflow.pathvars.apply().
RUN = config.get("run_name", "default")


pathvars:
    results=f"results/{RUN}",
    temp=f"temp/{RUN}",


def cfg_bool(key: str, default: bool) -> bool:
    """Read a bool config value. --config values arrive as strings; bool("false") is truthy."""
    return str(config.get(key, default)).lower() in ("true", "1", "yes")


PLOT_HTML = cfg_bool("plot_html", False)
PROJECT_HG38 = cfg_bool("project_hg38", True)
PLASMID = cfg_bool("plasmid", False)
FLANKED = int(config.get("flank", 0)) > 0

# Stage-1 outputs that Stage 2 reads: the contract between the two entry points.
# Stage 2 reads only these, never the per-peak aligned/ trees/ phylo_mm/ dirs.
HAP = "<results>/tables/haplotypes.parquet"
PEAKS = "<results>/fastas/peak_ids.txt"
ALIGNMENTS = "<results>/tables/alignments.parquet"
TREES = "<results>/tables/trees.parquet"
PHYLO_MM = "<results>/tables/phylo_mm.parquet"


def _hg38_ok() -> bool:
    """hg38 projection runnable? (union-BED mode + an AGC on disk)."""
    agc = config.get("agc")
    return bool(
        PROJECT_HG38
        and not config.get("preextracted_dir")
        and agc
        and Path(agc).exists()
    )


def _pclai_arg(_wildcards):
    p = config.get("pclai_bed")
    return f"--pclai-bed {p}" if p and Path(p).exists() else ""
