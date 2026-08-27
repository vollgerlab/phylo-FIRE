"""Stage 2 -- single-pass analysis on the persisted Stage-1 outputs.

Run AFTER workflow/Snakefile (Stage 1) has produced results/{run}/aligned,
trees, phylo_mm and tables/haplotypes.parquet. Each rule is ONE single-pass job
over the whole run (the build_plasmid_table.py pattern), so the DAG is a handful
of nodes and builds in seconds -- no per-peak fan-out, no tree recompute. Reads
only Stage-1 files + the AGC, so it also runs locally after downloading
results/{run}/{aligned,trees,phylo_mm} + tables/haplotypes.parquet.

    snakemake -s workflow/analysis.smk --config run_name=<run> --rerun-triggers mtime -j <N>
"""

import os
import sys
from pathlib import Path
from snakemake.utils import min_version

min_version("8.0")


configfile: "config/config.yaml"


RUN = config.get("run_name", "default")
RESULTS = f"results/{RUN}"

PLOT_HTML = str(config.get("plot_html", "false")).lower() in ("true", "1", "yes")
PROJECT_HG38 = str(config.get("project_hg38", "true")).lower() in ("true", "1", "yes")
PLASMID = str(config.get("plasmid", "false")).lower() in ("true", "1", "yes")

HAP = RESULTS + "/tables/haplotypes.parquet"
PEAKS = RESULTS + "/fastas/peak_ids.txt"
ALIGNED = RESULTS + "/aligned"
TREES = RESULTS + "/trees"
PHYLO_MM = RESULTS + "/phylo_mm"


def _hg38_ok() -> bool:
    """hg38 projection runnable? (union-BED mode + an AGC on disk)."""
    agc = config.get("agc")
    return bool(
        PROJECT_HG38
        and not config.get("preextracted_dir")
        and agc
        and Path(agc).exists()
    )


def _targets(_wildcards):
    out = [
        RESULTS + "/tables/re_variants.parquet",
        RESULTS + "/tables/re_variants.tsv.gz",
        RESULTS + "/plots/.done",
    ]
    if _hg38_ok():
        out += [
            RESULTS + "/tables/variants.hg38.tsv.gz",
            RESULTS + "/tables/variants.hg38.vcf.gz",
            RESULTS + "/tables/variants.hg38.bed",
            RESULTS + "/tables/elements.hg38.bed",
            RESULTS + "/tables/drivers.unlifted.tsv",
        ]
        if config.get("dbsnp"):
            out.append(RESULTS + "/tables/variants.hg38.rsid.vcf.gz")
    if PLASMID:
        out += [
            RESULTS + "/plasmid_fiberseq/plasmid_constructs.tsv",
            RESULTS + "/plasmid_fiberseq/plasmid_constructs.fa",
        ]
    return out


rule all:
    input:
        _targets,


rule re_variants_all:
    """Every peak's re_variants (seqs from aligned/, actuation from the parquet,
drivers from phylo_mm) -> unified parquet + tsv.gz + the plot manifests."""
    input:
        haplotypes=HAP,
        peaks=PEAKS,
    output:
        parquet=RESULTS + "/tables/re_variants.parquet",
        tsv=RESULTS + "/tables/re_variants.tsv.gz",
        drivers=RESULTS + "/tables/plot_drivers.txt",
        heritable=RESULTS + "/tables/plot_heritable_no_driver.txt",
        random=RESULTS + "/tables/plot_random.txt",
    threads: config.get("stage2_threads", 16)
    resources:
        mem_mb=32768,
        runtime=240,
    params:
        script=os.path.join(workflow.basedir, "scripts/re_variants_all.py"),
        aligned=ALIGNED,
        phylo_mm=PHYLO_MM,
        manifest_dir=RESULTS + "/tables",
        random_n=config.get("random_n", 100),
        seed=config.get("random_seed", 42),
    shell:
        "POLARS_MAX_THREADS={threads} python {params.script} "
        "--haplotypes {input.haplotypes} --aligned-dir {params.aligned} "
        "--phylo-mm-dir {params.phylo_mm} --peaks {input.peaks} "
        "--out-parquet {output.parquet} --out-tsv {output.tsv} "
        "--manifest-dir {params.manifest_dir} --random-n {params.random_n} "
        "--seed {params.seed} --threads {threads}"


rule peak_coords:
    """Peak -> CHM13 + hg38 coordinate table, straight from the union BED."""
    input:
        src=config.get("input_table", ""),
    output:
        RESULTS + "/tables/peak_coords.tsv",
    threads: config.get("bgzip_threads", 8)
    resources:
        mem_mb=4096,
        runtime=120,
    params:
        script=os.path.join(workflow.basedir, "scripts/build_peak_coords.py"),
        peak_ids=config.get("peak_ids", ""),
        peak_ids_arg=lambda wc: (
            f"--peak-ids {config['peak_ids']}" if config.get("peak_ids") else ""
        ),
    shell:
        "bgzip -cd -@ {threads} {input.src} | "
        "python {params.script} --out {output} {params.peak_ids_arg}"


rule project_all:
    """Project every variable site to CHM13 + hg38 + assembly (single pass)."""
    input:
        haplotypes=HAP,
        peaks=PEAKS,
        coords=RESULTS + "/tables/peak_coords.tsv",
    output:
        variants=RESULTS + "/tables/variants.hg38.tsv.gz",
        vcf=RESULTS + "/tables/variants.hg38.vcf.gz",
        tbi=RESULTS + "/tables/variants.hg38.vcf.gz.tbi",
        variants_bed=RESULTS + "/tables/variants.hg38.bed",
        elements_bed=RESULTS + "/tables/elements.hg38.bed",
        unlifted=RESULTS + "/tables/drivers.unlifted.tsv",
    # OWN knob (not stage2_threads): this rule's parallelism = concurrent `agc`
    # getctg processes, which abort (SIGABRT) under heavy contention on the shared
    # 3.3 GB archive. Keep it modest (agc-safe); raise with --set-threads
    # project_all=N only if agc tolerates it. Plots/re_variants keep stage2_threads.
    threads: config.get("project_threads", 32)
    resources:
        mem_mb=32768,
        runtime=480,
    params:
        script=os.path.join(workflow.basedir, "scripts/project_variants_all.py"),
        aligned=ALIGNED,
        phylo_mm=PHYLO_MM,
        outdir=RESULTS + "/tables",
        agc=config.get("agc", ""),
        flank=config.get("flank", 0),
    shell:
        "python {params.script} --haplotypes {input.haplotypes} "
        "--aligned-dir {params.aligned} "
        "--phylo-mm-dir {params.phylo_mm} --coords {input.coords} "
        "--peaks {input.peaks} --agc {params.agc} --flank {params.flank} "
        "--outdir {params.outdir} --threads {threads}"


rule annotate_rsids:
    """Add dbSNP rs IDs to the hg38 VCF (bcftools annotate -c ID; allele-aware).
Opt-in: set config['dbsnp'] to a tabix-indexed dbSNP VCF."""
    input:
        vcf=RESULTS + "/tables/variants.hg38.vcf.gz",
    output:
        vcf=RESULTS + "/tables/variants.hg38.rsid.vcf.gz",
    threads: config.get("bgzip_threads", 8)
    resources:
        mem_mb=8192,
        runtime=120,
    params:
        script=os.path.join(workflow.basedir, "scripts/annotate_rsids.py"),
        dbsnp=config.get("dbsnp", ""),
    shell:
        "python {params.script} --vcf {input.vcf} --dbsnp {params.dbsnp} "
        "--out-vcf {output.vcf} --threads {threads}"


rule plasmid_all:
    """Plasmid construct panel for every single-driver peak (single pass)."""
    input:
        haplotypes=HAP,
        drivers=RESULTS + "/tables/plot_drivers.txt",
    output:
        constructs=RESULTS + "/plasmid_fiberseq/plasmid_constructs.tsv",
        fa=RESULTS + "/plasmid_fiberseq/plasmid_constructs.fa",
    resources:
        mem_mb=16384,
        runtime=120,
    params:
        script=os.path.join(workflow.basedir, "scripts/build_plasmid_table.py"),
        aligned=ALIGNED,
        phylo_mm=PHYLO_MM,
    shell:
        "python {params.script} --aligned-dir {params.aligned} "
        "--haplotypes {input.haplotypes} --phylo-mm-dir {params.phylo_mm} "
        "--peaks {input.drivers} --out-constructs {output.constructs} "
        "--out-fa {output.fa}"


def _pclai_arg(_wildcards):
    p = config.get("pclai_bed")
    return f"--pclai-bed {p}" if p and Path(p).exists() else ""


rule plot_all:
    """Render the categorized plots (drivers/heritable/random) named in the
manifests -- one job over the subset, not per-peak."""
    input:
        haplotypes=HAP,
        drivers=RESULTS + "/tables/plot_drivers.txt",
        heritable=RESULTS + "/tables/plot_heritable_no_driver.txt",
        random=RESULTS + "/tables/plot_random.txt",
    output:
        done=touch(RESULTS + "/plots/.done"),
    threads: config.get("stage2_threads", 16)
    resources:
        mem_mb=16384,
        runtime=480,
    params:
        script=os.path.join(workflow.basedir, "scripts/plot_all.py"),
        aligned=ALIGNED,
        trees=TREES,
        phylo_mm=PHYLO_MM,
        manifest_dir=RESULTS + "/tables",
        outdir=RESULTS + "/plots",
        flank=config.get("flank", 0),
        pclai_arg=_pclai_arg,
        plasmid_arg=("--plasmid" if PLASMID else ""),
        html_arg=("--html" if PLOT_HTML else ""),
    shell:
        'LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}" '
        "python {params.script} --haplotypes {input.haplotypes} "
        "--aligned-dir {params.aligned} --trees-dir {params.trees} "
        "--phylo-mm-dir {params.phylo_mm} --manifest-dir {params.manifest_dir} "
        "--outdir {params.outdir} --flank {params.flank} --threads {threads} "
        "{params.pclai_arg} {params.plasmid_arg} {params.html_arg}"
