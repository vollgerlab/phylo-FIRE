"""Stage 2 -- single-pass analysis on the persisted Stage-1 outputs.

Run AFTER workflow/Snakefile (Stage 1) has produced results/{run}/tables/
{haplotypes,alignments,trees,phylo_mm}.parquet. Each rule is ONE single-pass job
over the whole run (the build_plasmid_table.py pattern), so the DAG is a handful
of nodes and builds in seconds -- no per-peak fan-out, no tree recompute. Reads
only those tables + the AGC, so it also runs locally after downloading
results/{run}/tables/ (~40 MB) and fastas/peak_ids.txt.

    snakemake -s workflow/analysis.smk --config run_name=<run> --rerun-triggers mtime -j <N>
"""

from snakemake.utils import min_version

min_version("8.0")


configfile: "config/config.yaml"


include: "rules/common.smk"


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
            RESULTS + "/tables/variants.hg38.bed.gz",
            RESULTS + "/tables/variants.chm13.bed.gz",
            RESULTS + "/tables/elements.hg38.bed.gz",
            RESULTS + "/tables/drivers.unlifted.tsv",
        ]
        if config.get("dbsnp"):
            out.append(RESULTS + "/tables/variants.hg38.rsid.vcf.gz")
    if PLASMID:
        out += [
            RESULTS + "/plasmid_fiberseq/plasmid_constructs.tsv",
            RESULTS + "/plasmid_fiberseq/plasmid_constructs.fa",
        ]
    if config.get("significant_loci") and config.get("dbsnp") and _hg38_ok():
        out += [
            RESULTS + "/phewas/significant_driver_rsids.tsv",
            RESULTS + "/phewas/significant_driver_rsids.txt",
        ]
    return out


rule all:
    input:
        _targets,


rule re_variants_all:
    """Every peak's re_variants (MSAs + phylo_mm from the Stage-1 summary tables,
actuation from haplotypes.parquet) -> unified parquet + tsv.gz + the plot manifests."""
    input:
        haplotypes=HAP,
        alignments=ALIGNMENTS,
        phylo_mm=PHYLO_MM,
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
        tables=TABLES,
        script=os.path.join(workflow.basedir, "scripts/re_variants_all.py"),
        manifest_dir=RESULTS + "/tables",
        random_n=config.get("random_n", 100),
        seed=config.get("random_seed", 42),
    shell:
        "POLARS_MAX_THREADS={threads} python {params.script} "
        "--haplotypes {input.haplotypes} --tables-dir {params.tables} "
        "--peaks {input.peaks} "
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
        alignments=ALIGNMENTS,
        phylo_mm=PHYLO_MM,
        peaks=PEAKS,
        coords=RESULTS + "/tables/peak_coords.tsv",
    output:
        variants=RESULTS + "/tables/variants.hg38.tsv.gz",
        vcf=RESULTS + "/tables/variants.hg38.vcf.gz",
        tbi=RESULTS + "/tables/variants.hg38.vcf.gz.tbi",
        variants_bed=RESULTS + "/tables/variants.hg38.bed.gz",
        variants_bed_tbi=RESULTS + "/tables/variants.hg38.bed.gz.tbi",
        variants_chm13_bed=RESULTS + "/tables/variants.chm13.bed.gz",
        variants_chm13_bed_tbi=RESULTS + "/tables/variants.chm13.bed.gz.tbi",
        elements_bed=RESULTS + "/tables/elements.hg38.bed.gz",
        elements_bed_tbi=RESULTS + "/tables/elements.hg38.bed.gz.tbi",
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
        tables=TABLES,
        script=os.path.join(workflow.basedir, "scripts/project_variants_all.py"),
        outdir=RESULTS + "/tables",
        agc=config.get("agc", ""),
        flank=config.get("flank", 0),
    shell:
        "python {params.script} --haplotypes {input.haplotypes} "
        "--tables-dir {params.tables} "
        "--coords {input.coords} "
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



rule significant_driver_rsids:
    """Map significant plasmid Fiber-seq loci to their driver variants' rs IDs
    (the PheWAS lookup list). Opt-in: set config['significant_loci'] to the
    plasmid significance table (locus, driver_id, beta, ...); needs config['dbsnp']."""
    input:
        sig=config.get("significant_loci") or "significant_loci_not_set",
        tsv=RESULTS + "/tables/variants.hg38.tsv.gz",
        vcf=RESULTS + "/tables/variants.hg38.rsid.vcf.gz",
        unlifted=RESULTS + "/tables/drivers.unlifted.tsv",
    output:
        tsv=RESULTS + "/phewas/significant_driver_rsids.tsv",
        txt=RESULTS + "/phewas/significant_driver_rsids.txt",
    params:
        script=os.path.join(workflow.basedir, "scripts/significant_driver_rsids.py"),
    shell:
        "python {params.script} --sig {input.sig} --variants {input.tsv} "
        "--rsid-vcf {input.vcf} --unlifted {input.unlifted} "
        "--out-tsv {output.tsv} --out-rsids {output.txt}"

rule plasmid_all:
    """Plasmid construct panel for every single-driver peak (single pass)."""
    input:
        haplotypes=HAP,
        alignments=ALIGNMENTS,
        phylo_mm=PHYLO_MM,
        drivers=RESULTS + "/tables/plot_drivers.txt",
    output:
        constructs=RESULTS + "/plasmid_fiberseq/plasmid_constructs.tsv",
        fa=RESULTS + "/plasmid_fiberseq/plasmid_constructs.fa",
    resources:
        mem_mb=16384,
        runtime=120,
    params:
        tables=TABLES,
        script=os.path.join(workflow.basedir, "scripts/build_plasmid_table.py"),
    shell:
        "python {params.script} --tables-dir {params.tables} "
        "--haplotypes {input.haplotypes} "
        "--peaks {input.drivers} --out-constructs {output.constructs} "
        "--out-fa {output.fa}"


rule plot_all:
    """Render the categorized plots (drivers/heritable/random) named in the
manifests -- one job over the subset, not per-peak."""
    input:
        haplotypes=HAP,
        alignments=ALIGNMENTS,
        phylo_mm=PHYLO_MM,
        trees=TREES,
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
        tables=TABLES,
        script=os.path.join(workflow.basedir, "scripts/plot_all.py"),
        manifest_dir=RESULTS + "/tables",
        outdir=RESULTS + "/plots",
        flank=config.get("flank", 0),
        pclai_arg=_pclai_arg,
        plasmid_arg=("--plasmid" if PLASMID else ""),
        html_arg=("--html" if PLOT_HTML else ""),
    shell:
        'LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}" '
        "python {params.script} --haplotypes {input.haplotypes} "
        "--tables-dir {params.tables} "
        "--manifest-dir {params.manifest_dir} "
        "--outdir {params.outdir} --flank {params.flank} --threads {threads} "
        "{params.pclai_arg} {params.plasmid_arg} {params.html_arg}"
