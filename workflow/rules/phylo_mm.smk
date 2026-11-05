rule phylo_mm:
    """MCMCglmm H²ₚ (3 priors) in the TSV header + one row per candidate variant with
forward-selection joint WLS stats (`is_split` = driver). Schema in phylo_mm.R."""
    input:
        tree=RESULTS + "/trees/{peak_id}.rooted.nwk",
        actuation=RESULTS + "/fastas/{peak_id}.tsv",
        alignment=RESULTS + "/aligned/{peak_id}.fa",
        re_tree=_re_tree,
        re_alignment=_re_alignment,
    output:
        tsv=RESULTS + "/phylo_mm/{peak_id}.tsv",
    group:
        "peak_chain"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=90,
    params:
        script=script_path("phylo_mm.R"),
        min_cov=config.get("min_coverage", 4),
        flank=config.get("flank", 0),
        h2_width_max=config.get("h2_hpd_width_max", 0.5),
        h2_lower_min=config.get("h2_hpd_lower_min", 0.1),
    shell:
        """
        OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        Rscript {params.script} \
            --tree      {input.tree} \
            --actuation {input.actuation} \
            --alignment {input.alignment} \
            --re-tree      {input.re_tree} \
            --re-alignment {input.re_alignment} \
            --tsv       {output.tsv} \
            --peak-id   {wildcards.peak_id} \
            --min-cov   {params.min_cov} \
            --flank     {params.flank} \
            --h2-width-max {params.h2_width_max} \
            --h2-lower-min {params.h2_lower_min}
        """


rule stage1_tables:
    """Combine every peak's MSA, rooted tree and phylo_mm TSV into three parquet
tables (single pass), so Stage 2 can run from ``tables/`` alone, e.g. locally
after downloading it. See workflow/scripts/stage1_tables.py."""
    input:
        aligned=expand(RESULTS + "/aligned/{peak_id}.fa", peak_id=PEAK_IDS),
        trees=expand(RESULTS + "/trees/{peak_id}.rooted.nwk", peak_id=PEAK_IDS),
        phylo_mm=expand(RESULTS + "/phylo_mm/{peak_id}.tsv", peak_id=PEAK_IDS),
        peaks=PEAKS,
    output:
        ALIGNMENTS,
        TREES,
        PHYLO_MM,
    params:
        script=os.path.join(workflow.basedir, "scripts/stage1_tables.py"),
    shell:
        "python {params.script} --results {RESULTS} --peaks {input.peaks}"
