rule phylo_mm:
    """MCMCglmm H²ₚ (3 priors) in the TSV header + one row per candidate variant with
forward-selection joint WLS stats (`is_split` = driver). Schema in phylo_mm.R."""
    input:
        tree=rules.root_tree.output.rooted,
        # extract_fastas declares the expanded list, so there is no pattern to reference.
        actuation="<results>/fastas/{peak_id}.tsv",
        alignment=rules.mafft.output.fa,
        re_tree=_re_tree,
        re_alignment=_re_alignment,
    output:
        tsv="<results>/phylo_mm/{peak_id}.tsv",
    log:
        "<results>/logs/phylo_mm/{peak_id}.log",
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
            --h2-lower-min {params.h2_lower_min} > {log} 2>&1
        """


rule stage1_tables:
    """Combine every peak's MSA, rooted tree and phylo_mm TSV into three parquet
tables (single pass), so Stage 2 can run from ``tables/`` alone, e.g. locally
after downloading it. See workflow/scripts/stage1_tables.py."""
    input:
        aligned=expand(rules.mafft.output.fa, peak_id=PEAK_IDS),
        trees=expand(rules.root_tree.output.rooted, peak_id=PEAK_IDS),
        phylo_mm=expand(rules.phylo_mm.output.tsv, peak_id=PEAK_IDS),
        peaks=rules.extract_fastas.output.manifest,
    output:
        alignments=ALIGNMENTS,
        trees=TREES,
        phylo_mm=PHYLO_MM,
    log:
        "<results>/logs/stage1_tables.log",
    params:
        script=os.path.join(workflow.basedir, "scripts/stage1_tables.py"),
        # results/{run}: the script reads aligned/ trees/ phylo_mm/ and writes tables/.
        results=lambda wc, output: os.path.dirname(os.path.dirname(output.alignments)),
    shell:
        "python {params.script} --results {params.results} --peaks {input.peaks} > {log} 2>&1"
