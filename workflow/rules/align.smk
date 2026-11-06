rule add_flanks:
    """Add `config["flank"]` bp of AGC context per haplotype (PanSN asm_chr from the TSV).
Only used when flank > 0; otherwise `_mafft_input` sends the RE-only FASTA to mafft.
"""
    input:
        fa="<results>/fastas/{peak_id}.fa",
        tsv="<results>/fastas/{peak_id}.tsv",
    output:
        fa=temp("<temp>/flanked/{peak_id}.fa"),
    group:
        "peak_chain"
    threads: 4
    resources:
        mem_mb=4096,
        runtime=30,
    params:
        script=script_path("add_flanks.py"),
        flank=config.get("flank", 0),
        agc=config.get("agc", ""),
    shell:
        """
        python {params.script} \
            --fa      {input.fa} \
            --tsv     {input.tsv} \
            --agc     {params.agc} \
            --flank   {params.flank} \
            --threads {threads} \
            --out     {output.fa}
        """


rule mafft:
    """Align sequences for a single peak with MAFFT."""
    input:
        fa=_mafft_input,
    output:
        # Persisted: Stage 2's sequence source (re_variants, plasmid, plots).
        fa="<results>/aligned/{peak_id}.fa",
    group:
        "peak_chain"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=60,
    params:
        maxiterate=config.get("mafft_maxiterate", 1000),
    shell:
        """
        mafft --maxiterate {params.maxiterate} --reorder --auto --adjustdirection \
            --thread {threads} {input.fa} > {output.fa}
        """
