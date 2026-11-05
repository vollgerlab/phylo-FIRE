rule extract_fastas:
    """Stream the union BED into per-peak FASTA/TSV, haplotypes.parquet and the survivor
manifest, or symlink a preextracted_dir (e.g. PARE paralog groups). Not a checkpoint.
"""
    input:
        src=_extract_fastas_input,
        peak_ids=optional_file("peak_ids"),
        highlight=optional_file("highlight_file"),
    output:
        # Every requested peak is declared so the DAG is static. A peak that yields
        # nothing leaves a declared output missing; extract_fastas.py exits listing
        # those ids. fastas/ is temp; haplotypes.parquet is the persisted copy.
        fa=temp(expand(RESULTS + "/fastas/{peak_id}.fa", peak_id=PEAK_IDS)),
        tsv=temp(expand(RESULTS + "/fastas/{peak_id}.tsv", peak_id=PEAK_IDS)),
        manifest=PEAKS,
        haplotypes=HAP,
    threads: config.get("bgzip_threads", 16)
    resources:
        mem_mb=32768,
        runtime=480,
    params:
        script=script_path("extract_fastas.py"),
        ingest=script_path("ingest_preextracted.py"),
        outdir=RESULTS + "/fastas",
        preextracted=config.get("preextracted_dir") or "",
        min_seqs=config.get("min_seqs", 2),
        opt_args=_extract_opt_args,
    # The union BED is asm-coord-sorted; group rows by consensus_peak_id (col 5)
    # for the streaming extractor. `read` keeps the header first (BSD `sed -u 1q`
    # reads past line 1 on macOS and silently drops rows).
    shell:
        """
        if [ -n "{params.preextracted}" ]; then
            python {params.ingest:q} --link-from {input.src:q} \
                --fastas-dir {params.outdir:q} --out {output.haplotypes:q} \
                --threads {threads}
        else
            bgzip -cd -@ {threads} {input.src:q} \
                | (IFS= read -r hdr; printf '%s\\n' "$hdr"; sort -t$'\\t' -k5,5 -S 50%) \
                | python {params.script:q} --outdir {params.outdir:q} \
                    --min-seqs {params.min_seqs} \
                    --haplotypes-parquet {output.haplotypes:q} {params.opt_args}
        fi
        """
