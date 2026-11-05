# OpenBLAS/MKL ignore the rule's threads and grab every core; cap them.
BLAS_CAP = "OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} MKL_NUM_THREADS={threads} "


rule encode_indels:
    """Simmons-Ochoterena indel coding so IQ-TREE models gaps (GTR2+ASC) instead of
treating them as missing; coding rules are in scripts/encode_indels.py."""
    input:
        fa=RESULTS + "/aligned/{peak_id}.fa",
    output:
        nex=temp(TEMP + "/coded/{peak_id}/{peak_id}.nex"),
        nuc=temp(TEMP + "/coded/{peak_id}/{peak_id}_nucleic.phy"),
        indel=temp(TEMP + "/coded/{peak_id}/{peak_id}_indel.phy"),
    group:
        "peak_chain"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=10,
    params:
        script=script_path("encode_indels.py"),
        peak_id=lambda wc: wc.peak_id,
        out_dir=lambda wc: f"{TEMP}/coded/{wc.peak_id}",
    shell:
        """
        python {params.script} \
            --aligned {input.fa} \
            --peak-id {params.peak_id} \
            --out-dir {params.out_dir}
        """


rule iqtree:
    """ML tree from the DNA + indel partition. Models are fixed by the .nex charpartition
(GTR+I+G on DNA, GTR2+ASC on indels), so IQ-TREE runs no model selection."""
    input:
        nex=TEMP + "/coded/{peak_id}/{peak_id}.nex",
        nuc=TEMP + "/coded/{peak_id}/{peak_id}_nucleic.phy",
        indel=TEMP + "/coded/{peak_id}/{peak_id}_indel.phy",
    output:
        # unrooted ML tree — temp: superseded by the rooted .nwk (root_tree keeps
        # branch lengths + bootstrap); only root_tree + re_iqtree-fallback read it.
        treefile=temp(TEMP + "/trees/{peak_id}.treefile"),
    group:
        "peak_chain"
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120,
    params:
        prefix=lambda wc: os.path.abspath(f"{TEMP}/trees/{wc.peak_id}"),
        peak_id=lambda wc: wc.peak_id,
        bootstrap=iqtree_bootstrap,
    shell:
        # iqtree runs beside the .nex (relative .phy refs) with an absolute --prefix
        # into TEMP/trees; shadow: conflicts with that prefix, so aux files are deleted
        # after. --redo overwrites partial runs.
        """
        (cd $(dirname {input.nex}) && iqtree -p {params.peak_id}.nex \
            --prefix {params.prefix} -T AUTO --threads-max {threads} {params.bootstrap} --redo)
        # Delete every IQ-TREE aux file sharing the prefix except the treefile,
        # so new aux extensions across IQ-TREE versions don't leak.
        find "$(dirname {params.prefix})" -maxdepth 1 \
            -name "$(basename {params.prefix}).*" \
            ! -name "$(basename {params.prefix}).treefile" -delete
        """


rule root_tree:
    """Root once at the ancestral proxy (haps closest to the MSA consensus; midpoint if
that set is not monophyletic). phylo_mm and the plots never re-root, so clades match.
"""
    input:
        treefile=TEMP + "/trees/{peak_id}.treefile",
        alignment=RESULTS + "/aligned/{peak_id}.fa",
    output:
        rooted=RESULTS + "/trees/{peak_id}.rooted.nwk",
    group:
        "peak_chain"
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    params:
        script=script_path("root_tree.R"),
    shell:
        BLAS_CAP
        + "Rscript {params.script} {input.treefile} {input.alignment} {output.rooted}"


# ── RE-only tree ─────────────────────────────────────────────────────────────
# The MSA trimmed to the RE, through the same encode → iqtree → root chain, so
# H²ₚ(RE) and H²ₚ(region) differ only in scope. Built only when flank > 0.


rule re_trim:
    """Slice the flanked MSA down to the RE columns."""
    input:
        fa=RESULTS + "/aligned/{peak_id}.fa",
    output:
        fa=temp(TEMP + "/re/aligned/{peak_id}.fa"),
    group:
        "peak_chain"
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    params:
        script=script_path("trim_to_re.py"),
        flank=config.get("flank", 0),
    shell:
        "python {params.script} --aln {input.fa} --flank {params.flank} --out {output.fa}"


use rule encode_indels as re_encode_indels with:
    """encode_indels on the RE-only alignment (mirrors rule encode_indels)."""
    input:
        fa=TEMP + "/re/aligned/{peak_id}.fa",
    output:
        nex=temp(TEMP + "/re/coded/{peak_id}/{peak_id}.nex"),
        nuc=temp(TEMP + "/re/coded/{peak_id}/{peak_id}_nucleic.phy"),
        indel=temp(TEMP + "/re/coded/{peak_id}/{peak_id}_indel.phy"),
    params:
        out_dir=lambda wc: f"{TEMP}/re/coded/{wc.peak_id}",


rule re_iqtree:
    """RE-only ML tree (mirrors rule iqtree). A monomorphic RE makes IQ-TREE fail; the main
tree is copied instead, which only needs the right tips since H²ₚ(RE) is then NA."""
    input:
        nex=TEMP + "/re/coded/{peak_id}/{peak_id}.nex",
        nuc=TEMP + "/re/coded/{peak_id}/{peak_id}_nucleic.phy",
        indel=TEMP + "/re/coded/{peak_id}/{peak_id}_indel.phy",
        main_tree=TEMP + "/trees/{peak_id}.treefile",
    output:
        # internal: only re_root_tree → phylo_mm read it, nothing downstream.
        treefile=temp(TEMP + "/re/trees/{peak_id}.treefile"),
    group:
        "peak_chain"
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120,
    params:
        prefix=lambda wc: os.path.abspath(f"{TEMP}/re/trees/{wc.peak_id}"),
        peak_id=lambda wc: wc.peak_id,
        bootstrap=iqtree_bootstrap,
    shell:
        """
        if (cd $(dirname {input.nex}) && iqtree -p {params.peak_id}.nex \
                --prefix {params.prefix} -T AUTO --threads-max {threads} {params.bootstrap} --redo); then
            find "$(dirname {params.prefix})" -maxdepth 1 \
                -name "$(basename {params.prefix}).*" \
                ! -name "$(basename {params.prefix}).treefile" -delete
        else
            echo "[re_iqtree] {params.peak_id}: IQ-TREE failed (monomorphic RE?); using main tree" >&2
            cp {input.main_tree} {output.treefile}
        fi
        """


use rule root_tree as re_root_tree with:
    """Ancestral-proxy root the RE-only tree (mirrors rule root_tree)."""
    input:
        treefile=TEMP + "/re/trees/{peak_id}.treefile",
        alignment=TEMP + "/re/aligned/{peak_id}.fa",
    output:
        # internal: only phylo_mm (H²_re scope) reads this, then it's done.
        rooted=temp(TEMP + "/re/trees/{peak_id}.rooted.nwk"),
