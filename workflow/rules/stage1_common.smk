"""Stage-1 helpers: wildcard constraints, the parse-time peak set and the input functions."""


# Keep peak_id from swallowing a slash.
wildcard_constraints:
    peak_id=r"[^/]+",


def script_path(name: str):
    """workflow.source_path() for workflow/scripts/<name>. Joining onto basedir keeps
    the path free of `rules/../`, so params and shell commands match the old
    single-file Snakefile."""
    return workflow.source_path(os.path.join(workflow.basedir, "scripts", name))


def count_unique_seqs_phylip(phy_path: str) -> int:
    """Count unique sequences in a one-line-per-tip PHYLIP file."""
    seqs: set[str] = set()
    with open(phy_path) as fh:
        next(fh, None)  # skip the header `<n_taxa> <n_cols>`
        for line in fh:
            parts = line.strip().split(None, 1)
            if len(parts) == 2:
                seqs.add(parts[1])
    return len(seqs)


# ── static peak set (resolved at parse time, no checkpoint) ───────────────────


def _sanitize_peak(p: str) -> str:
    """Mirror extract_fastas.sanitize: peak_ids become filenames, `#`→`_`."""
    return p.replace("#", "_")


def _run_peak_ids():
    """Peak set resolved at parse time (a checkpoint's re-evaluation hangs at 100K+ peaks).
    Priority: preextracted_dir *.fa, the survivor manifest fastas/peak_ids.txt (resume),
    config peak_ids (fresh run), then a glob of fastas/."""
    pre = config.get("preextracted_dir")
    if pre:
        src = Path(pre)
        return sorted(
            _sanitize_peak(fa.stem)
            for fa in src.glob("*.fa")
            if (src / f"{fa.stem}.tsv").exists()
        )
    manifest = Path(PEAKS)
    if manifest.exists():
        return [ln.strip() for ln in manifest.read_text().splitlines() if ln.strip()]
    pi = config.get("peak_ids")
    if pi and Path(pi).exists():
        return [
            _sanitize_peak(ln.strip())
            for ln in Path(pi).read_text().splitlines()
            if ln.strip()
        ]
    return sorted(_sanitize_peak(fa.stem) for fa in Path(RESULTS).glob("fastas/*.fa"))


def _stage1_targets():
    """Persisted Stage-1 outputs that analysis.smk consumes; Stage 2 never re-enters this DAG."""
    peaks = _run_peak_ids()
    return (
        expand(RESULTS + "/aligned/{peak_id}.fa", peak_id=peaks)
        + expand(RESULTS + "/trees/{peak_id}.rooted.nwk", peak_id=peaks)
        + expand(RESULTS + "/phylo_mm/{peak_id}.tsv", peak_id=peaks)
        + [HAP]
        + rules.stage1_tables.output
    )


def _extract_fastas_input(wildcards):
    """Either the bgzipped union BED, or the pre-extracted directory."""
    pre = config.get("preextracted_dir")
    if pre:
        return str(Path(pre).resolve())
    return config["input_table"]


def _peak_ids_input(wildcards):
    """Track the peak_ids file so switching it re-runs extraction."""
    pi = config.get("peak_ids")
    return [pi] if pi else []


def _highlight_input(wildcards):
    """Track the highlight TSV so changing it re-runs extraction."""
    hf = config.get("highlight_file")
    return [hf] if hf else []


def _mafft_input(wildcards):
    """RE-only FASTA, or the flanked FASTA when flanking is enabled."""
    if int(config.get("flank", 0)) > 0:
        return f"{TEMP}/flanked/{wildcards.peak_id}.fa"
    return f"{RESULTS}/fastas/{wildcards.peak_id}.fa"


def _re_tree(wildcards):
    """RE-only rooted tree: the trimmed-then-rebuilt tree when flank>0, else
    the main tree (with no flank the RE *is* the whole alignment)."""
    if config.get("flank", 0) > 0:
        return f"{TEMP}/re/trees/{wildcards.peak_id}.rooted.nwk"
    return f"{RESULTS}/trees/{wildcards.peak_id}.rooted.nwk"


def _re_alignment(wildcards):
    """RE-only alignment paired with `_re_tree`."""
    if config.get("flank", 0) > 0:
        return f"{TEMP}/re/aligned/{wildcards.peak_id}.fa"
    return f"{RESULTS}/aligned/{wildcards.peak_id}.fa"
