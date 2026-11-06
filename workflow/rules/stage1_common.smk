"""Stage-1 helpers: wildcard constraints, the parse-time peak set and the input functions."""

import shlex


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
    manifest = Path(workflow.pathvars.apply(PEAKS))
    if manifest.exists():
        return [ln.strip() for ln in manifest.read_text().splitlines() if ln.strip()]
    pi = config.get("peak_ids")
    if pi and Path(pi).exists():
        return [
            _sanitize_peak(ln.strip())
            for ln in Path(pi).read_text().splitlines()
            if ln.strip()
        ]
    return sorted(_sanitize_peak(fa.stem) for fa in Path(workflow.pathvars.apply("<results>/fastas")).glob("*.fa"))


PEAK_IDS = _run_peak_ids()


def _extract_fastas_input(wildcards):
    """Either the bgzipped union BED, or the pre-extracted directory."""
    pre = config.get("preextracted_dir")
    if pre:
        return str(Path(pre).resolve())
    return config["input_table"]


def optional_file(key: str):
    """Input function tracking config[key] when set, so changing that file re-runs the rule."""
    return lambda wildcards: [config[key]] if config.get(key) else []


def _extract_opt_args(wildcards, input) -> str:
    """Optional extract_fastas.py flags from config and the optional inputs."""
    args = []
    if config.get("primary_only", True):
        args.append("--primary-only")
    if input.peak_ids:
        args.append(f"--peak-ids {shlex.quote(str(input.peak_ids[0]))}")
    if config.get("haplotype_distinct"):
        args.append(
            f"--haplotype-distinct {shlex.quote(str(config['haplotype_distinct']))}"
        )
    if input.highlight:
        args.append(f"--highlight-file {shlex.quote(str(input.highlight[0]))}")
    return " ".join(args)


def iqtree_bootstrap(wildcards, input) -> str:
    """Ultrafast bootstrap only with >=4 unique sequences; IQ-TREE deduplicates
    internally and refuses bootstrap on fewer distinct taxa."""
    return "-B 1000" if count_unique_seqs_phylip(str(input.nuc)) >= 4 else ""


def _mafft_input(wildcards):
    """RE-only FASTA, or the flanked FASTA when flanking is enabled."""
    if FLANKED:
        return expand(rules.add_flanks.output.fa, peak_id=wildcards.peak_id)
    # extract_fastas declares the expanded list, so there is no pattern to reference.
    return f"<results>/fastas/{wildcards.peak_id}.fa"


def _re_tree(wildcards):
    """RE-only rooted tree: the trimmed-then-rebuilt tree when flank>0, else
    the main tree (with no flank the RE *is* the whole alignment)."""
    out = rules.re_root_tree.output if FLANKED else rules.root_tree.output
    return expand(out.rooted, peak_id=wildcards.peak_id)


def _re_alignment(wildcards):
    """RE-only alignment paired with `_re_tree`."""
    out = rules.re_trim.output if FLANKED else rules.mafft.output
    return expand(out.fa, peak_id=wildcards.peak_id)
