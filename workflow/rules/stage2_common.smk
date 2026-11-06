"""Stage-2 helpers (workflow/analysis.smk): opt-in flags, the target list and
input/param functions. Include after common.smk."""

PLOT_HTML = cfg_bool("plot_html", False)
PROJECT_HG38 = cfg_bool("project_hg38", True)
PLASMID = cfg_bool("plasmid", False)


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
    """plot_all's --pclai-bed flag, or "" when config pclai_bed is unset or missing."""
    p = config.get("pclai_bed")
    return f"--pclai-bed {p}" if p and Path(p).exists() else ""


def _targets(_wildcards):
    """Stage-2 targets (re_variants and plots) plus the opt-in hg38, rsID, plasmid and
    PheWAS outputs whose config/inputs are present."""
    out = [
        rules.re_variants_all.output.parquet,
        rules.re_variants_all.output.tsv,
        rules.plot_all.output.done,
    ]
    if _hg38_ok():
        p = rules.project_all.output
        out += [
            p.variants,
            p.vcf,
            p.variants_bed,
            p.variants_chm13_bed,
            p.elements_bed,
            p.unlifted,
        ]
        if config.get("dbsnp"):
            out.append(rules.annotate_rsids.output.vcf)
    if PLASMID:
        out += [rules.plasmid_all.output.constructs, rules.plasmid_all.output.fa]
    if config.get("significant_loci") and config.get("dbsnp") and _hg38_ok():
        out += [
            rules.significant_driver_rsids.output.tsv,
            rules.significant_driver_rsids.output.txt,
        ]
    return out
