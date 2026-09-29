"""Post-phylo-FIRE -- local follow-up on a finished cluster run.

Downloads results/{run}/tables from Redwood, annotates the hg38 driver VCF with
dbSNP rs IDs (dbSNP subset streamed from Redwood at our sites only), and maps
the significant plasmid fiber-seq loci to their driver variants for PheWAS.
Needs an open ssh ControlMaster to `remote_host` (log in first).

    pixi run post -- --config significant_loci=<plasmid_significant_loci.tsv>
"""

import sys
from pathlib import Path
from snakemake.utils import min_version

min_version("8.0")


configfile: "config/post.yaml"


RUN = config["run_name"]
T = f"results/{RUN}/tables"
HOST = config["remote_host"]
REMOTE_T = f"{config['remote_repo']}/{T}"
PULLED = ["variants.hg38.tsv.gz", "variants.hg38.vcf.gz", "variants.hg38.vcf.gz.tbi", "drivers.unlifted.tsv"]

localrules: all, pull_tables, dbsnp_subset


rule all:
    input:
        f"results/{RUN}/phewas/significant_driver_rsids.tsv",
        f"results/{RUN}/phewas/significant_driver_rsids.txt",


rule pull_tables:
    """rsync the Stage-2 tables this workflow reads (download only)."""
    output:
        [f"{T}/{f}" for f in PULLED],
    params:
        src=" ".join(f"{HOST}:{REMOTE_T}/{f}" for f in PULLED),
    shell:
        "rsync -a {params.src} {T}/"


rule dbsnp_subset:
    """Stream dbSNP records at our variant sites from Redwood (RefSeq contig names)."""
    input:
        vcf=f"{T}/variants.hg38.vcf.gz",
    output:
        vcf=f"resources/dbsnp.{RUN}.nc.vcf.gz",
        tbi=f"resources/dbsnp.{RUN}.nc.vcf.gz.tbi",
    params:
        remote=f"{config['remote_bcftools']} view -R <(cat) -Oz {config['remote_dbsnp']}",
    run:
        sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
        from annotate_rsids import CHR2NC

        sites = shell(f"bcftools query -f '%CHROM\\t%POS\\n' {input.vcf} | sort -u", read=True)
        regions = Path(output.vcf).with_suffix(".regions.tsv")
        rows = sorted((CHR2NC[c], int(p)) for c, p in (l.split("\t") for l in sites.splitlines()) if c in CHR2NC)
        regions.write_text("".join(f"{c}\t{p}\t{p}\n" for c, p in rows))
        shell(f"ssh -o BatchMode=yes {HOST} '{params.remote}' < {regions} > {output.vcf}")
        shell(f"tabix -f -p vcf {output.vcf}")
        regions.unlink()


rule annotate_rsids:
    input:
        vcf=f"{T}/variants.hg38.vcf.gz",
        dbsnp=rules.dbsnp_subset.output.vcf,
    output:
        vcf=f"{T}/variants.hg38.rsid.vcf.gz",
    threads: 4
    params:
        script=Path(workflow.basedir) / "scripts/annotate_rsids.py",
    shell:
        "python {params.script} --vcf {input.vcf} --dbsnp {input.dbsnp} "
        "--out-vcf {output.vcf} --threads {threads}"


rule significant_driver_rsids:
    input:
        sig=config["significant_loci"] or "significant_loci_not_set",
        tsv=f"{T}/variants.hg38.tsv.gz",
        vcf=rules.annotate_rsids.output.vcf,
        unlifted=f"{T}/drivers.unlifted.tsv",
    output:
        tsv=f"results/{RUN}/phewas/significant_driver_rsids.tsv",
        txt=f"results/{RUN}/phewas/significant_driver_rsids.txt",
    params:
        script=Path(workflow.basedir) / "scripts/significant_driver_rsids.py",
    shell:
        "python {params.script} --sig {input.sig} --variants {input.tsv} "
        "--rsid-vcf {input.vcf} --unlifted {input.unlifted} "
        "--out-tsv {output.tsv} --out-rsids {output.txt}"
