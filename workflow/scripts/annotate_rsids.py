#!/usr/bin/env python3
"""Annotate the hg38 variants VCF with dbSNP rs IDs (Stage 2).

`bcftools annotate -c ID` matches on CHROM+POS+REF+ALT (allele-aware), so it
assigns the correct rs (no coordinate-only / off-by-one issues). Reconciles
chr <-> RefSeq contig naming so either an NCBI dbSNP (NC_*.*, e.g.
GCF_000001405.40.gz, local or a remote-indexed subset) or a chr-named dbSNP
works. Emits variants.hg38.rsid.vcf.gz (+ .tbi).

    pixi run python workflow/scripts/annotate_rsids.py \
        --vcf results/<run>/tables/variants.hg38.vcf.gz \
        --dbsnp <dbsnp.vcf.gz> --out-vcf results/<run>/tables/variants.hg38.rsid.vcf.gz
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

# GRCh38 RefSeq chromosome accessions (dbSNP from NCBI uses these).
CHR2NC = {
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
    "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
    "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
}


def _run(cmd):
    subprocess.run([str(c) for c in cmd], check=True)


def _dbsnp_is_refseq(dbsnp) -> bool:
    hdr = subprocess.run(["bcftools", "view", "-h", str(dbsnp)],
                         capture_output=True, text=True).stdout
    return "ID=NC_0000" in hdr or "NC_000001" in hdr


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path, help="chr-named variants.hg38.vcf.gz")
    ap.add_argument("--dbsnp", required=True, type=Path, help="tabix-indexed dbSNP VCF")
    ap.add_argument("--out-vcf", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    work = Path(tempfile.mkdtemp(prefix="rsids_"))
    vcf = args.vcf
    back = None
    if _dbsnp_is_refseq(args.dbsnp):
        c2n = work / "chr2nc.txt"
        c2n.write_text("\n".join(f"{k}\t{v}" for k, v in CHR2NC.items()) + "\n")
        nc = work / "in.nc.vcf.gz"
        _run(["bcftools", "annotate", "--rename-chrs", c2n, args.vcf, "-Oz", "-o", nc])
        _run(["tabix", "-f", "-p", "vcf", nc])
        vcf = nc
        back = work / "nc2chr.txt"
        back.write_text("\n".join(f"{v}\t{k}" for k, v in CHR2NC.items()) + "\n")

    ann = work / "ann.vcf.gz"
    _run(["bcftools", "annotate", "-a", args.dbsnp, "-c", "ID", vcf,
          "-Oz", "-o", ann, "--threads", args.threads])
    args.out_vcf.parent.mkdir(parents=True, exist_ok=True)
    if back:
        _run(["bcftools", "annotate", "--rename-chrs", back, ann, "-Oz", "-o", args.out_vcf])
    else:
        args.out_vcf.write_bytes(ann.read_bytes())
    _run(["tabix", "-f", "-p", "vcf", args.out_vcf])

    lines = subprocess.run(["bcftools", "view", "-H", str(args.out_vcf)],
                           capture_output=True, text=True).stdout.splitlines()
    n_rs = sum(1 for ln in lines if ln.split("\t")[2].startswith("rs"))
    print(f"[annotate_rsids] {n_rs}/{len(lines)} variants got an rs ID -> "
          f"{args.out_vcf}", file=sys.stderr)


if __name__ == "__main__":
    main()
