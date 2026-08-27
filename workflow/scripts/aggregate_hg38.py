#!/usr/bin/env python3
"""Concatenate per-peak project_variants TSVs into the hg38-projection deliverables:

  * variants.hg38.tsv.gz  — every variable site in every peak (all columns).
  * variants.hg38.vcf     — sites-only VCF in hg38 coords (REF = hg38 reference
                            base, ALT = cohort allele(s); DRIVER/BETA in INFO).
  * variants.hg38.bed     — BED of the same sites (name = peak:RE:REF>ALT), with
                            is_driver + beta (the GWAS foreground/background list).
  * elements.hg38.bed     — one interval per peak (hg38), with has_driver +
                            n_variants + n_drivers (the element-level list).

Sites without an hg38 coordinate (peak not in hg38, or the site sits in an
insertion relative to hg38) are kept in the TSV but omitted from the VCF/BED.
"""
import argparse
import gzip
from pathlib import Path

VCOLS = ["peak_id", "re_pos", "cons_chr", "cons_pos", "cons_ref",
         "hg38_chr", "hg38_pos", "hg38_ref", "hap_major", "hap_minor",
         "driver_allele", "n_carriers", "is_driver", "beta", "carriers",
         "asm_coords"]

VCF_HEADER = """##fileformat=VCFv4.2
##source=phylo-FIRE/project_variants
##INFO=<ID=PEAK,Number=1,Type=String,Description="consensus peak id">
##INFO=<ID=RE,Number=1,Type=Integer,Description="1-based position within the RE/peak">
##INFO=<ID=DRIVER,Number=0,Type=Flag,Description="attribution driver variant">
##INFO=<ID=DALLELE,Number=1,Type=Character,Description="driver/effect allele (the phylo_mm clade allele BETA is for; may be REF or ALT, major or minor)">
##INFO=<ID=BETA,Number=1,Type=Float,Description="joint driver effect (actuation shift) of DALLELE">
##INFO=<ID=NCAR,Number=1,Type=Integer,Description="haplotype carriers of the cohort minor allele">
##INFO=<ID=CONS,Number=1,Type=String,Description="CHM13 coordinate chr:pos">
##INFO=<ID=CONSREF,Number=1,Type=Character,Description="CHM13 reference base">
##INFO=<ID=HAPMAJ,Number=1,Type=Character,Description="cohort major allele">
##INFO=<ID=HAPMIN,Number=1,Type=Character,Description="cohort minor (variant) allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def _alt(ref: str, *alleles: str) -> str:
    """ALT alleles = the cohort alleles that differ from the reference base.
    Pass the driver allele too: at a >2-allele column the driver's clade can carry
    a THIRD base that isn't the column major/minor, and it must appear in ALT so
    the reported variant actually contains the driver."""
    seen, out = set(), []
    for a in alleles:
        if a and a != ref and a not in ("-", "N", "NA") and a not in seen:
            seen.add(a); out.append(a)
    return ",".join(out) if out else "."


def _chrom_key(c: str):
    """Natural-ish chrom sort: chr1..chr22, chrX, chrY, chrM, then the rest."""
    s = c[3:] if c.startswith("chr") else c
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, int(s)) if s.isdigit() else (1, order.get(s, 99))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--coords", required=True, type=Path, help="peak_coords.tsv")
    ap.add_argument("--out-variants", required=True, type=Path)
    ap.add_argument("--out-vcf", required=True, type=Path)
    ap.add_argument("--out-variants-bed", required=True, type=Path)
    ap.add_argument("--out-elements-bed", required=True, type=Path)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    rows, missing = [], 0
    for peak in peaks:
        f = args.indir / f"{peak}.tsv"
        if not f.exists():
            missing += 1
            continue
        for ln in f.read_text().splitlines()[1:]:       # skip header
            if ln.strip():
                rows.append(dict(zip(VCOLS, ln.split("\t"))))
    if missing:
        print(f"[aggregate_hg38] {missing}/{len(peaks)} peak projections missing", flush=True)

    args.out_variants.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.out_variants, "wt") as fh:      # full table
        fh.write("\t".join(VCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in VCOLS) + "\n")

    # sites with a usable hg38 coordinate + single-base reference (SNV)
    sites = [r for r in rows
             if r["hg38_chr"] != "NA" and r["hg38_pos"].lstrip("-").isdigit()
             and r["hg38_ref"] in ("A", "C", "G", "T")]
    # Include the exact chrom string: many hg38 alt/random contigs share the same
    # _chrom_key (1, 99), so without it distinct contigs interleave by position and
    # tabix rejects the file ("Chromosome blocks not continuous").
    sites.sort(key=lambda r: (_chrom_key(r["hg38_chr"]), r["hg38_chr"], int(r["hg38_pos"])))

    # VCF
    n_vcf = 0
    with args.out_vcf.open("w") as fh:
        fh.write(VCF_HEADER)
        for r in sites:
            alt = _alt(r["hg38_ref"], r["hap_major"], r["hap_minor"], r.get("driver_allele", "NA"))
            info = [f"PEAK={r['peak_id']}", f"RE={r['re_pos']}", f"NCAR={r['n_carriers']}",
                    f"CONS={r['cons_chr']}:{r['cons_pos']}", f"CONSREF={r['cons_ref']}",
                    f"HAPMAJ={r['hap_major']}", f"HAPMIN={r['hap_minor']}"]
            if r["is_driver"] == "TRUE":
                info.append("DRIVER")
                if r.get("driver_allele", "NA") not in ("NA", ""):
                    info.append(f"DALLELE={r['driver_allele']}")
                if r["beta"] not in ("NA", ""):
                    info.append(f"BETA={r['beta']}")
            vid = f"{r['peak_id']}:RE{r['re_pos']}"
            fh.write("\t".join([r["hg38_chr"], r["hg38_pos"], vid, r["hg38_ref"],
                                alt, ".", "PASS", ";".join(info)]) + "\n")
            n_vcf += 1

    # variant BED (same sites)
    with args.out_variants_bed.open("w") as fh:
        for r in sites:
            pos = int(r["hg38_pos"])
            alt = _alt(r["hg38_ref"], r["hap_major"], r["hap_minor"], r.get("driver_allele", "NA"))
            name = f"{r['peak_id']}:RE{r['re_pos']}:{r['hg38_ref']}>{alt}"
            fh.write("\t".join([r["hg38_chr"], str(pos - 1), str(pos), name,
                                r["is_driver"], r["beta"], r["n_carriers"], r["peak_id"],
                                r.get("driver_allele", "NA")]) + "\n")

    # per-peak counts (over ALL rows, not just hg38-mapped)
    n_var, n_drv = {}, {}
    for r in rows:
        p = r["peak_id"]
        n_var[p] = n_var.get(p, 0) + 1
        if r["is_driver"] == "TRUE":
            n_drv[p] = n_drv.get(p, 0) + 1

    # element BED from the coords table (peaks that exist in hg38)
    n_elem = 0
    clines = args.coords.read_text().splitlines()
    ci = {h: i for i, h in enumerate(clines[0].split("\t"))}
    with args.out_elements_bed.open("w") as fh:
        for ln in clines[1:]:
            c = ln.split("\t")
            if c[ci["exists_in_hg38"]].upper() not in ("TRUE", "1"):
                continue
            pid = c[ci["peak_id"]]
            fh.write("\t".join([
                c[ci["hg38_chr"]], c[ci["hg38_start"]], c[ci["hg38_end"]], pid,
                "TRUE" if n_drv.get(pid, 0) else "FALSE",
                str(n_var.get(pid, 0)), str(n_drv.get(pid, 0)),
            ]) + "\n")
            n_elem += 1

    print(f"[aggregate_hg38] {len(rows)} sites "
          f"({sum(1 for r in rows if r['is_driver']=='TRUE')} drivers) over {len(n_var)} "
          f"peaks → {n_vcf} in VCF/BED, {n_elem} elements")


if __name__ == "__main__":
    main()
