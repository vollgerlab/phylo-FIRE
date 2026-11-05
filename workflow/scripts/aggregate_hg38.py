#!/usr/bin/env python3
"""Concatenate per-peak project_variants TSVs into the hg38-projection deliverables:

  * variants.hg38.tsv.gz  — every variable site in every peak (all columns).
  * variants.hg38.vcf     — sites-only VCF in hg38 coords (REF = hg38 reference
                            base, ALT = cohort allele(s); DRIVER/BETA in INFO).
  * variants.hg38.bed     — BED of the same sites (name = peak:RE:REF>ALT), with
                            is_driver + beta + divergent (the GWAS foreground/background list).
  * variants.chm13.bed    — the same BED in CHM13 coordinates (cons_chr/cons_pos/
                            cons_ref; sites on CHM13 chromosomes only).
  * elements.hg38.bed     — one interval per peak (hg38), with has_driver +
                            n_variants + n_drivers (the element-level list).

Sites without an hg38 coordinate (peak not in hg38, or the site sits in an
insertion relative to hg38) are kept in the TSV but omitted from the VCF/BED.
Called in-process by project_variants_all.py.
"""
import gzip
import re
from pathlib import Path

VCOLS = ["peak_id", "re_pos", "cons_chr", "cons_pos", "cons_ref",
         "hg38_chr", "hg38_pos", "hg38_ref", "hap_major", "hap_minor",
         "driver_allele", "n_carriers", "is_driver", "beta", "carriers",
         "asm_coords", "divergent"]

VCF_HEADER = """##fileformat=VCFv4.2
##source=phylo-FIRE/project_variants
##INFO=<ID=PEAK,Number=1,Type=String,Description="consensus peak id">
##INFO=<ID=RE,Number=1,Type=Integer,Description="1-based position within the RE/peak">
##INFO=<ID=DRIVER,Number=0,Type=Flag,Description="attribution driver variant">
##INFO=<ID=DIVERGENT,Number=0,Type=Flag,Description="driver set covers >5% of the attributed band (divergent haplotype, not a single clonable variant); flagged, not filtered">
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
            seen.add(a)
            out.append(a)
    return ",".join(out) if out else "."


_CHM13_CHR = re.compile(r"^chr([0-9]+|X|Y|M)$")


def _write_variant_bed(path: Path, rows: list[dict], chr_k: str, pos_k: str, ref_k: str) -> int:
    """Write the variant BED for one coordinate system (hg38 or CHM13 columns).

    Keeps SNV sites with a coordinate and a single-base reference, sorted by
    chrom then position; returns the number of sites written."""
    sites = [r for r in rows
             if r[chr_k] != "NA" and r[pos_k].lstrip("-").isdigit()
             and r[ref_k] in ("A", "C", "G", "T")]
    sites.sort(key=lambda r: (_chrom_key(r[chr_k]), r[chr_k], int(r[pos_k])))
    with path.open("w") as fh:
        for r in sites:
            pos = int(r[pos_k])
            alt = _alt(r[ref_k], r["hap_major"], r["hap_minor"], r.get("driver_allele", "NA"))
            name = f"{r['peak_id']}:RE{r['re_pos']}:{r[ref_k]}>{alt}"
            fh.write("\t".join([r[chr_k], str(pos - 1), str(pos), name,
                                r["is_driver"], r["beta"], r["n_carriers"], r["peak_id"],
                                r.get("driver_allele", "NA"), r.get("divergent", "FALSE")]) + "\n")
    return len(sites)


def _chrom_key(c: str):
    """Natural-ish chrom sort: chr1..chr22, chrX, chrY, chrM, then the rest."""
    s = c[3:] if c.startswith("chr") else c
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (0, int(s)) if s.isdigit() else (1, order.get(s, 99))


def load_rows(indir: Path, peaks: list[str]) -> list[dict]:
    """Every row of ``{indir}/{peak}.tsv`` (header skipped) as {VCOLS: str}."""
    rows, missing = [], 0
    for peak in peaks:
        f = indir / f"{peak}.tsv"
        if not f.exists():
            missing += 1
            continue
        for ln in f.read_text().splitlines()[1:]:       # skip header
            if ln.strip():
                rows.append(dict(zip(VCOLS, ln.split("\t"))))
    if missing:
        print(f"[aggregate_hg38] {missing}/{len(peaks)} peak projections missing", flush=True)
    return rows


def write_variants_tsv(rows: list[dict], path: Path) -> None:
    """The full gzipped variant table (every site, every column)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(VCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in VCOLS) + "\n")


def write_vcf(rows: list[dict], path: Path) -> int:
    """Sites-only hg38 VCF of the SNV sites with an hg38 coordinate; returns the
    number of records."""
    sites = [r for r in rows
             if r["hg38_chr"] != "NA" and r["hg38_pos"].lstrip("-").isdigit()
             and r["hg38_ref"] in ("A", "C", "G", "T")]
    # Include the exact chrom string: many hg38 alt/random contigs share the same
    # _chrom_key (1, 99), so without it distinct contigs interleave by position and
    # tabix rejects the file ("Chromosome blocks not continuous").
    sites.sort(key=lambda r: (_chrom_key(r["hg38_chr"]), r["hg38_chr"], int(r["hg38_pos"])))
    with path.open("w") as fh:
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
                if r.get("divergent") == "TRUE":
                    info.append("DIVERGENT")
            vid = f"{r['peak_id']}:RE{r['re_pos']}"
            fh.write("\t".join([r["hg38_chr"], r["hg38_pos"], vid, r["hg38_ref"],
                                alt, ".", "PASS", ";".join(info)]) + "\n")
    return len(sites)


def count_per_peak(rows: list[dict]) -> tuple[dict[str, int], dict[str, int]]:
    """(variants per peak, drivers per peak) over ALL rows, not just hg38-mapped."""
    n_var, n_drv = {}, {}
    for r in rows:
        p = r["peak_id"]
        n_var[p] = n_var.get(p, 0) + 1
        if r["is_driver"] == "TRUE":
            n_drv[p] = n_drv.get(p, 0) + 1
    return n_var, n_drv


def write_elements_bed(coords: Path, n_var: dict, n_drv: dict, path: Path) -> int:
    """One hg38 interval per peak that exists in hg38 (from the coords table),
    with has_driver, n_variants and n_drivers; returns the number written."""
    n_elem = 0
    clines = coords.read_text().splitlines()
    ci = {h: i for i, h in enumerate(clines[0].split("\t"))}
    with path.open("w") as fh:
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
    return n_elem


def aggregate(indir: Path, peaks: list[str], coords: Path, out_variants: Path,
              out_vcf: Path, out_variants_bed: Path, out_elements_bed: Path,
              out_variants_chm13_bed: Path | None = None) -> None:
    """Write the variant table, plain VCF and BEDs from the per-peak TSVs."""
    rows = load_rows(Path(indir), peaks)
    write_variants_tsv(rows, Path(out_variants))
    n_vcf = write_vcf(rows, Path(out_vcf))
    _write_variant_bed(Path(out_variants_bed), rows, "hg38_chr", "hg38_pos", "hg38_ref")
    if out_variants_chm13_bed:
        chm13_rows = [r for r in rows if _CHM13_CHR.match(r["cons_chr"])]
        n_chm13 = _write_variant_bed(Path(out_variants_chm13_bed), chm13_rows,
                                     "cons_chr", "cons_pos", "cons_ref")
        print(f"[aggregate_hg38] {n_chm13} sites in CHM13 BED", flush=True)
    n_var, n_drv = count_per_peak(rows)
    n_elem = write_elements_bed(Path(coords), n_var, n_drv, Path(out_elements_bed))
    print(f"[aggregate_hg38] {len(rows)} sites "
          f"({sum(1 for r in rows if r['is_driver'] == 'TRUE')} drivers) over {len(n_var)} "
          f"peaks → {n_vcf} in VCF/BED, {n_elem} elements")
