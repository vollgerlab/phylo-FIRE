#!/usr/bin/env python3
"""Concatenate per-peak project_variants TSVs into the hg38-projection deliverables:

  * variants.hg38.tsv.gz  — every variable site in every peak (all columns).
  * variants.hg38.bed     — BED of sites that carry an hg38 coordinate, with
                            is_driver + beta (the GWAS foreground/background list).
  * elements.hg38.bed     — one interval per peak (hg38), with has_driver +
                            n_variants + n_drivers (the element-level list).

Sites without an hg38 coordinate (peak not in hg38, or the site sits in an
insertion relative to hg38) are kept in the TSV but omitted from the BED.
"""
import argparse
import gzip
from pathlib import Path

VCOLS = ["peak_id", "re_pos", "cons_chr", "cons_pos", "hg38_chr", "hg38_pos",
         "ref_allele", "alt_allele", "n_carriers", "is_driver", "beta",
         "carriers", "asm_coords"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--coords", required=True, type=Path, help="peak_coords.tsv")
    ap.add_argument("--out-variants", required=True, type=Path)
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
        lines = f.read_text().splitlines()
        for ln in lines[1:]:                       # skip header
            if ln.strip():
                rows.append(dict(zip(VCOLS, ln.split("\t"))))
    if missing:
        print(f"[aggregate_hg38] {missing}/{len(peaks)} peak projections missing", flush=True)

    args.out_variants.parent.mkdir(parents=True, exist_ok=True)
    # full variant table (gzipped TSV)
    with gzip.open(args.out_variants, "wt") as fh:
        fh.write("\t".join(VCOLS) + "\n")
        for r in rows:
            fh.write("\t".join(r[c] for c in VCOLS) + "\n")

    # variant BED (hg38-coord sites only)
    n_bed = 0
    with args.out_variants_bed.open("w") as fh:
        for r in rows:
            if r["hg38_chr"] == "NA" or not r["hg38_pos"].lstrip("-").isdigit():
                continue
            pos = int(r["hg38_pos"])                # 1-based
            name = f"{r['peak_id']}:RE{r['re_pos']}:{r['ref_allele']}>{r['alt_allele']}"
            fh.write("\t".join([
                r["hg38_chr"], str(pos - 1), str(pos), name,
                r["is_driver"], r["beta"], r["n_carriers"], r["peak_id"],
            ]) + "\n")
            n_bed += 1

    # per-peak counts
    n_var: dict[str, int] = {}
    n_drv: dict[str, int] = {}
    for r in rows:
        p = r["peak_id"]
        n_var[p] = n_var.get(p, 0) + 1
        if r["is_driver"] == "TRUE":
            n_drv[p] = n_drv.get(p, 0) + 1

    # element BED from the coords table (peaks that exist in hg38)
    n_elem = 0
    clines = args.coords.read_text().splitlines()
    chdr = clines[0].split("\t")
    ci = {h: i for i, h in enumerate(chdr)}
    with args.out_elements_bed.open("w") as fh:
        for ln in clines[1:]:
            c = ln.split("\t")
            if c[ci["exists_in_hg38"]].upper() not in ("TRUE", "1"):
                continue
            pid = c[ci["peak_id"]]
            has_drv = "TRUE" if n_drv.get(pid, 0) > 0 else "FALSE"
            fh.write("\t".join([
                c[ci["hg38_chr"]], c[ci["hg38_start"]], c[ci["hg38_end"]], pid,
                has_drv, str(n_var.get(pid, 0)), str(n_drv.get(pid, 0)),
            ]) + "\n")
            n_elem += 1

    print(f"[aggregate_hg38] {len(rows)} sites ({sum(1 for r in rows if r['is_driver']=='TRUE')} "
          f"drivers) over {len(n_var)} peaks → {n_bed} in variant BED, {n_elem} elements")


if __name__ == "__main__":
    main()
