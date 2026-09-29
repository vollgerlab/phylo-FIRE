#!/usr/bin/env python3
"""Map significant plasmid fiber-seq loci to their driver variants' hg38 coords + rs IDs.

One output row per driver variant. The variant tables carry no driver_id, so a
driver is keyed on peak + its full-precision beta from variants.hg38.tsv.gz
(drivers within a peak have distinct betas; the VCF's BETA is rounded, so it
is not used as a key). rs IDs come from the annotated VCF via PEAK + RE.
Use `driver_allele` as the effect allele: it is often the hg38 REF.

    python workflow/scripts/significant_driver_rsids.py --sig <loci.tsv> \
        --variants variants.hg38.tsv.gz --rsid-vcf variants.hg38.rsid.vcf.gz \
        --unlifted drivers.unlifted.tsv --out-tsv out.tsv --out-rsids out.txt
"""
import argparse
import io
import subprocess
import sys
from pathlib import Path

import polars as pl

SIG_COLS = ["locus", "driver_id", "category", "beta", "driver_variant", "only_driver", "p", "q", "delta_fire"]
VAR_COLS = ["peak_id", "re_pos", "hg38_chr", "hg38_pos", "hg38_ref", "driver_allele", "n_carriers", "is_driver", "beta"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig", required=True, type=Path)
    ap.add_argument("--variants", required=True, type=Path)
    ap.add_argument("--rsid-vcf", required=True, type=Path)
    ap.add_argument("--unlifted", required=True, type=Path)
    ap.add_argument("--out-tsv", required=True, type=Path)
    ap.add_argument("--out-rsids", required=True, type=Path)
    a = ap.parse_args()

    sig = pl.read_csv(a.sig, separator="\t", infer_schema=False, columns=SIG_COLS)
    drv = (pl.read_csv(a.variants, separator="\t", infer_schema=False, columns=VAR_COLS)
             .filter(pl.col("is_driver") == "TRUE").drop("is_driver").rename({"peak_id": "locus"}))
    q = subprocess.run(["bcftools", "query", "-f", r"%PEAK\t%RE\t%ID\t%REF\t%ALT\n", str(a.rsid_vcf)],
                       capture_output=True, text=True, check=True).stdout
    rs = (pl.read_csv(io.StringIO(q), separator="\t", has_header=False, infer_schema=False,
                      new_columns=["locus", "re_pos", "rsid", "vcf_ref", "vcf_alt"])
            .with_columns(pl.col("rsid").replace(".", None)))
    unl = (pl.read_csv(a.unlifted, separator="\t", infer_schema=False, columns=["peak_id", "beta", "reason"])
             .rename({"peak_id": "locus", "reason": "unlifted_reason"}))

    out = (sig.join(drv, on=["locus", "beta"], how="left")
              .join(rs, on=["locus", "re_pos"], how="left")
              .join(unl, on=["locus", "beta"], how="left"))
    assert out.select("locus", "driver_id").unique().height == sig.select("locus", "driver_id").unique().height

    a.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out.write_csv(a.out_tsv, separator="\t", null_value="NA")
    rsids = out["rsid"].drop_nulls().unique().sort().to_list()
    a.out_rsids.write_text("".join(f"{r}\n" for r in rsids))

    per = out.group_by("locus", "driver_id").agg(
        n_var=pl.col("re_pos").drop_nulls().len(),
        n_hg38=pl.col("hg38_pos").replace("NA", None).drop_nulls().len(),
        n_rs=pl.col("rsid").drop_nulls().len(),
        unl=pl.col("unlifted_reason").first())
    status = (pl.when(pl.col("n_rs") > 0).then(pl.lit("has rsID"))
                .when(pl.col("n_var") == 0).then(pl.concat_str(pl.lit("no driver row: "), pl.col("unl").fill_null("?")))
                .when(pl.col("n_hg38") == 0).then(pl.lit("no hg38 position"))
                .otherwise(pl.lit("hg38 position, no rsID")))
    counts = per.with_columns(status=status).group_by("status").len().sort("len", descending=True)
    print(f"[significant_driver_rsids] {sig.height} drivers, {out.height} variant rows, {len(rsids)} rsIDs",
          file=sys.stderr)
    for s, n in counts.iter_rows():
        print(f"  {s}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
