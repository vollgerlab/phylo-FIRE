"""Funnel counts for the processing figure, read from a finished run's tables.

    pixi run python docs/figures/processing_funnel/funnel_counts.py [run_name]

writes counts.json (the snapshot the draft concepts read); the final Sankey
imports counts() directly so it always reflects the current tables.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import polars as pl

HSCA_INDEX = Path("resources/hsca/hsca-index.tsv")
N_INDIVIDUALS = 39


def tables(run: str) -> Path:
    return Path("results") / run / "tables"


def counts(run: str = "hsca_hapsel") -> dict[str, int]:
    """Funnel counts from the HSCA index and the run's result tables."""
    hsca = pl.read_csv(HSCA_INDEX, separator="\t")
    cls = dict(hsca["haplotype_distinct"].value_counts().iter_rows())

    rv = pl.read_parquet(tables(run) / "re_variants.parquet")
    pk = rv.group_by("peak_id").agg(
        pl.first("attribution_mode"), pl.first("n_variants_selected")
    )
    mode = dict(pk["attribution_mode"].value_counts().iter_rows())
    n_drv_peaks = pk.filter(pl.col("n_variants_selected") > 0).height

    var = pl.read_csv(tables(run) / "variants.hg38.tsv.gz", separator="\t",
                      infer_schema_length=0, null_values=["NA", ""],
                      columns=["peak_id", "is_driver", "hg38_pos"])
    drv_snv = var.filter(pl.col("is_driver") == "TRUE").height
    placed = var.filter((pl.col("is_driver") == "TRUE") & pl.col("hg38_pos").is_not_null())
    drv_snv_hg38 = placed.height
    drv_peaks_hg38 = placed["peak_id"].n_unique()
    unlifted = pl.read_csv(tables(run) / "drivers.unlifted.tsv", separator="\t").height

    with gzip.open(tables(run) / "variants.hg38.rsid.vcf.gz", "rt") as fh:
        drv_rsid = sum(1 for ln in fh if not ln.startswith("#") and ";DRIVER" in ln and ln.split("\t", 3)[2].startswith("rs"))

    heritable = mode.get("RE", 0) + mode.get("region", 0)
    return {
        "hsca": hsca.height,
        "hapsel": cls["Haplotype-selective"],
        "epi": cls.get("Epigenetically haplotype-selective", 0),
        "gen": cls.get("Genetically haplotype-specific", 0),
        "complex": cls.get("Complex", 0),
        "peaks": pk.height,
        "unique_seqs": rv.height,
        "var_sites": var.height,
        "heritable": heritable,
        "her_re": mode.get("RE", 0),
        "her_region": mode.get("region", 0),
        "not_heritable": mode.get("none", 0),
        "drv_peaks": n_drv_peaks,
        "her_no_drv": heritable - n_drv_peaks,
        "drv_sets": int(pk["n_variants_selected"].sum()),
        "multi_drv_peaks": pk.filter(pl.col("n_variants_selected") > 1).height,
        "drv_snv": drv_snv,
        "drv_snv_hg38": drv_snv_hg38,
        "drv_peaks_hg38": drv_peaks_hg38,
        "drv_rsid": drv_rsid,
        "unlifted": unlifted,
        "n_individuals": N_INDIVIDUALS,
        "n_haplotypes": 2 * N_INDIVIDUALS,
    }


if __name__ == "__main__":
    c = counts(sys.argv[1] if len(sys.argv) > 1 else "hsca_hapsel")
    (Path(__file__).resolve().parent / "counts.json").write_text(json.dumps(c, indent=1))
    print(c)
