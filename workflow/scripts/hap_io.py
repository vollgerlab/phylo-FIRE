"""Shared readers for the consolidated haplotypes.parquet (the Stage-2 source of
per-haplotype actuation + coords + RE sequence, replacing the fastas/ dir)."""
import polars as pl


def load_parts(parquet) -> dict:
    """peak_id -> per-peak sub-DataFrame. Call once; per-peak lookups are O(1)."""
    return pl.read_parquet(parquet).partition_by("peak_id", as_dict=True)


def peak_actuation(parts: dict, peak: str) -> dict:
    """sample -> (actuation, coverage) for one peak."""
    sub = parts.get((peak,))
    if sub is None:
        return {}
    return {
        s: (a, c)
        for s, a, c in zip(
            sub["sample"].to_list(),
            sub["actuation"].to_list(),
            sub["coverage"].to_list(),
        )
    }
