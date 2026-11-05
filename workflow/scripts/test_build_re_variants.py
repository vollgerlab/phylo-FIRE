#!/usr/bin/env python3
"""Self-checks for build_re_variants (driver_delta, clean_pair, divergence).

    python workflow/scripts/test_build_re_variants.py
"""
import tempfile
from pathlib import Path

from build_re_variants import build, is_divergent

PHYLO_MM_COLS = [
    "branch_id", "n_below", "descendants", "mean_in", "mean_out", "t_stat",
    "p_raw", "p_fdr", "is_split", "linked_sites", "selected_by_lasso",
    "lasso_coef", "beta", "se", "ci_lower", "ci_upper", "p",
]


def test_driver_delta_and_clean_pair():
    # 3 seqs over a 5bp RE (cols 1..5). col 3 is the driver: carriers have 'T'.
    # driver_delta = sum of carried betas; clean_pair finds the non-carrier that
    # differs only at the driver column.
    with tempfile.TemporaryDirectory() as d:
        aln = Path(d) / "a.fa"
        aln.write_text(">s_carrier\nACTGA\n>s_clean\nACAGA\n>s_dirty\nACAGG\n")
        act = Path(d) / "a.tsv"
        act.write_text("sample\tactuation\tcoverage\n"
                       "s_carrier\t0.1\t10\ns_clean\t0.9\t10\ns_dirty\t0.8\t10\n")
        pm = Path(d) / "pm.tsv"
        pm.write_text(
            "# peak_id=chr1_100_105\n# attribution_mode=RE\n"
            "# re_first_col=1\n# re_last_col=5\n# flank=0\n"
            "# h2_re_mean=0.9\n# h2_re_lower=0.5\n# h2_re_upper=1\n"
            "# h2_region_mean=0.9\n# h2_region_lower=0.5\n# h2_region_upper=1\n"
            "# n_unique_seqs=3\n# n_variants_selected=1\n"
            + "\t".join(PHYLO_MM_COLS) + "\n"
            + "\t".join(["1", "1", "s_carrier", "0.1", "0.85", "-5", "0.01",
                         "0.01", "TRUE", "3", "TRUE", "NA", "-0.5", "0.1",
                         "-0.7", "-0.3", "0.01"]) + "\n"
        )
        rows = build(str(aln), str(act), str(pm))
    carrier = next(r for r in rows if r["has_driver"])
    assert abs(carrier["driver_delta"] - (-0.5)) < 1e-9, carrier["driver_delta"]
    assert carrier["clean_pair"] is True, carrier
    assert carrier["minimal_edit"] == "RE3 A->T", carrier["minimal_edit"]
    assert carrier["n_other_diffs"] == 0
    assert carrier["driver_positions"] == "RE3"
    assert carrier["driver_divergent"] == "TRUE"   # 1 linked site / 5 bp band > 5%


def test_is_divergent():
    assert not is_divergent(5, 200)
    assert is_divergent(11, 200)


if __name__ == "__main__":
    test_driver_delta_and_clean_pair()
    test_is_divergent()
    print("build_re_variants self-checks passed")
