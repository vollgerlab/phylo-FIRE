#!/usr/bin/env python3
"""Self-checks for project_variants' column bridge + driver flagging.

    python workflow/scripts/test_project_variants.py
"""
from project_variants import _col_bridge, _col_carriers, _driver_cols, cons_agc_contig


def test_cons_agc_contig():
    # CHM13 chromosomes get the CHM13#0# prefix; a consensus on an assembly is
    # queried on that assembly directly (PanSN passthrough / HG002 T2T mapping).
    assert cons_agc_contig("chr6") == "CHM13#0#chr6"
    assert cons_agc_contig("chrX") == "CHM13#0#chrX"
    assert cons_agc_contig("HG02257#1#JAGYVH020000072.1") == "HG02257#1#JAGYVH020000072.1"
    assert cons_agc_contig("chr6_PATERNAL") == "HG002#1#chr6"
    assert cons_agc_contig("chr12_MATERNAL") == "HG002#2#chr12"
    assert cons_agc_contig("weird_name") is None


def test_col_bridge():
    # aln_full = aln_haps with one inserted all-hap-gap column (a ref insertion)
    # between haps-col 1 and 2. Each haps col must map to the full col holding the
    # same base, computed purely from ungapped offsets.
    aln_haps = {"h1": "ACGT", "h2": "A-GT"}
    aln_full = {"h1": "AC-GT", "h2": "A--GT", "__ref__": "ACXGT"}
    b = _col_bridge(aln_haps, aln_full)
    assert b == {0: 0, 1: 1, 2: 3, 3: 4}, b


def test_col_carriers_exact_tie():
    # 2 G / 2 C is an exact count tie: major and minor must stay DISTINCT
    # (the old min/max both returned G, dropping C -> a bogus X>. projection).
    aln = {"a": "G", "b": "G", "c": "C", "d": "C"}
    (minor, major, carr) = _col_carriers(aln, min_carriers=1)[0]
    assert {major, minor} == {"G", "C"}, (major, minor)
    assert carr == frozenset({n for n in aln if aln[n] == minor})


def test_driver_flag_by_linked_site():
    # col 1 is the only substitution (h2/h3 carry T). A driver whose linked_sites
    # names column 2 (1-based -> col 1) must flag col 1, not the invariant cols.
    aln = {"a": "AAAA", "b": "ATAA", "c": "ATAA"}
    colcar = _col_carriers(aln, min_carriers=1)
    assert set(colcar) == {1}, colcar
    drivers = {frozenset({"b", "c"}): ("0.5", [2])}   # linked_sites 1-based
    col_beta, col_desc, unmatched = _driver_cols(colcar, drivers, re_band=None)
    assert col_beta == {1: "0.5"}, col_beta
    assert col_desc == {1: frozenset({"b", "c"})}, col_desc
    assert unmatched == [], unmatched


def test_indel_only_driver_is_unmatched():
    # A driver whose only linked site is an all-gap (indel) column and whose clade
    # matches no substitution column must be returned unmatched (-> unlifted).
    aln = {"a": "AC-A", "b": "AT-A", "c": "AT-A", "d": "AC-A", "e": "AC-A"}
    colcar = _col_carriers(aln, min_carriers=1)   # col1 carriers = {b,c}; col2 all-gap
    drivers = {frozenset({"a"}): ("0.3", [3])}    # linked col3(1-based)=col2 gap; clade {a} != {b,c}
    col_beta, col_desc, unmatched = _driver_cols(colcar, drivers, re_band=None)
    assert col_beta == {}, col_beta
    assert len(unmatched) == 1, unmatched


def test_driver_allele_is_clade_not_minor():
    # A biallelic driver whose clade carries the MAJORITY allele (G, 3 haps) while
    # the count-minor allele is C (2 haps). major/minor must stay FREQUENCY labels
    # (major=G, minor=C); driver_allele must be the clade's allele (G), and
    # carriers must be the clade -- so beta and carriers describe the same haps.
    from project_variants import project_peak
    import tempfile, os
    aln = {"h1": "G", "h2": "G", "h3": "G", "h4": "C", "h5": "C"}
    fd, pm = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w") as fh:
        fh.write("# attribution_mode=RE\n# re_first_col=1\n# re_last_col=1\n")
        fh.write("is_split\tdescendants\tbeta\tlinked_sites\n")
        fh.write("TRUE\th1,h2,h3\t0.4\t1\n")
    try:
        rows, _ = project_peak("pk", aln, {}, {}, pm, "chrX", 100, 101,
                               agc="", flank=0, min_carriers=1)
    finally:
        os.unlink(pm)
    drv = [r for r in rows if r["is_driver"] == "TRUE"]
    assert len(drv) == 1, rows
    r = drv[0]
    assert (r["hap_major"], r["hap_minor"]) == ("G", "C"), r   # frequency labels
    assert r["driver_allele"] == "G", r                        # clade allele (the major!)
    assert r["n_carriers"] == 3, r
    assert set(r["carriers"].split(",")) == {"h1", "h2", "h3"}, r


if __name__ == "__main__":
    test_cons_agc_contig()
    test_col_bridge()
    test_col_carriers_exact_tie()
    test_driver_flag_by_linked_site()
    test_indel_only_driver_is_unmatched()
    test_driver_allele_is_clade_not_minor()
    print("project_variants self-checks passed")
