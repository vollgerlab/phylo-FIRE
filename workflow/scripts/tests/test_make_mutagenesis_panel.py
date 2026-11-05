"""Tests for tools/make_mutagenesis_panel (dedup + natural-haplotype annotation)."""
from make_mutagenesis_panel import build


def test_panel_dedups_and_annotates_natural(tmp_path):
    # Driver d1 bundles RE2-4 (perfectly co-occurring in car1/car2). The panel
    # is the background, one single edit per change, and the full carrier.
    aln = tmp_path / "a.fa"
    aln.write_text(">car1\nAGTTAC\n>car2\nAGTTGC\n>non1\nAACAAC\n>non2\nAACAAC\n")
    pm = tmp_path / "pm.tsv"
    pm.write_text(
        "# peak_id=chrT_1_6\n# attribution_mode=RE\n# re_first_col=1\n# re_last_col=6\n"
        "branch_id\tis_split\tlinked_sites\tdescendants\tbeta\n"
        "1\tTRUE\t2,3,4\tcar1,car2\t-0.4\n"
    )
    rv = tmp_path / "rv.tsv"
    # the 'full' (carrier) AGTTAC and the 'background' AACAAC are both natural
    rv.write_text(
        "seq_id\tsamples\tn_haps\tact_mean\tre_seq\n"
        "chrT_1_6_s00\tcar1\t1\t0.10\tAGTTAC\n"
        "chrT_1_6_s01\tnon1,non2\t2\t0.80\tAACAAC\n"
    )
    recs, _ = build(str(aln), str(pm), str(rv))
    seqs = [r["seq"] for r in recs]
    assert len(seqs) == len(set(seqs)), "not unique"
    byseq = {r["seq"]: r for r in recs}
    # full = carrier AGTTAC, annotated natural=s00 AND role d1:full
    full = byseq["AGTTAC"]
    assert "d1:full" in full["roles"] and full["natural"][0] == "chrT_1_6_s00", full
    # background = AACAAC, natural=s01 (2 haps) AND role d1:background
    bg = byseq["AACAAC"]
    assert "d1:background" in bg["roles"] and bg["natural"][0] == "chrT_1_6_s01", bg
    assert bg["natural"][1] == "2", bg["natural"]
    # the single edits are synthetic (no natural match)
    syn = [r for r in recs if r["natural"] is None]
    assert syn
    assert all(any(":RE" in role for role in r["roles"]) for r in syn), syn
