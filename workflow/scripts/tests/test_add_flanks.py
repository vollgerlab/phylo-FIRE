"""Tests for add_flanks contig mapping."""
from add_flanks import agc_contig


def test_pansn_passthrough():
    assert agc_contig("HG02886_1", "HG02886#1#CM089178.1") == "HG02886#1#CM089178.1"


def test_hg002_maternal_is_hap2():
    # AGC truth: HG002#2#chrN = maternal. asm_chr suffix (the real assembly)
    # drives the mapping, regardless of the BED's flipped Haplotype label.
    assert agc_contig("HG002_1", "chr12_MATERNAL") == "HG002#2#chr12"


def test_hg002_paternal_is_hap1():
    assert agc_contig("HG002_2", "chr12_PATERNAL") == "HG002#1#chr12"


def test_highlight_tagged_name():
    # Highlighted leaves carry a 'tag__' prefix; the individual is still HG002.
    assert agc_contig("highlight__HG002_1", "chr8_MATERNAL") == "HG002#2#chr8"


def test_unknown_nonpansn_is_none():
    assert agc_contig("X_1", "weird_contig_name") is None
