"""Tests for aggregate_hg38._alt."""
from aggregate_hg38 import _alt


def test_alt():
    # biallelic: ALT = the one non-ref allele
    assert _alt("G", "G", "C") == "C"
    assert _alt("G", "G", "C", "C") == "C"          # driver allele = the ALT, dedup
    # ref-driver: driver == REF must NOT be added to ALT
    assert _alt("G", "G", "C", "G") == "C"
    # >2-allele column: the driver's 3rd base must appear in ALT
    assert _alt("C", "C", "A", "T") == "A,T"
    # gaps / N / NA are never alleles
    assert _alt("A", "A", "-", "NA") == "."
    assert _alt("A", "A", "A") == "."               # monomorphic vs ref
