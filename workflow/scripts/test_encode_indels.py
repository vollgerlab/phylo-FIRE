#!/usr/bin/env python3
"""Self-checks for encode_indels simple-coding logic. Run: python test_encode_indels.py"""
from encode_indels import simple_indel_code


def test_internal_gap_coded():
    # Internal gap [4,7) shared by t1,t2; t3 has sequence there -> column [1,1,0].
    seqs = {"t1": "ACGT---ACGT", "t2": "ACGT---ACGT", "t3": "ACGTAAAACGT"}
    _, codes = simple_indel_code(seqs)
    assert codes["t1"] == "1" and codes["t2"] == "1" and codes["t3"] == "0", codes


def test_terminal_gap_not_coded():
    # A shared TRAILING gap is a window/length artifact, not an indel event.
    seqs = {"t1": "ACGTACGT---", "t2": "ACGTACGT---", "t3": "ACGTACGTACG"}
    _, codes = simple_indel_code(seqs)
    assert all(c == "" for c in codes.values()), codes


def test_invariant_column_dropped():
    # Internal gap present in ALL taxa -> all '1', no '0' -> not variable -> dropped.
    seqs = {"t1": "ACGT---ACGT", "t2": "ACGT---ACGT", "t3": "ACGT---ACGT"}
    _, codes = simple_indel_code(seqs)
    assert all(c == "" for c in codes.values()), codes


def test_outside_span_is_missing():
    # t3's sequence starts after the gap -> the gap is in its unobserved 5'
    # region -> '?', not '0'. t4 spans the gap (state 0) so the column stays
    # variable and is not dropped.
    seqs = {
        "t1": "ACGT---ACGT", "t2": "ACGT---ACGT",
        "t3": "-------ACGT", "t4": "ACGTAAAACGT",
    }
    _, codes = simple_indel_code(seqs)
    assert codes["t3"] == "?", codes
    assert codes["t4"] == "0", codes


if __name__ == "__main__":
    test_internal_gap_coded()
    test_terminal_gap_not_coded()
    test_invariant_column_dropped()
    test_outside_span_is_missing()
    print("ok")
