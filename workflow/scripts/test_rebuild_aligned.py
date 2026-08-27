#!/usr/bin/env python3
"""Self-check for rebuild_aligned's skip-integrity check.

    python workflow/scripts/test_rebuild_aligned.py
"""
import tempfile
from pathlib import Path

from rebuild_aligned import _valid_aligned


def test_valid_aligned():
    d = Path(tempfile.mkdtemp())
    (d / "in.fa").write_text(">a\nACGT\n>b\nACGT\n>c\nAGGT\n")   # 3 input haplotypes
    cases = {
        "valid": (">a\nAC-GT\n>b\nAC-GT\n>c\nAGGGT\n", True),
        "truncated (2 of 3 records)": (">a\nAC-GT\n>b\nAC-GT\n", False),
        "ragged (unequal lengths)": (">a\nAC-GT\n>b\nACGT\n>c\nAGGGT\n", False),
        "empty": ("", False),
        "headers only (zero-len)": (">a\n>b\n>c\n", False),
        "single record": (">a\nACGT\n", False),
    }
    for name, (content, expected) in cases.items():
        (d / "out.fa").write_text(content)
        assert _valid_aligned(d / "out.fa", d / "in.fa") is expected, name


if __name__ == "__main__":
    test_valid_aligned()
    print("rebuild_aligned self-checks passed")
