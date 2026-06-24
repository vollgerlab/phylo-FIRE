#!/usr/bin/env python
"""Trim a flanked MSA down to just the regulatory-element (RE) columns.

The RE tree (for the RE-only heritability test) is built from the SAME mafft
alignment as the whole-region tree, just sliced to the RE band — so the two
H2_p estimates differ only in scope, not in alignment/tree method. The RE band
is the columns spanned once `flank` non-gap bases are trimmed off each end of
every sequence (the same rule phylo_mm's `re_band` uses). With flank=0 (or if
the band can't be found) the alignment is copied through unchanged.

Usage: trim_to_re.py --aln aligned.fa --flank 250 --out re_aligned.fa
"""
import argparse
from pathlib import Path


def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def re_band(seqs: dict[str, str], flank: int) -> tuple[int, int] | None:
    """0-based inclusive [lo, hi] alignment columns of the RE, or None."""
    if flank <= 0 or not seqs:
        return None
    lo, hi = float("inf"), -1
    for s in seqs.values():
        nongap = [i for i, c in enumerate(s) if c != "-"]
        if len(nongap) <= 2 * flank:
            continue
        lo = min(lo, nongap[flank])
        hi = max(hi, nongap[len(nongap) - flank - 1])
    if lo == float("inf") or hi < lo:
        return None
    return int(lo), int(hi)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln", required=True)
    ap.add_argument("--flank", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.aln)
    band = re_band(seqs, args.flank)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as out:
        for name, s in seqs.items():
            sub = s if band is None else s[band[0]:band[1] + 1]
            out.write(f">{name}\n{sub}\n")


if __name__ == "__main__":
    main()
