"""Small readers shared by the Stage-2 scripts: FASTA, reverse complement, peak
ids and per-peak actuation TSVs. Stage-1 rule scripts keep their own copies."""

from pathlib import Path

_COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def revcomp(s: str) -> str:
    """Reverse complement of an A/C/G/T/N sequence (either case; gaps kept)."""
    return s.translate(_COMP)[::-1]


def parse_fasta(text: str, strip_r: bool = False) -> dict[str, str]:
    """FASTA text -> {name: seq}; the name is the header up to the first space.

    ``strip_r`` drops the ``_R_`` prefix MAFFT ``--adjustdirection`` puts on
    reverse-complemented records. Callers pass it explicitly: some readers have
    always stripped it and some have not, and changing either changes outputs.
    """
    out: dict[str, str] = {}
    name, parts = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
            if strip_r and name.startswith("_R_"):
                name = name[3:]
        else:
            parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def read_fasta(path, strip_r: bool = False) -> dict[str, str]:
    """FASTA file -> {name: seq}; see ``parse_fasta`` for ``strip_r``."""
    return parse_fasta(Path(path).read_text(), strip_r=strip_r)


def split_peak_id(peak_id: str) -> tuple[str, str, str]:
    """``{chrom}_{start}_{end}`` -> (chrom, start, end) as strings. The chrom may
    itself contain underscores (``chr8_MATERNAL``), so split off the last two."""
    chrom, _, end = peak_id.rpartition("_")
    chrom, _, start = chrom.rpartition("_")
    return chrom, start, end


def _float(x: str) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load_actuation_tsv(path) -> dict[str, tuple[float, float]]:
    """sample -> (actuation, coverage) from a TSV with sample/actuation/coverage
    columns. Rows with a non-numeric value or coverage <= 0 are left out."""
    act: dict[str, tuple[float, float]] = {}
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        si, ai, ci = cols.index("sample"), cols.index("actuation"), cols.index("coverage")
        for line in fh:
            p = line.rstrip("\n").split("\t")
            a, c = _float(p[ai]), _float(p[ci])
            if a is not None and c is not None and c > 0:
                act[p[si]] = (a, c)
    return act
