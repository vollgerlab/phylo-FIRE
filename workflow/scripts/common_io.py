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
    reverse-complemented records, so names match the actuation/asm tables and the
    phylo_mm carriers (phylo_mm.R strips it from tree tips). Use
    ``read_alignment`` when the caller also needs to know which were reversed.
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


def read_alignment(path) -> tuple[dict[str, str], frozenset[str]]:
    """Stage-1 haplotype MSA -> ({name: gapped seq}, names MAFFT reversed).

    Names have the ``_R_`` prefix stripped; the second value lists the stripped
    names of the records MAFFT ``--adjustdirection`` reverse-complemented, whose
    aligned sequence is on the opposite strand of the haplotype's assembly."""
    raw = read_fasta(path)
    seqs = {(n[3:] if n.startswith("_R_") else n): s for n, s in raw.items()}
    return seqs, frozenset(n[3:] for n in raw if n.startswith("_R_"))


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
