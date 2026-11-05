"""Reader for the per-peak phylo_mm TSV written by phylo_mm.R.

The file is ``# key=value`` header lines (H², attribution mode, RE band, ...)
then a tab-separated table with one row per candidate variant block. Rows with
``is_split`` TRUE are the selected drivers. Values stay raw strings here, so a
caller can write ``beta`` back out verbatim.
"""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PhyloMM:
    header: dict[str, str]
    columns: list[str]
    records: list[list[str]]  # tab-split data rows, blank lines dropped

    def rows(self) -> list[dict[str, str]]:
        """Data rows as {column: value}; a short row lacks its trailing keys."""
        return [dict(zip(self.columns, f)) for f in self.records]


@dataclass(frozen=True)
class Driver:
    branch_id: str
    cols: tuple[int, ...]  # linked_sites: 1-based alignment columns
    carriers: frozenset[str]  # descendants
    beta_str: str  # joint beta exactly as written ("NA" when absent)
    p_str: str

    @property
    def id(self) -> str:
        return f"d{self.branch_id}"

    @property
    def beta(self) -> float | None:
        return to_float(self.beta_str)

    @property
    def p(self) -> float | None:
        return to_float(self.p_str)


def to_float(x) -> float | None:
    """float(x), or None for NA / blank / missing."""
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_phylo_mm(path) -> PhyloMM:
    """Parse a phylo_mm TSV in one pass. A missing or unreadable file is empty."""
    header: dict[str, str] = {}
    columns: list[str] = []
    records: list[list[str]] = []
    try:
        text = Path(path).read_text() if path else ""
    except OSError:
        text = ""
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            key, sep, val = line[1:].partition("=")
            if sep:
                header[key.strip()] = val.strip()
        elif not columns:
            columns = line.split("\t")
        else:
            records.append(line.split("\t"))
    return PhyloMM(header, columns, records)


def is_true(value: str) -> bool:
    """phylo_mm booleans: R writes TRUE/FALSE; 1 is accepted too."""
    return value.upper() in ("TRUE", "1")


def selected_drivers(pm: PhyloMM) -> list[Driver]:
    """The is_split rows in file order. Rows without a descendants field are
    skipped, and linked_sites keeps only its integer entries."""
    out = []
    for row in pm.rows():
        if not is_true(row.get("is_split", "")) or "descendants" not in row:
            continue
        out.append(Driver(
            branch_id=row.get("branch_id", "?"),
            cols=tuple(int(x) for x in row.get("linked_sites", "").split(",")
                       if x.strip().isdigit()),
            carriers=frozenset(d for d in row["descendants"].split(",") if d),
            beta_str=row.get("beta", "NA"),
            p_str=row.get("p", "NA"),
        ))
    return out


def _col(header: dict[str, str], key: str, default: int) -> int:
    v = header.get(key, "NA")
    return int(v) if v not in ("NA", "") else default


def re_cols(header: dict[str, str], aln_len: int) -> tuple[int, int]:
    """1-based inclusive RE columns; an NA side (flank=0) is the alignment edge."""
    return _col(header, "re_first_col", 1), _col(header, "re_last_col", aln_len)


def attributed_band(header: dict[str, str], aln_len: int) -> tuple[int, int]:
    """1-based inclusive columns the drivers are attributed over: the whole
    alignment in region mode (flank drivers), else the RE band. An empty band
    falls back to the whole alignment."""
    if header.get("attribution_mode", "none") == "region":
        return 1, aln_len
    lo, hi = re_cols(header, aln_len)
    if hi < lo or hi == 0:
        return 1, aln_len
    return lo, hi
