#!/usr/bin/env python3
"""Project every variable site in a peak's alignment to CHM13 + hg38 (+ each
carrier's assembly), flagging which sites are attribution drivers.

Projects onto the Stage-1 flanked alignment (aligned/{peak}.fa) so driver columns
match phylo_mm's `linked_sites` exactly, including homoplastic and flank drivers.
CHM13/hg38 are added with `mafft --add`, which only inserts all-hap-gap columns,
so hap columns are preserved; alignment columns then map to genomic coordinates.

A driver is emitted as a bed/VCF SNV when at least one of its `linked_sites`
columns (or a column whose carriers equal its clade) is a real substitution.
Drivers defined only by an indel (a contiguous gap run) can't be a single-base
SNV, so they are returned separately as "unlifted" with a reason, not fabricated
into a position.

Emits one row per variable site (schema `_COLS`) plus a list of unlifted-driver
rows (schema `_UNLIFTED_COLS`).
"""
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass
from itertools import accumulate
from pathlib import Path

from build_re_variants import is_divergent
from common_io import parse_fasta, revcomp
from phylo_mm_io import PhyloMM, attributed_band, read_phylo_mm, selected_drivers

_CHM13_CHR = re.compile(r"^chr([0-9]+|X|Y|M)$")


def cons_agc_contig(cons_chr: str) -> str | None:
    """AGC contig key for a peak's consensus reference. A standard CHM13
    chromosome gets the `CHM13#0#` sample prefix; a peak whose consensus sits on
    an assembly instead of CHM13 is queried on that assembly directly -- a PanSN
    contig (`HG02257#1#JAGYVH...`) passes through, and HG002's T2T
    `chrN_PATERNAL`/`_MATERNAL` maps to `HG002#1#chrN` / `HG002#2#chrN` (hap1=pat,
    hap2=mat, matching add_flanks). None for an unresolvable name (-> cons NA)."""
    if _CHM13_CHR.match(cons_chr):
        return f"CHM13#0#{cons_chr}"
    if "#" in cons_chr:
        return cons_chr
    for suffix, hap in (("_MATERNAL", "2"), ("_PATERNAL", "1")):
        if cons_chr.endswith(suffix):
            return f"HG002#{hap}#{cons_chr[: -len(suffix)]}"
    return None

def _base(ch: str):
    """Upper-cased A/C/G/T, or None for a gap/N/other."""
    ch = ch.upper()
    return ch if ch in ("A", "C", "G", "T") else None


def parse_agc_fasta(text: str) -> str:
    return "".join(ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">"))


def pull_agc(agc: str, contig: str, start0: int, end0: int, retries: int = 4) -> str | None:
    """agc getctg <agc> <contig>:<start>-<end>. agc coords are 0-BASED (the base
    returned for N is the standard 1-based base N+1), so callers pass a 0-based
    start DIRECTLY -- not BED-start+1.

    Retries transient failures with backoff: under high concurrency agc can abort
    (SIGABRT) from memory/contention on the shared archive. A silent None here
    drops a peak's coords (it looks un-liftable), so we retry before giving up --
    but also keep the worker count modest to avoid the aborts in the first place."""
    region = f"{contig}:{start0}-{end0}"
    for attempt in range(retries):
        try:
            res = subprocess.run(
                ["agc", "getctg", agc, region, "-p"],
                capture_output=True, text=True, check=True,
            )
            seq = parse_agc_fasta(res.stdout).upper()
            if seq:
                return seq
        except FileNotFoundError:
            print("[project_variants] agc not found on PATH", file=sys.stderr)
            return None
        except subprocess.CalledProcessError as ex:
            # returncode < 0 => killed by a signal (SIGABRT under archive contention):
            # transient, worth retrying. returncode > 0 => a clean agc error, e.g. the
            # contig isn't in the archive (asm-only peak whose consensus is on a
            # sample contig, not CHM13/GRCh38): deterministic, so fail fast.
            if ex.returncode is None or ex.returncode >= 0:
                return None
        if attempt < retries - 1:
            time.sleep(0.5 * (attempt + 1))       # 0.5s, 1s, 1.5s backoff
    print(f"[project_variants] agc getctg SIGABRT-failed for {region} after {retries} "
          f"tries (concurrency too high -> this peak lost coords)", file=sys.stderr)
    return None


def add_refs(aln_haps: dict[str, str], refs: dict[str, str], threads: int) -> dict[str, str]:
    """Add ungapped `refs` onto the fixed `aln_haps` alignment via `mafft --add`
    (no --keeplength: ref insertions get their own all-hap-gap columns, so no ref
    base is dropped and every hap column is preserved). Strips `_R_` prefixes."""
    if not refs:
        return dict(aln_haps)
    with tempfile.NamedTemporaryFile("w", suffix=".aln.fa", delete=False) as ah:
        for n, s in aln_haps.items():
            ah.write(f">{n}\n{s}\n")
        base = ah.name
    with tempfile.NamedTemporaryFile("w", suffix=".ref.fa", delete=False) as rf:
        for n, s in refs.items():
            rf.write(f">{n}\n{s}\n")
        addf = rf.name
    try:
        res = subprocess.run(
            ["mafft", "--add", addf, "--adjustdirection", "--thread", str(threads), base],
            capture_output=True, text=True, check=True,
        )
    finally:
        Path(base).unlink(missing_ok=True)
        Path(addf).unlink(missing_ok=True)
    return parse_fasta(res.stdout, strip_r=True)


def orient_and_map(aln_ref: str, fwd_seq: str, start0: int):
    """(orientation, col->pos0) for a ref's aligned string vs its FORWARD pull
    starting at 0-based `start0`. pos0 = 0-based genomic position at a column."""
    ungapped = aln_ref.replace("-", "").upper()
    end0 = start0 + len(fwd_seq)
    if ungapped == fwd_seq.upper():
        orient = "fwd"
    elif ungapped == revcomp(fwd_seq).upper():
        orient = "rev"
    else:
        fwd_id = sum(a == b for a, b in zip(ungapped, fwd_seq.upper()))
        rev_id = sum(a == b for a, b in zip(ungapped, revcomp(fwd_seq).upper()))
        orient = "fwd" if fwd_id >= rev_id else "rev"
    col2pos: dict[int, int] = {}
    k = 0
    for c, ch in enumerate(aln_ref):
        if ch != "-":
            col2pos[c] = start0 + k if orient == "fwd" else end0 - 1 - k
            k += 1
    return orient, col2pos


def _col_bridge(aln_haps: dict[str, str], aln_full: dict[str, str]) -> dict[int, int]:
    """Map each aln_haps column -> the aln_full column that holds the same
    alignment position. `mafft --add` only inserts all-hap-gap columns, so every
    hap's k-th non-gap base occupies the same relative slot in both; we walk both
    by ungapped offset (unambiguous, unlike counting all-gap columns)."""
    names = [n for n in aln_haps if n in aln_full]
    full_pos = {n: [i for i, ch in enumerate(aln_full[n]) if ch != "-"] for n in names}
    L = len(next(iter(aln_haps.values())))
    off = {n: 0 for n in names}
    bridge: dict[int, int] = {}
    for c in range(L):
        for n in names:
            if aln_haps[n][c] != "-":
                bridge[c] = full_pos[n][off[n]]
                off[n] += 1
    return bridge


def _col_carriers(aln_haps: dict[str, str], min_carriers: int):
    """Per column: (minor_allele, major_allele, frozenset(minor carriers)) for
    columns that are a real substitution (>=2 alleles). Small-carrier columns are
    kept (a driver may have few carriers); the emit loop applies min_carriers."""
    names = list(aln_haps)
    L = len(next(iter(aln_haps.values())))
    out = {}
    for c in range(L):
        bases = {n: aln_haps[n][c].upper() for n in names}
        cnt = Counter(b for b in bases.values() if b not in ("-", "N"))
        if len(cnt) < 2:
            continue
        # most_common breaks count ties by insertion order and keeps the two
        # alleles DISTINCT; min/max both return the first key on an exact tie
        # (e.g. 39 G / 39 C -> major==minor==G, dropping C -> a bogus X>. call).
        mc = cnt.most_common()
        major, minor = mc[0][0], mc[-1][0]
        carr = frozenset(n for n, b in bases.items() if b == minor)
        out[c] = (minor, major, carr)
    return out


def _driver_cols(colcar, drivers, re_band):
    """Map each driver to the SNV columns that represent it.

    Primary: its `linked_sites` columns that are substitutions (what the plot
    annotates). Fallback, only when none is: an RE-band column whose minor
    carriers equal the driver's clade. `drivers` is {frozenset(desc): (beta,
    [linked_1based])}; `re_band` is 0-based inclusive (lo, hi) or None (all).
    Returns (col->beta, col->clade, unmatched drivers for the unlifted list);
    col->clade pins the reported allele to the clade, not to allele count."""
    lo, hi = re_band if re_band else (0, 10**9)
    by_carrier = {}
    for c, (_, _, carr) in colcar.items():
        if lo <= c <= hi:
            by_carrier.setdefault(carr, []).append(c)
    col_beta: dict[int, str] = {}
    col_desc: dict[int, frozenset] = {}
    unmatched = []
    for desc, (beta, linked) in drivers.items():
        cols = {ls - 1 for ls in linked if (ls - 1) in colcar}   # 1-based -> 0-based
        if not cols:
            cols = set(by_carrier.get(desc, []))                 # homoplasy fallback
        if cols:
            for c in cols:
                col_beta[c] = beta
                col_desc[c] = desc
        else:
            unmatched.append((desc, beta, linked))
    return col_beta, col_desc, unmatched


def _re_band(header: dict[str, str]) -> tuple[int, int] | None:
    """0-based inclusive RE columns from the phylo_mm header, or None if absent.
    int() raises on the 'NA' written for flank=0 or an empty fit, as it always has."""
    lo = int(header["re_first_col"]) if "re_first_col" in header else None
    hi = int(header["re_last_col"]) if "re_last_col" in header else None
    return (lo - 1, hi - 1) if lo is not None and hi is not None else None


def _load_drivers(pm: PhyloMM) -> dict[frozenset, tuple]:
    """frozenset(clade leaves) -> (beta string, [linked_sites]) for is_split rows."""
    return {d.carriers: (d.beta_str, list(d.cols)) for d in selected_drivers(pm)}


def _left_flank(aln_hap_gapped: str, re_seq: str, reverse: bool = False) -> int | None:
    """Length of assembly context before the RE within a hap's flanked sequence
    (0 for RE-only fallback haps), so an alignment column maps to an assembly
    offset. None if the RE can't be located (asm coord then reported NA).

    ``reverse``: MAFFT reverse-complemented this hap (``_R_``), so its aligned
    sequence is revcomp(assembly-strand flanked seq). The RE is then found as its
    reverse complement and the flank is still counted on the assembly strand."""
    if not re_seq:
        return None
    ung = aln_hap_gapped.replace("-", "").upper()
    if reverse:
        i = ung.find(revcomp(re_seq).upper())
        return len(ung) - i - len(re_seq) if i >= 0 else None
    i = ung.find(re_seq.upper())
    return i if i >= 0 else None


@dataclass(frozen=True)
class RefMap:
    """A forward-strand reference pull placed in the alignment by mafft --add."""
    fwd: str                  # forward-strand sequence
    s0: int                   # 0-based genomic start of ``fwd``
    orient: str               # "fwd" or "rev" relative to the alignment
    col2pos: dict[int, int]   # aln_full column -> 0-based genomic position

    def pos0(self, fcol: int | None) -> int | None:
        return self.col2pos.get(fcol) if (self.col2pos and fcol is not None) else None

    def base(self, pos0: int) -> str:
        return self.fwd[pos0 - self.s0].upper()


def _pull_refs(agc, span, cons_chr, cons_start, cons_end, hg38_chr, hg38_start, hg38_end):
    """{"__CHM13__"/"__hg38__": (forward seq, 0-based start)} over the FLANKED
    extent (peak +/- span) so flank columns also map to coordinates. A reference
    that does not lift or fails to pull is left out."""
    pulls = {}
    cons_s0 = max(0, cons_start - span)
    cons_contig = cons_agc_contig(cons_chr)      # CHM13#0#chrN, or the assembly contig
    cons_fwd = pull_agc(agc, cons_contig, cons_s0, cons_end + span) if cons_contig else None
    if cons_fwd:
        pulls["__CHM13__"] = (cons_fwd, cons_s0)
    if hg38_chr != "NA" and str(hg38_start).isdigit():
        hg38_s0 = max(0, int(hg38_start) - span)
        hg38_fwd = pull_agc(agc, f"GRCh38#0#{hg38_chr}", hg38_s0, int(hg38_end) + span)
        if hg38_fwd:
            pulls["__hg38__"] = (hg38_fwd, hg38_s0)
    return pulls


def _ref_maps(aln_haps, pulls, threads):
    """Add the pulled references to the alignment. Returns (bridge from aln_haps
    to aln_full columns, CHM13 RefMap or None, hg38 RefMap or None)."""
    aln_full = add_refs(aln_haps, {n: fwd for n, (fwd, _) in pulls.items()}, threads)
    bridge = _col_bridge(aln_haps, aln_full)
    maps = {}
    for name, (fwd, s0) in pulls.items():
        if name in aln_full:
            orient, col2pos = orient_and_map(aln_full[name], fwd, s0)
            maps[name] = RefMap(fwd, s0, orient, col2pos)
    return bridge, maps.get("__CHM13__"), maps.get("__hg38__")


def _out_is_rev(cons: RefMap | None, hg38: RefMap | None) -> bool:
    """Whether to reverse-complement hap alleles for the output genome.

    The reported REF is the FORWARD base of whichever genome the bed/VCF keys on
    (hg38 when the peak lifts, else CHM13); the cohort major/minor come from the
    haps in the MSA's orientation, so flip them onto THAT genome's forward strand.
    Use hg38's own orientation, NOT CHM13's -- CHM13 and GRCh38 sit on opposite
    strands at some loci (a reference inversion, e.g. chr1_145862159), and keying
    the flip on CHM13 there put the forward-hg38 REF next to reverse-strand ALTs."""
    if hg38 is not None and hg38.col2pos:
        return hg38.orient == "rev"
    return cons is not None and cons.orient == "rev"


def _ungapped_offsets(seq: str) -> list[int]:
    """offsets[c] = number of non-gap bases before column c (one O(L) pass)."""
    return list(accumulate((ch != "-" for ch in seq), initial=0))


def _driver_allele(aln_haps, col, clade, out_is_rev):
    """(reported allele, carriers) for a driver column: the phylo_mm clade's base,
    which is INDEPENDENT of frequency (it can be the major allele) and of the
    reference (it can be the REF base). Carriers become every hap with that base,
    so beta and carriers describe the SAME haplotypes. None if no clade member
    has a base at ``col``."""
    hap_names = list(aln_haps)
    clade = clade & set(hap_names)
    cb = Counter(b for b in (_base(aln_haps[n][col]) for n in clade) if b)
    if not cb:
        return None
    drv = cb.most_common(1)[0][0]                      # clade allele (aln)
    carriers = frozenset(n for n in hap_names if _base(aln_haps[n][col]) == drv)
    return (revcomp(drv) if out_is_rev else drv), carriers


def _asm_hits(carriers, asm, left_flank, offsets, col, reversed_haps=frozenset()) -> str:
    """``sample=chr:pos`` assembly coordinates of ``col`` for each carrier.
    A hap in ``reversed_haps`` is reverse-complemented in the MSA, so its k-th
    aligned base is base (n - 1 - k) of its assembly-strand sequence."""
    hits = []
    for nm in sorted(carriers):
        a = asm.get(nm)
        lf = left_flank.get(nm)
        if a and lf is not None and str(a["start"]).isdigit():
            k = offsets[nm][col]
            if nm in reversed_haps:
                k = offsets[nm][-1] - 1 - k
            hits.append(f"{nm}={a['chr']}:{int(a['start']) + k - lf}")
    return ";".join(hits) if hits else "NA"


def project_peak(peak_id, aln_haps, re_seqs, asm, phylo_mm, cons_chr, cons_start, cons_end,
                 hg38_chr="NA", hg38_start="NA", hg38_end="NA", agc="",
                 flank=0, threads=2, pad=25, min_carriers=2, reversed_haps=frozenset()):
    """Project a peak's variable sites (from its Stage-1 flanked alignment) to
    CHM13 + hg38 + assembly, flagging drivers by their phylo_mm column identity.

    ``aln_haps``: {sample -> gapped aligned seq} (aligned/{peak}.fa);
    ``re_seqs``: {sample -> RE-only seq} (for the assembly-offset lookup);
    ``asm``: {sample -> {chr,start,end}}; ``flank``: bp of context each side;
    ``reversed_haps``: samples MAFFT reverse-complemented (``_R_``, already
    stripped from the ``aln_haps`` keys), for the assembly-coordinate math.
    Returns ``(rows, unlifted)``. ``rows`` is the variable-site schema (``_COLS``);
    ``unlifted`` is indel-only drivers that can't be a single-base SNV."""
    if len(aln_haps) < 2:
        return [], []
    pm = read_phylo_mm(phylo_mm)
    drivers = _load_drivers(pm)
    hap_names = list(aln_haps)
    L = len(next(iter(aln_haps.values())))

    colcar = _col_carriers(aln_haps, min_carriers)
    re_band = _re_band(pm.header)
    region_mode = pm.header.get("attribution_mode", "") == "region"
    # Same band build_re_variants attributes over: whole alignment in region mode.
    band_lo, band_hi = attributed_band(pm.header, L)
    divergent = {desc: is_divergent(len(linked), band_hi - band_lo + 1)
                 for desc, (_, linked) in drivers.items()}
    # Driver carrier-set fallback: whole alignment in region mode, else RE band.
    col_beta, col_desc, unmatched = _driver_cols(
        colcar, drivers, None if region_mode else re_band)

    pulls = _pull_refs(agc, int(flank) + pad, cons_chr, cons_start, cons_end,
                       hg38_chr, hg38_start, hg38_end)
    bridge, cons, hg38 = _ref_maps(aln_haps, pulls, threads)
    out_is_rev = _out_is_rev(cons, hg38)
    left_flank = {n: _left_flank(aln_haps[n], re_seqs.get(n, ""), n in reversed_haps)
                  for n in hap_names}
    offsets = {n: _ungapped_offsets(s) for n, s in aln_haps.items()}

    rows = []
    for col in range(L):
        if col not in colcar:
            continue
        minor, major, carriers = colcar[col]      # major/minor are FREQUENCY labels
        is_driver = col in col_beta
        # Emit non-driver variants only within the regulatory element; the flank is
        # assembly context added to resolve the tree, not RE variation. Drivers are
        # always emitted (incl. region-mode drivers that sit in the flank).
        if not is_driver and re_band is not None and not (re_band[0] <= col <= re_band[1]):
            continue
        if not is_driver and min(len(carriers), len(hap_names) - len(carriers)) < min_carriers:
            continue

        driver_allele = "NA"
        if is_driver and col in col_desc:
            hit = _driver_allele(aln_haps, col, col_desc[col], out_is_rev)
            if hit:
                driver_allele, carriers = hit

        fcol = bridge.get(col)
        cons_pos0 = cons.pos0(fcol) if cons else None
        hg38_pos0 = hg38.pos0(fcol) if hg38 else None
        maj_f, min_f = (revcomp(major), revcomp(minor)) if out_is_rev else (major, minor)
        rows.append({
            "peak_id": peak_id,
            "re_pos": (cons_pos0 - cons_start + 1) if cons_pos0 is not None else "NA",
            "cons_chr": cons_chr,
            "cons_pos": (cons_pos0 + 1) if cons_pos0 is not None else "NA",
            "cons_ref": cons.base(cons_pos0) if cons_pos0 is not None else "NA",
            "hg38_chr": hg38_chr if hg38_pos0 is not None else "NA",
            "hg38_pos": (hg38_pos0 + 1) if hg38_pos0 is not None else "NA",
            "hg38_ref": hg38.base(hg38_pos0) if hg38_pos0 is not None else "NA",
            "hap_major": maj_f,
            "hap_minor": min_f,
            "driver_allele": driver_allele,
            "n_carriers": len(carriers),
            "is_driver": "TRUE" if is_driver else "FALSE",
            "beta": col_beta.get(col, "NA"),
            "carriers": ",".join(sorted(carriers)),
            "asm_coords": _asm_hits(carriers, asm, left_flank, offsets, col, reversed_haps),
            "divergent": "TRUE" if is_driver and divergent.get(col_desc.get(col)) else "FALSE",
        })

    unlifted = [{
        "peak_id": peak_id,
        "n_carriers": len(desc),
        "beta": beta,
        "reason": "indel_only_no_snv_column",
        "carriers": ",".join(sorted(desc)),
        "divergent": "TRUE" if divergent.get(desc) else "FALSE",
    } for desc, beta, _ in unmatched]
    return rows, unlifted


_COLS = ["peak_id", "re_pos", "cons_chr", "cons_pos", "cons_ref",
         "hg38_chr", "hg38_pos", "hg38_ref", "hap_major", "hap_minor",
         "driver_allele", "n_carriers", "is_driver", "beta", "carriers",
         "asm_coords", "divergent"]
_UNLIFTED_COLS = ["peak_id", "n_carriers", "beta", "reason", "carriers", "divergent"]


def _header() -> str:
    return "\t".join(_COLS) + "\n"


def _unlifted_header() -> str:
    return "\t".join(_UNLIFTED_COLS) + "\n"
