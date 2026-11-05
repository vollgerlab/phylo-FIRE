#!/usr/bin/env python3
"""Project every variable site in a peak's alignment to CHM13 + hg38 (+ each
carrier's assembly), flagging which sites are attribution drivers.

Projects onto the PIPELINE'S OWN Stage-1 flanked alignment (aligned/{peak}.fa),
NOT a fresh re-alignment. This matters: phylo_mm's driver `linked_sites` and
`descendants` are defined on that alignment's columns, so reusing it lets us flag
a driver at the SAME column the tree/MSA plot annotates -- including homoplastic
drivers (whose carriers are NOT a clean column partition) and flank drivers,
which a fresh RE-only re-alignment silently dropped. The CHM13/hg38 reference
sequences are added onto that fixed alignment with `mafft --add` (which only
inserts all-hap-gap columns, so hap columns -- and thus phylo_mm's -- are
preserved), then alignment columns map to genomic coordinates.

A driver is emitted as a bed/VCF SNV when at least one of its `linked_sites`
columns (or a column whose carriers equal its clade) is a real substitution.
Drivers defined only by an indel (a contiguous gap run) can't be a single-base
SNV, so they are returned separately as "unlifted" with a reason, not fabricated
into a position.

Emits one row per variable site (schema `_COLS`) plus a list of unlifted-driver
rows (schema `_UNLIFTED_COLS`).
"""
import argparse
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from build_re_variants import is_divergent

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

_COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _base(ch: str):
    """Upper-cased A/C/G/T, or None for a gap/N/other."""
    ch = ch.upper()
    return ch if ch in ("A", "C", "G", "T") else None


def read_fasta(path: Path) -> dict[str, str]:
    return read_fasta_str(Path(path).read_text())


def read_fasta_str(text: str) -> dict[str, str]:
    out, name, parts = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(parts)
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def parse_agc_fasta(text: str) -> str:
    return "".join(l.strip() for l in text.splitlines() if l and not l.startswith(">"))


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
    aln = read_fasta_str(res.stdout)
    return {n[3:] if n.startswith("_R_") else n: s for n, s in aln.items()}


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
    """For each driver, the SNV columns that represent it. Primary: its
    `linked_sites` columns that are substitutions -- exactly what the tree/MSA
    plot annotates, so the bed can't disagree with the plot. Fallback (only when
    NONE of a driver's linked columns is a substitution -- e.g. a mono column
    from a 1-col alignment nudge): a column whose minor-carriers equal its clade,
    restricted to the RE band so a perfectly-linked flank hitchhiker isn't
    mislabelled the causal variant. Returns (col->beta, list-of-unmatched).

    `drivers`: {frozenset(desc): (beta, [linked_1based])}; `re_band`: (lo,hi)
    0-based inclusive RE columns, or None to allow all. Unmatched drivers
    (indel-only / unresolvable) are returned for the unlifted list. `col_desc`
    carries each driver column's clade (descendants) so the emit loop can pin the
    reported driver allele + carriers to the phylo_mm clade, not to allele count
    (which ties on balanced sites and can even disagree with the clade)."""
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


def _re_header(phylo_mm):
    """((lo, hi) 0-based inclusive RE columns or None, attribution_mode string).
    `linked_sites` and re_first/last_col are 1-based columns -> subtract 1."""
    band, mode = None, ""
    if phylo_mm and Path(phylo_mm).exists():
        lo = hi = None
        for line in Path(phylo_mm).read_text().splitlines():
            if not line.startswith("#"):
                break
            if line.startswith("# attribution_mode="):
                mode = line.split("=", 1)[1].strip()
            elif line.startswith("# re_first_col="):
                lo = int(line.split("=", 1)[1])
            elif line.startswith("# re_last_col="):
                hi = int(line.split("=", 1)[1])
        if lo is not None and hi is not None:
            band = (lo - 1, hi - 1)
    return band, mode


def _load_drivers(phylo_mm) -> dict[frozenset, tuple]:
    """frozenset(clade leaves) -> (beta, [linked_sites]) for is_split rows."""
    out = {}
    if not phylo_mm or not Path(phylo_mm).exists():
        return out
    lines = [l for l in Path(phylo_mm).read_text().splitlines() if not l.startswith("#")]
    if not lines:
        return out
    idx = {h: i for i, h in enumerate(lines[0].split("\t"))}
    if "is_split" not in idx or "descendants" not in idx:
        return out
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= idx["descendants"] or f[idx["is_split"]].upper() not in ("TRUE", "1"):
            continue
        desc = frozenset(d for d in f[idx["descendants"]].split(",") if d)
        beta = f[idx["beta"]] if "beta" in idx and len(f) > idx["beta"] else "NA"
        linked = []
        if "linked_sites" in idx and len(f) > idx["linked_sites"]:
            linked = [int(x) for x in f[idx["linked_sites"]].split(",") if x.strip().isdigit()]
        out[desc] = (beta, linked)
    return out


def _left_flank(aln_hap_gapped: str, re_seq: str) -> int | None:
    """Length of assembly context before the RE within a hap's flanked sequence
    (0 for RE-only fallback haps), so an alignment column maps to an assembly
    offset. None if the RE can't be located (asm coord then reported NA)."""
    if not re_seq:
        return None
    ung = aln_hap_gapped.replace("-", "").upper()
    i = ung.find(re_seq.upper())
    return i if i >= 0 else None


def project_peak(peak_id, aln_haps, re_seqs, asm, phylo_mm, cons_chr, cons_start, cons_end,
                 hg38_chr="NA", hg38_start="NA", hg38_end="NA", agc="",
                 flank=0, threads=2, pad=25, min_carriers=2):
    """Project a peak's variable sites (from its Stage-1 flanked alignment) to
    CHM13 + hg38 + assembly, flagging drivers by their phylo_mm column identity.

    ``aln_haps``: {sample -> gapped aligned seq} (aligned/{peak}.fa);
    ``re_seqs``: {sample -> RE-only seq} (for the assembly-offset lookup);
    ``asm``: {sample -> {chr,start,end}}; ``flank``: bp of context each side.
    Returns ``(rows, unlifted)``. ``rows`` is the variable-site schema (``_COLS``);
    ``unlifted`` is indel-only drivers that can't be a single-base SNV."""
    if len(aln_haps) < 2:
        return [], []
    drivers = _load_drivers(phylo_mm)
    hap_names = list(aln_haps)
    L = len(next(iter(aln_haps.values())))

    colcar = _col_carriers(aln_haps, min_carriers)
    re_band, attr_mode = _re_header(phylo_mm)
    # Same band build_re_variants attributes over: whole alignment in region mode.
    band_len = L if (attr_mode == "region" or re_band is None) else re_band[1] - re_band[0] + 1
    divergent = {desc: is_divergent(len(linked), band_len)
                 for desc, (_, linked) in drivers.items()}
    # Driver carrier-set fallback: whole alignment in region mode, else RE band.
    col_beta, col_desc, unmatched = _driver_cols(
        colcar, drivers, None if attr_mode == "region" else re_band)

    # Reference pulls span the FLANKED extent (peak +/- flank +/- pad) so flank
    # columns also map to coordinates.
    span = int(flank) + pad
    refs: dict[str, str] = {}
    cons_s0 = max(0, cons_start - span)
    cons_contig = cons_agc_contig(cons_chr)      # CHM13#0#chrN, or the assembly contig
    cons_fwd = pull_agc(agc, cons_contig, cons_s0, cons_end + span) if cons_contig else None
    if cons_fwd:
        refs["__CHM13__"] = cons_fwd
    hg38_ok = hg38_chr != "NA" and str(hg38_start).isdigit()
    hg38_s0 = hg38_fwd = None
    if hg38_ok:
        hg38_s0 = max(0, int(hg38_start) - span)
        hg38_fwd = pull_agc(agc, f"GRCh38#0#{hg38_chr}", hg38_s0, int(hg38_end) + span)
        if hg38_fwd:
            refs["__hg38__"] = hg38_fwd

    aln_full = add_refs(aln_haps, refs, threads)
    bridge = _col_bridge(aln_haps, aln_full)
    cons_map, cons_orient = None, "fwd"
    if "__CHM13__" in aln_full and cons_fwd:
        cons_orient, cons_map = orient_and_map(aln_full["__CHM13__"], cons_fwd, cons_s0)
    hg38_map, hg38_orient = None, "fwd"
    if "__hg38__" in aln_full and hg38_fwd and hg38_s0 is not None:
        hg38_orient, hg38_map = orient_and_map(aln_full["__hg38__"], hg38_fwd, hg38_s0)
    # The reported REF is the FORWARD base of whichever genome the bed/VCF keys on
    # (hg38 when the peak lifts, else CHM13); the cohort major/minor come from the
    # haps in the MSA's orientation, so flip them onto THAT genome's forward strand.
    # Use hg38's own orientation, NOT CHM13's -- CHM13 and GRCh38 sit on opposite
    # strands at some loci (a reference inversion, e.g. chr1_145862159), and keying
    # the flip on CHM13 there put the forward-hg38 REF next to reverse-strand ALTs.
    out_is_rev = (hg38_orient == "rev") if (hg38_fwd and hg38_map) else (cons_orient == "rev")

    left_flank = {n: _left_flank(aln_haps[n], re_seqs.get(n, "")) for n in hap_names}

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
        n_car = len(carriers)
        n_ref = len(hap_names) - n_car
        if not is_driver and min(n_car, n_ref) < min_carriers:
            continue

        fcol = bridge.get(col)
        cons_pos0 = cons_map.get(fcol) if (cons_map and fcol is not None) else None
        hg38_pos0 = hg38_map.get(fcol) if (hg38_map and fcol is not None) else None
        cons_ref = cons_fwd[cons_pos0 - cons_s0].upper() if cons_pos0 is not None else "NA"
        hg38_ref = (hg38_fwd[hg38_pos0 - hg38_s0].upper()
                    if (hg38_fwd and hg38_pos0 is not None) else "NA")
        maj_f, min_f = (revcomp(major), revcomp(minor)) if out_is_rev else (major, minor)

        # The driver/effect allele is the phylo_mm clade's allele -- reported
        # explicitly because it is INDEPENDENT of frequency (it can be the major
        # allele) and of the reference (it can be the REF base). For a driver we
        # also set carriers/n_carriers to that clade, so beta and carriers describe
        # the SAME haplotypes; hap_major/hap_minor stay pure frequency labels.
        driver_allele = "NA"
        if is_driver and col in col_desc:
            clade = col_desc[col] & set(hap_names)
            cb = Counter(b for b in (_base(aln_haps[n][col]) for n in clade) if b)
            if cb:
                drv = cb.most_common(1)[0][0]                      # clade allele (aln)
                driver_allele = revcomp(drv) if out_is_rev else drv
                carriers = frozenset(
                    n for n in hap_names if _base(aln_haps[n][col]) == drv)
                n_car = len(carriers)

        asm_hits = []
        for nm in sorted(carriers):
            a = asm.get(nm)
            lf = left_flank.get(nm)
            if a and lf is not None and str(a["start"]).isdigit():
                # ungapped offset of this hap at `col`
                off = sum(1 for ch in aln_haps[nm][:col] if ch != "-")
                asm_hits.append(f"{nm}={a['chr']}:{int(a['start']) + off - lf}")

        rows.append({
            "peak_id": peak_id,
            "re_pos": (cons_pos0 - cons_start + 1) if cons_pos0 is not None else "NA",
            "cons_chr": cons_chr,
            "cons_pos": (cons_pos0 + 1) if cons_pos0 is not None else "NA",
            "cons_ref": cons_ref,
            "hg38_chr": hg38_chr if hg38_pos0 is not None else "NA",
            "hg38_pos": (hg38_pos0 + 1) if hg38_pos0 is not None else "NA",
            "hg38_ref": hg38_ref,
            "hap_major": maj_f,
            "hap_minor": min_f,
            "driver_allele": driver_allele,
            "n_carriers": n_car,
            "is_driver": "TRUE" if is_driver else "FALSE",
            "beta": col_beta.get(col, "NA"),
            "carriers": ",".join(sorted(carriers)),
            "asm_coords": ";".join(asm_hits) if asm_hits else "NA",
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


def _lookup_coords(table: Path, peak_id: str) -> dict | None:
    lines = Path(table).read_text().splitlines()
    if not lines:
        return None
    hdr = lines[0].split("\t")
    for ln in lines[1:]:
        f = ln.split("\t")
        if f and f[0] == peak_id:
            return dict(zip(hdr, f))
    return None


def _load_asm(tsv: Path) -> dict[str, dict]:
    out = {}
    lines = Path(tsv).read_text().splitlines()
    if not lines:
        return out
    idx = {h: i for i, h in enumerate(lines[0].split("\t"))}
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= max(idx.get("asm_end", 0), idx.get("sample", 0)):
            continue
        out[f[idx["sample"]]] = {
            "chr": f[idx["asm_chr"]], "start": f[idx["asm_start"]], "end": f[idx["asm_end"]],
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--aligned", required=True, type=Path, help="aligned/{peak}.fa (flanked MSA)")
    ap.add_argument("--fa", required=True, type=Path, help="fastas/{peak}.fa (RE-only, for asm offset)")
    ap.add_argument("--tsv", required=True, type=Path, help="fastas/{peak}.tsv (asm coords)")
    ap.add_argument("--phylo-mm", type=Path, default=None, help="phylo_mm/{peak}.tsv")
    ap.add_argument("--peak-id", required=True)
    ap.add_argument("--coords-table", type=Path, default=None)
    ap.add_argument("--cons-chr")
    ap.add_argument("--cons-start", type=int, help="0-based (BED)")
    ap.add_argument("--cons-end", type=int)
    ap.add_argument("--hg38-chr", default="NA")
    ap.add_argument("--hg38-start", default="NA")
    ap.add_argument("--hg38-end", default="NA")
    ap.add_argument("--agc", required=True)
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--min-carriers", type=int, default=2)
    ap.add_argument("--pad", type=int, default=25)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-unlifted", type=Path, default=None)
    args = ap.parse_args()

    if args.coords_table:
        row = _lookup_coords(args.coords_table, args.peak_id)
        if row is None:
            sys.exit(f"project_variants: {args.peak_id} not in {args.coords_table}")
        args.cons_chr, args.cons_start, args.cons_end = row["cons_chr"], int(row["cons_start"]), int(row["cons_end"])
        if str(row.get("exists_in_hg38", "")).upper() in ("TRUE", "1"):
            args.hg38_chr, args.hg38_start, args.hg38_end = row["hg38_chr"], row["hg38_start"], row["hg38_end"]
    if args.cons_chr is None or args.cons_start is None or args.cons_end is None:
        sys.exit("project_variants: need --coords-table or --cons-chr/--cons-start/--cons-end")

    aln_haps = read_fasta(args.aligned)
    re_seqs = read_fasta(args.fa)
    asm = _load_asm(args.tsv)
    rows, unlifted = project_peak(
        args.peak_id, aln_haps, re_seqs, asm, args.phylo_mm,
        args.cons_chr, args.cons_start, args.cons_end,
        args.hg38_chr, args.hg38_start, args.hg38_end,
        args.agc, args.flank, args.threads, args.pad, args.min_carriers,
    )
    with args.out.open("w") as fh:
        fh.write(_header())
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in _COLS) + "\n")
    if args.out_unlifted:
        with args.out_unlifted.open("w") as fh:
            fh.write(_unlifted_header())
            for r in unlifted:
                fh.write("\t".join(str(r[c]) for c in _UNLIFTED_COLS) + "\n")
    print(f"[project_variants] {args.peak_id}: {len(rows)} sites "
          f"({sum(r['is_driver']=='TRUE' for r in rows)} driver rows, "
          f"{len(unlifted)} indel-only drivers unlifted)", file=sys.stderr)


if __name__ == "__main__":
    main()
