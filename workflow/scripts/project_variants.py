#!/usr/bin/env python3
"""Project every variable site in a peak's RE alignment to CHM13 + hg38 (+ each
carrier's assembly), flagging which sites are attribution drivers.

Self-contained: consumes only PERSISTED pipeline outputs (the RE-only FASTA and
its per-hap TSV, and the phylo_mm TSV) plus the AGC, so it can be added to a
finished run without re-running trees/phylo_mm. It builds its OWN small
alignment of [haplotypes + CHM13 ref + hg38 ref] rather than reusing the
pipeline's temp alignment.

Coordinate logic: the CHM13/hg38 reference sequences are pulled FORWARD from the
AGC. mafft --adjustdirection may flip the whole alignment; we detect each ref's
orientation by comparing its ungapped aligned sequence to the forward pull, then
map alignment columns to genomic coordinates (forward or reverse) so an internal
CHM13↔hg38 indel or a minus-strand peak is handled correctly.

Emits one TSV row per (variable site), with CHM13 + hg38 coords, alleles,
carrier count, is_driver + beta, and a comma-sep of per-carrier assembly
coordinates.
"""
import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

_COMP = str.maketrans("ACGTNacgtn-", "TGCANtgcan-")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


def read_fasta(path: Path) -> dict[str, str]:
    out, name, parts = {}, None, []
    for line in Path(path).read_text().splitlines():
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
    """Return the single concatenated sequence from `agc getctg` stdout."""
    return "".join(l.strip() for l in text.splitlines() if l and not l.startswith(">"))


def pull_agc(agc: str, contig: str, start0: int, end0: int) -> str | None:
    """agc getctg <agc> <contig>:<start>-<end>. NOTE: agc coords are 0-BASED
    (verified on the archive: `agc getctg …chrN:N-N` returns the standard 1-based
    base N+1), so callers pass a 0-based start DIRECTLY -- not BED-start+1."""
    region = f"{contig}:{start0}-{end0}"
    try:
        res = subprocess.run(
            ["agc", "getctg", agc, region, "-p"],
            capture_output=True, text=True, check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as ex:
        print(f"[project_variants] agc getctg failed for {region}: {ex}", file=sys.stderr)
        return None
    seq = parse_agc_fasta(res.stdout).upper()
    return seq or None


def mafft_align(seqs: dict[str, str], threads: int) -> dict[str, str]:
    """mafft --adjustdirection over an ordered dict; returns aligned {name: seq}.
    Strips mafft's `_R_` reversed-prefix from names."""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        for nm, sq in seqs.items():
            fh.write(f">{nm}\n{sq}\n")
        inp = fh.name
    try:
        res = subprocess.run(
            ["mafft", "--auto", "--adjustdirection", "--thread", str(threads), inp],
            capture_output=True, text=True, check=True,
        )
    finally:
        Path(inp).unlink(missing_ok=True)
    aln = read_fasta_str(res.stdout)
    return {nm[3:] if nm.startswith("_R_") else nm: sq for nm, sq in aln.items()}


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


def orient_and_map(aln_ref: str, fwd_seq: str, start0: int):
    """Given a ref's aligned string and the FORWARD pulled sequence starting at
    0-based genomic `start0`, return (orientation, col->pos0 map). pos0 is the
    0-based genomic position of the ref base at a column (absent where gapped).

    orientation 'fwd': aligned ungapped == fwd_seq. 'rev': == revcomp(fwd_seq)."""
    ungapped = aln_ref.replace("-", "").upper()
    end0 = start0 + len(fwd_seq)
    if ungapped == fwd_seq.upper():
        orient = "fwd"
    elif ungapped == revcomp(fwd_seq).upper():
        orient = "rev"
    else:
        # Imperfect (SNVs between ref pull and alignment consensus); fall back to
        # whichever direction shares more bases, so a few mismatches don't break it.
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fa", required=True, type=Path, help="fastas/{peak}.fa (RE-only)")
    ap.add_argument("--tsv", required=True, type=Path, help="fastas/{peak}.tsv (asm coords)")
    ap.add_argument("--phylo-mm", type=Path, default=None, help="phylo_mm/{peak}.tsv")
    ap.add_argument("--peak-id", required=True)
    ap.add_argument("--coords-table", type=Path, default=None,
                    help="peak_coords.tsv from build_peak_coords.py; if given, "
                         "cons/hg38 coords are looked up by peak_id from it")
    ap.add_argument("--cons-chr")
    ap.add_argument("--cons-start", type=int, help="0-based (BED)")
    ap.add_argument("--cons-end", type=int)
    ap.add_argument("--hg38-chr", default="NA")
    ap.add_argument("--hg38-start", default="NA")
    ap.add_argument("--hg38-end", default="NA")
    ap.add_argument("--agc", required=True)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--min-carriers", type=int, default=2)
    ap.add_argument("--pad", type=int, default=25,
                    help="bp to pad each ref pull so an RE that overhangs the "
                         "exact peak coords still projects fully")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    # Resolve coords from the table (by peak_id) when supplied.
    if args.coords_table:
        row = _lookup_coords(args.coords_table, args.peak_id)
        if row is None:
            sys.exit(f"project_variants: {args.peak_id} not in {args.coords_table}")
        args.cons_chr, args.cons_start, args.cons_end = row["cons_chr"], int(row["cons_start"]), int(row["cons_end"])
        if str(row.get("exists_in_hg38", "")).upper() in ("TRUE", "1"):
            args.hg38_chr, args.hg38_start, args.hg38_end = row["hg38_chr"], row["hg38_start"], row["hg38_end"]
    if args.cons_chr is None or args.cons_start is None or args.cons_end is None:
        sys.exit("project_variants: need --coords-table or --cons-chr/--cons-start/--cons-end")

    haps = read_fasta(args.fa)
    if len(haps) < 2:
        args.out.write_text(_header() + "")
        print(f"[project_variants] {args.peak_id}: <2 haps, no variants", file=sys.stderr)
        return

    # Per-hap assembly coords (sample -> asm_chr/start/end) from the peak TSV.
    asm = _load_asm(args.tsv)
    # Driver carrier sets from phylo_mm (frozenset(sample) -> beta).
    drivers = _load_drivers(args.phylo_mm) if args.phylo_mm else {}

    # ── pull refs + align [CHM13, hg38, haps] (refs first) ────────────────────
    pad = args.pad
    refs: dict[str, str] = {}
    cons_s0 = max(0, args.cons_start - pad)          # 0-based padded ref start
    # agc getctg is 0-based, so pass cons_s0 directly. The historical `+ 1` here
    # assumed agc was 1-based and shifted every pulled base one position right,
    # making reported CHM13/hg38 coords 1 too low (REF landed at pos+1).
    cons_fwd = pull_agc(args.agc, f"CHM13#0#{args.cons_chr}", cons_s0, args.cons_end + pad)
    if cons_fwd:
        refs["__CHM13__"] = cons_fwd
    hg38_ok = args.hg38_chr != "NA" and str(args.hg38_start).isdigit()
    hg38_s0 = hg38_fwd = None
    if hg38_ok:
        hg38_s0 = max(0, int(args.hg38_start) - pad)
        hg38_fwd = pull_agc(args.agc, f"GRCh38#0#{args.hg38_chr}", hg38_s0, int(args.hg38_end) + pad)
        if hg38_fwd:
            refs["__hg38__"] = hg38_fwd

    combined = {**refs, **haps}
    aln = mafft_align(combined, args.threads)
    aln_len = len(next(iter(aln.values())))

    cons_map = None
    if "__CHM13__" in aln and cons_fwd:
        _, cons_map = orient_and_map(aln["__CHM13__"], cons_fwd, cons_s0)
    hg38_map = None
    if "__hg38__" in aln and hg38_fwd and hg38_s0 is not None:
        _, hg38_map = orient_and_map(aln["__hg38__"], hg38_fwd, hg38_s0)

    hap_names = [n for n in aln if not n.startswith("__")]
    # Alignment orientation vs the genome. Refs are pulled FORWARD; if the CHM13
    # ref aligned as reverse-complement, the whole column set is minus-strand, so
    # haplotype alleles get revcomp'd to report forward-genomic (matching REF).
    aln_is_rev = _rev_orient(aln.get("__CHM13__", ""), cons_fwd) if cons_fwd else False

    rows = []
    for col in range(aln_len):
        col_bases = {n: aln[n][col].upper() for n in hap_names}
        # SNV sites only for now (>=2 distinct NON-gap bases). Indels (gap vs
        # base) are handled separately — TODO.
        counts: dict[str, int] = {}
        for b in col_bases.values():
            if b not in ("-", "N"):
                counts[b] = counts.get(b, 0) + 1
        if len(counts) < 2:
            continue
        minor = min(counts, key=counts.get)
        major = max(counts, key=counts.get)
        carriers = sorted(n for n, b in col_bases.items() if b == minor)
        if min(len(carriers), sum(counts.values()) - len(carriers)) < args.min_carriers:
            continue  # singleton/near-private — mirrors phylo_mm MIN_CARRIERS
        cset = frozenset(carriers)
        is_driver = cset in drivers
        beta = drivers.get(cset, "NA")

        cons_pos0 = cons_map.get(col) if cons_map else None
        hg38_pos0 = hg38_map.get(col) if hg38_map else None
        # REF = the reference base straight from the forward pull at the mapped
        # position (== the true CHM13/hg38 reference base; also self-validates the
        # coordinate — re-fetching that position returns this base).
        cons_ref = cons_fwd[cons_pos0 - cons_s0].upper() if cons_pos0 is not None else "NA"
        hg38_ref = (hg38_fwd[hg38_pos0 - hg38_s0].upper()
                    if (hg38_fwd and hg38_pos0 is not None) else "NA")
        maj_f, min_f = (revcomp(major), revcomp(minor)) if aln_is_rev else (major, minor)

        # per-carrier assembly coordinates (asm_start + ungapped offset in that hap)
        asm_hits = []
        for nm in carriers:
            off = _ungapped_offset(aln[nm], col)
            a = asm.get(nm)
            if a and off is not None and a["start"].isdigit():
                asm_hits.append(f"{nm}={a['chr']}:{int(a['start'])+off}")

        rows.append({
            "peak_id": args.peak_id,
            # position within the peak (1-based, from CHM13 genomic); NA/≤0 for
            # sites in the extraction overhang outside the peak interval.
            "re_pos": (cons_pos0 - args.cons_start + 1) if cons_pos0 is not None else "NA",
            "cons_chr": args.cons_chr,
            "cons_pos": (cons_pos0 + 1) if cons_pos0 is not None else "NA",   # 1-based
            "cons_ref": cons_ref,
            "hg38_chr": args.hg38_chr if hg38_pos0 is not None else "NA",
            "hg38_pos": (hg38_pos0 + 1) if hg38_pos0 is not None else "NA",   # 1-based
            "hg38_ref": hg38_ref,
            "hap_major": maj_f,
            "hap_minor": min_f,
            "n_carriers": len(carriers),
            "is_driver": "TRUE" if is_driver else "FALSE",
            "beta": beta,
            "carriers": ",".join(carriers),
            "asm_coords": ";".join(asm_hits) if asm_hits else "NA",
        })

    with args.out.open("w") as fh:
        fh.write(_header())
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in _COLS) + "\n")
    print(f"[project_variants] {args.peak_id}: {len(rows)} variable sites "
          f"({sum(r['is_driver']=='TRUE' for r in rows)} drivers)", file=sys.stderr)


_COLS = ["peak_id", "re_pos", "cons_chr", "cons_pos", "cons_ref",
         "hg38_chr", "hg38_pos", "hg38_ref", "hap_major", "hap_minor",
         "n_carriers", "is_driver", "beta", "carriers", "asm_coords"]


def _header() -> str:
    return "\t".join(_COLS) + "\n"


def _rev_orient(aln_ref: str, fwd_seq: str | None) -> bool:
    if not aln_ref or not fwd_seq:
        return False
    ung = aln_ref.replace("-", "").upper()
    return ung != fwd_seq.upper() and ung == revcomp(fwd_seq).upper()


def _ungapped_offset(aln_row: str, col: int) -> int | None:
    """0-based count of non-gap chars in aln_row[:col] (the hap's own coord
    offset at that column), or None if the hap is gapped at col."""
    if aln_row[col] == "-":
        return None
    return sum(1 for ch in aln_row[:col] if ch != "-")


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
    hdr = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= max(idx.get("asm_end", 0), idx.get("sample", 0)):
            continue
        out[f[idx["sample"]]] = {
            "chr": f[idx["asm_chr"]], "start": f[idx["asm_start"]], "end": f[idx["asm_end"]],
        }
    return out


def _load_drivers(phylo_mm: Path) -> dict[frozenset, str]:
    """frozenset(carrier leaves) -> beta, for is_split rows in the phylo_mm TSV."""
    out = {}
    if not phylo_mm or not Path(phylo_mm).exists():
        return out
    lines = [l for l in Path(phylo_mm).read_text().splitlines() if not l.startswith("#")]
    if not lines:
        return out
    hdr = lines[0].split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    if "is_split" not in idx or "descendants" not in idx:
        return out
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) <= idx["descendants"] or f[idx["is_split"]].upper() not in ("TRUE", "1"):
            continue
        desc = [d for d in f[idx["descendants"]].split(",") if d]
        beta = f[idx["beta"]] if "beta" in idx and len(f) > idx["beta"] else "NA"
        out[frozenset(desc)] = beta
    return out


if __name__ == "__main__":
    main()
