#!/usr/bin/env python
"""Build the per-peak ``re_variants`` table: one row per unique RE sequence.

The clonable / validatable unit for plasmid work is a regulatory-element
sequence plus (a) its actuation, (b) which driver variant(s) distinguish it,
and (c) how to recreate it from a near-identical non-carrier. This script
assembles that from three already-produced per-peak files — it does NOT redo
any modeling:

  * the ``phylo_mm`` TSV   — drivers (rows with ``is_split=TRUE``: carrier set,
                             linked alignment columns, joint beta, p), plus the
                             header's RE band, attribution mode and H2.
  * the mafft alignment    — the actual sequences (RE band sliced out).
  * the actuation TSV       — per-haplotype % FIRE + coverage.

Grain = one row per unique RE sequence (haplotypes identical over the attributed
band collapse). ``driver_delta`` is the additive model's predicted actuation
shift for that sequence (sum of the joint betas of the driver alleles it
carries). Plasmid helpers (``nearest_noncarrier_seq_id`` / ``n_other_diffs`` /
``clean_pair`` / ``minimal_edit``) name the cheapest single-variable validation
contrast for each driver-carrying sequence.

Library only: re_variants_all.py calls ``build`` for every peak.
"""
import statistics
from pathlib import Path

# A driver set whose linked sites cover more than this fraction of the attributed
# band is flagged `divergent`: one clade differing at that many positions is a
# divergent haplotype (paralog / mismapped copy / misalignment), not a clonable
# variant. Flagged, never filtered: the drivers stay in every output.
DIVERGENT_FRAC = 0.05


def is_divergent(n_sites: int, band_len: int) -> bool:
    """True when a driver set's linked sites exceed DIVERGENT_FRAC of the band."""
    return band_len > 0 and n_sites / band_len > DIVERGENT_FRAC


# Output column order (denormalized: one table, peak context repeated per row).
COLUMNS = [
    "peak_id", "chrom", "start", "end",
    "seq_id", "re_seq", "re_len", "gc_frac",
    "n_haps", "act_mean", "act_sd", "act_min", "act_max", "total_cov",
    "act_values", "samples",
    "has_driver", "driver_ids", "driver_positions", "driver_betas",
    "driver_p", "driver_divergent", "driver_delta",
    "nearest_noncarrier_seq_id", "n_other_diffs", "clean_pair", "minimal_edit",
    "attribution_mode",
    "h2_re_mean", "h2_re_lower", "h2_re_upper", "h2_re_ok",
    "h2_region_mean", "h2_region_lower", "h2_region_upper", "h2_region_ok",
    "n_unique_seqs", "n_variants_selected",
]


def read_fasta(path: str) -> dict[str, str]:
    """Aligned FASTA → {leaf: seq}, stripping mafft's _R_ reverse-prefix."""
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                if name.startswith("_R_"):
                    name = name[3:]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def parse_phylo_mm(path: str) -> tuple[dict[str, str], list[dict]]:
    """Return (header dict, driver rows). Driver rows are the is_split=TRUE
    blocks with their carrier set, linked columns, beta and p."""
    header: dict[str, str] = {}
    col_names: list[str] = []
    drivers: list[dict] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                k, _, v = line[1:].strip().partition("=")
                header[k.strip()] = v.strip()
                continue
            parts = line.split("\t")
            if not col_names:
                col_names = parts
                continue
            row = dict(zip(col_names, parts))
            if row.get("is_split", "").upper() != "TRUE":
                continue
            drivers.append({
                "id": f"d{row.get('branch_id', '?')}",
                "cols": [int(c) for c in row["linked_sites"].split(",") if c],
                "carriers": set(
                    d for d in row.get("descendants", "").split(",") if d
                ),
                "beta": _num(row.get("beta")),
                "p": _num(row.get("p")),
            })
    return header, drivers


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _hnum(header: dict, key: str):
    return _num(header.get(key))


def gc_frac(seq: str) -> float:
    bases = [b for b in seq.upper() if b in "ACGT"]
    if not bases:
        return 0.0
    return round(sum(b in "GC" for b in bases) / len(bases), 4)


def build(aligned: str, actuation: str, phylo_mm: str) -> list[dict]:
    header, drivers = parse_phylo_mm(phylo_mm)
    peak_id = header.get("peak_id", Path(phylo_mm).stem)
    # peak_id = "{chrom}_{start}_{end}" (chrom may itself contain underscores,
    # e.g. chr8_MATERNAL_...), so split off the last two fields.
    chrom, _, rest = peak_id.rpartition("_")
    chrom, _, start = chrom.rpartition("_")
    end = rest
    # phylo_mm reports the chrom/start/end via peak_id; keep raw on parse failure.

    mode = header.get("attribution_mode", "none")
    re_first = header.get("re_first_col", "NA")
    re_first_col = int(re_first) if re_first not in ("NA", "") else 1
    seqs = read_fasta(aligned)
    aln_len = len(next(iter(seqs.values()))) if seqs else 0
    re_last = header.get("re_last_col", "NA")
    re_last_col = int(re_last) if re_last not in ("NA", "") else aln_len

    # Attributed band = the span the drivers live in and the clonable unit.
    # region mode attributes the whole alignment (flank drivers); otherwise the
    # RE band. 1-based inclusive [lo, hi] → 0-based slice [lo-1 : hi].
    if mode == "region":
        lo, hi = 1, aln_len
    else:
        lo, hi = re_first_col, re_last_col
    if hi < lo or hi == 0:
        lo, hi = 1, aln_len
    for d in drivers:
        d["divergent"] = is_divergent(len(d["cols"]), hi - lo + 1)

    # actuation: sample → (act, cov). Accept either a pre-loaded dict (Stage-2
    # single-pass, from haplotypes.parquet) or a per-peak TSV path (CLI).
    if isinstance(actuation, dict):
        act = {s: v for s, v in actuation.items() if v[1] and v[1] > 0}
    else:
        act = {}
        with open(actuation) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            si, ai, ci = cols.index("sample"), cols.index("actuation"), cols.index("coverage")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                a, c = _num(p[ai]), _num(p[ci])
                if a is not None and c is not None and c > 0:
                    act[p[si]] = (a, c)

    # Group tips (with actuation) by their attributed-band aligned substring.
    groups: dict[str, list[str]] = {}
    band: dict[str, str] = {}
    for tip, s in seqs.items():
        if tip not in act:
            continue
        sub = s[lo - 1:hi]
        groups.setdefault(sub, []).append(tip)
        band[tip] = sub

    rows = []
    for sub, members in groups.items():
        acts = [act[m][0] for m in members]
        covs = [act[m][1] for m in members]
        tot_cov = sum(covs)
        wmean = sum(a * c for a, c in zip(acts, covs)) / tot_cov if tot_cov else 0.0
        carried = [d for d in drivers if members[0] in d["carriers"]]
        rows.append({
            "members": members,
            "sub": sub,
            "re_seq": sub.replace("-", ""),
            "acts": acts, "covs": covs, "tot_cov": tot_cov, "wmean": wmean,
            "carried": carried,
        })
    # seq_id by descending actuation (high actuation first), stable.
    rows.sort(key=lambda r: (-r["wmean"], r["sub"]))
    for i, r in enumerate(rows):
        r["seq_id"] = f"{peak_id}_s{i:02d}"

    out = []
    for r in rows:
        carried = r["carried"]
        driver_cols = sorted({c for d in carried for c in d["cols"]})
        nn_id, n_other, clean, edit = _nearest_noncarrier(
            r, rows, lo, driver_cols, re_first_col)
        out.append({
            "peak_id": peak_id, "chrom": chrom, "start": start, "end": end,
            "seq_id": r["seq_id"], "re_seq": r["re_seq"],
            "re_len": len(r["re_seq"]), "gc_frac": gc_frac(r["re_seq"]),
            "n_haps": len(r["members"]),
            "act_mean": round(r["wmean"], 4),
            "act_sd": round(statistics.pstdev(r["acts"]), 4) if len(r["acts"]) > 1 else 0.0,
            "act_min": round(min(r["acts"]), 4),
            "act_max": round(max(r["acts"]), 4),
            "total_cov": r["tot_cov"],
            "act_values": ",".join(f"{a:.4g}" for a in r["acts"]),
            "samples": ",".join(r["members"]),
            "has_driver": bool(carried),
            "driver_ids": ",".join(d["id"] for d in carried),
            "driver_positions": ",".join(
                f"RE{min(d['cols']) - re_first_col + 1}" for d in carried),
            "driver_betas": ",".join(
                f"{d['beta']:.4g}" if d["beta"] is not None else "NA" for d in carried),
            "driver_p": ",".join(
                f"{d['p']:.3g}" if d["p"] is not None else "NA" for d in carried),
            "driver_divergent": ",".join(
                "TRUE" if d["divergent"] else "FALSE" for d in carried),
            "driver_delta": round(sum(d["beta"] or 0.0 for d in carried), 4) if carried else 0.0,
            "nearest_noncarrier_seq_id": nn_id,
            "n_other_diffs": n_other,
            "clean_pair": clean,
            "minimal_edit": edit,
            "attribution_mode": mode,
            "h2_re_mean": _hnum(header, "h2_re_mean"),
            "h2_re_lower": _hnum(header, "h2_re_lower"),
            "h2_re_upper": _hnum(header, "h2_re_upper"),
            "h2_re_ok": header.get("h2_re_ok", "").upper() == "TRUE",
            "h2_region_mean": _hnum(header, "h2_region_mean"),
            "h2_region_lower": _hnum(header, "h2_region_lower"),
            "h2_region_upper": _hnum(header, "h2_region_upper"),
            "h2_region_ok": header.get("h2_region_ok", "").upper() == "TRUE",
            "n_unique_seqs": _intornone(header.get("n_unique_seqs")),
            "n_variants_selected": _intornone(header.get("n_variants_selected")),
        })
    return out


def _intornone(x):
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _nearest_noncarrier(r, rows, lo, driver_cols, re_first_col):
    """For a driver-carrying sequence, find the nearest non-carrier and the
    minimal edit (at driver columns) that converts it into this sequence.

    Returns (nearest_seq_id, n_other_diffs, clean_pair, minimal_edit). n_other_diffs
    counts differences OUTSIDE the driver columns — 0 means a natural clean
    single-variable pair exists (clone both as-is, no mutagenesis)."""
    if not r["carried"]:
        return "", None, False, ""
    dcols0 = {c - lo for c in driver_cols}      # 0-based within the band slice
    best = None
    for o in rows:
        if o is r or o["carried"]:
            continue
        a, b = r["sub"], o["sub"]
        other = sum(
            1 for k in range(min(len(a), len(b)))
            if a[k] != b[k] and k not in dcols0
        )
        if best is None or other < best[1]:
            best = (o, other)
    if best is None:
        return "", None, False, ""
    o, other = best
    edits = []
    for c in driver_cols:
        k = c - lo
        if k < 0 or k >= len(o["sub"]) or k >= len(r["sub"]):
            continue
        frm, to = o["sub"][k], r["sub"][k]
        if frm != to:
            edits.append(f"RE{c - re_first_col + 1} {frm}->{to}")
    return o["seq_id"], other, other == 0, ";".join(edits)
