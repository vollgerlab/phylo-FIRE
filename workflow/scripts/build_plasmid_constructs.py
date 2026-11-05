#!/usr/bin/env python3
"""Plasmid fiber-seq construct table for a peak (per the experimental-design notes).

For each SINGLE-driver regulatory element (peaks with >1 distinct driver effect
are excluded), pick two representative haplotypes and engineer the driver on/off
both backgrounds:

  * driver rep   = modal RE sequence AMONG CARRIERS of the driver; tie -> highest
                   actuation.
  * outgroup rep = modal RE sequence among NON-carriers;           tie -> lowest
                   actuation.

Then up to 4 constructs (a driver may span several linked columns; ALL of them
are flipped together):
  driver_bg_with     driver rep, native                (+driver)
  driver_bg_without  driver rep, driver cols -> outgroup allele   (-driver)
  outgroup_bg_without outgroup rep, native             (-driver)
  outgroup_bg_with   outgroup rep, driver cols -> driver allele    (+driver)

If the driver is the ONLY difference between the reps the synthetics recreate the
natives -> only 2 unique sequences (flagged `collapsed`). Each rep's modal count
(`source_n_haps`) and mean actuation (`source_act_mean`) are recorded so the
plasmid readout can be checked against the native effect.

``build_peak`` returns the construct rows (build_plasmid_table.py writes them for
the whole run) and the aligned panel that plot_all draws as grey pseudo-leaf rows
under the MSA.
"""
from pathlib import Path


def read_fasta(path):
    out, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        out[name] = "".join(buf)
    return out


def parse_phylo_mm(path):
    """(header dict, [drivers]) — driver = {id, cols(0-based aln), carriers, beta}."""
    header, cols, drivers = {}, [], []
    if not path or not Path(path).exists():
        return header, drivers
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            k, _, v = line[1:].partition("=")
            header[k.strip()] = v.strip()
        elif not cols:
            cols = line.split("\t")
        else:
            r = dict(zip(cols, line.split("\t")))
            if r.get("is_split", "").upper() == "TRUE":
                drivers.append({
                    "id": f"d{r['branch_id']}",
                    "cols": [int(c) - 1 for c in r.get("linked_sites", "").split(",") if c],
                    "carriers": {d for d in r.get("descendants", "").split(",") if d},
                    "beta": r.get("beta", "NA"),
                })
    return header, drivers


def cov_mean(haps, act):
    """coverage-weighted mean actuation over a list of hap names."""
    num = den = 0.0
    for h in haps:
        if h in act:
            a, c = act[h]
            num += a * c
            den += c
    return (num / den) if den else 0.0


def pick_rep(names, aln, act, re_lo, re_hi, prefer_high):
    """Modal RE-band sequence among `names`; tie broken by cov-weighted mean
    actuation (highest if prefer_high else lowest). Returns
    (rep_leaf, n_haps, act_mean) — rep_leaf is a representative haplotype."""
    groups = {}
    for n in names:
        key = aln[n][re_lo:re_hi + 1]
        groups.setdefault(key, []).append(n)
    best = None
    for key, haps in groups.items():
        m = cov_mean(haps, act)
        score = (len(haps), m if prefer_high else -m)
        if best is None or score > best[0]:
            # rep leaf = the group's extreme-actuation member (deterministic)
            rep = max(haps, key=lambda h: act.get(h, (0, 0))[0]) if prefer_high \
                else min(haps, key=lambda h: act.get(h, (1, 0))[0])
            best = (score, rep, len(haps), m)
    return (best[1], best[2], best[3]) if best else (None, 0, 0.0)


COLS = ["peak_id", "driver_id", "driver_beta", "driver_re_cols", "role", "background",
        "carries_driver", "is_synthetic", "collapsed", "source_rep", "source_leaf",
        "source_n_haps", "source_act_mean", "edit", "re_seq", "re_len"]


def build_peak(aligned, act, phylo_mm, peak_id):
    """Build the plasmid constructs for one peak.

    ``aligned`` and ``phylo_mm`` are file paths; ``act`` is {sample: (actuation,
    coverage)}. Returns (rows, panel). rows == [] when the peak is excluded (no
    driver, >1 distinct driver effect, or a rep group is empty).
    """
    aln = read_fasta(aligned)
    header, drivers = parse_phylo_mm(phylo_mm)
    aln_len = len(next(iter(aln.values()))) if aln else 0
    re_lo = int(header.get("re_first_col", 1)) - 1
    re_hi = int(header.get("re_last_col", aln_len)) - 1

    rows, panel = [], []
    # Exclusion: keep only peaks with a SINGLE distinct driver EFFECT (one beta).
    betas = {d["beta"] for d in drivers}
    if len(drivers) == 0 or len(betas) > 1:
        return rows, panel
    drv = drivers[0]
    all_names = [n for n in aln]
    carriers = [n for n in all_names if n in drv["carriers"]]
    outgroup = [n for n in all_names if n not in drv["carriers"]]
    if not carriers or not outgroup:
        return rows, panel

    d_leaf, d_n, d_mean = pick_rep(carriers, aln, act, re_lo, re_hi, prefer_high=True)
    o_leaf, o_n, o_mean = pick_rep(outgroup, aln, act, re_lo, re_hi, prefer_high=False)
    d_seq, o_seq = list(aln[d_leaf]), list(aln[o_leaf])
    cols = [c for c in drv["cols"] if 0 <= c < aln_len]
    if not cols:
        # Flank-driven split with no RE-scoped variant: the driver isn't in a
        # column we can locate, so there's no clonable pair to build -- skip.
        return rows, panel
    # The clonable band spans the RE AND every driver column, so the sequence we
    # emit always CONTAINS the edit. A region-mode (or flank-region) driver sits
    # past re_last_col; slicing the RE band alone drops it, making all four
    # constructs identical over the clonable unit (spurious duplicates).
    band_lo = min(re_lo, min(cols))
    band_hi = max(re_hi, max(cols))
    d_alleles = {c: d_seq[c] for c in cols}   # driver allele
    o_alleles = {c: o_seq[c] for c in cols}   # background/outgroup allele

    def edited(base_seq, alleles):
        s = list(base_seq)
        for c, a in alleles.items():
            s[c] = a
        return s

    def band_seq(seq):
        return "".join(seq[band_lo:band_hi + 1]).replace("-", "")

    def edit_str(frm, to):
        parts = [f"RE{c - band_lo + 1}:{frm[c]}>{to[c]}" for c in cols]
        return ";".join(parts)

    # (role, background, carries, synthetic, aligned_seq, rep, leaf, n, mean, edit)
    specs = [
        ("driver_bg_with", "driver", "TRUE", "FALSE", d_seq, "driver", d_leaf, d_n, d_mean, ""),
        ("driver_bg_without", "driver", "FALSE", "TRUE", edited(d_seq, o_alleles), "driver",
         d_leaf, d_n, d_mean, edit_str(d_alleles, o_alleles)),
        ("outgroup_bg_without", "outgroup", "FALSE", "FALSE", o_seq, "outgroup", o_leaf, o_n, o_mean, ""),
        ("outgroup_bg_with", "outgroup", "TRUE", "TRUE", edited(o_seq, d_alleles), "outgroup",
         o_leaf, o_n, o_mean, edit_str(o_alleles, d_alleles)),
    ]
    # Dedup on the CLONABLE band sequence (== what we actually synthesize),
    # reserving the two natives first so a synthetic that merely recreates a
    # native is the NA'd row. Two constructs identical over the band are the SAME
    # plasmid even when they differ elsewhere in the (uncloned) flank -- so the
    # key is the band, not the full aligned sequence.
    seen = {band_seq(sp[4]) for sp in specs if sp[3] == "FALSE"}
    for role, bg, carries, synth, seq, rep, leaf, n, mean, edit in specs:
        bseq = band_seq(seq)
        collapsed = "TRUE" if (synth == "TRUE" and bseq in seen) else "FALSE"
        seen.add(bseq)
        if collapsed == "TRUE":
            re_seq, re_len = "NA", "NA"
        else:
            re_seq, re_len = bseq, len(bseq)
            panel.append((f"PL_{peak_id}_{drv['id']}_{role}", leaf, role, "".join(seq)))
        rows.append({
            "peak_id": peak_id, "driver_id": drv["id"], "driver_beta": drv["beta"],
            "driver_re_cols": ",".join(str(c - band_lo + 1) for c in cols),
            "role": role, "background": bg, "carries_driver": carries, "is_synthetic": synth,
            "collapsed": collapsed, "source_rep": rep, "source_leaf": leaf,
            "source_n_haps": n, "source_act_mean": f"{mean:.4f}", "edit": edit or "NA",
            "re_seq": re_seq, "re_len": re_len,
        })
    # A testable pair needs >=2 distinct clonable sequences; otherwise the driver
    # doesn't manifest in the clone (e.g. only gap columns) -- drop the peak.
    if len({r["re_seq"] for r in rows if r["re_seq"] != "NA"}) < 2:
        return [], []
    return rows, panel
