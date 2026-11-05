#!/usr/bin/env python
"""Mutagenesis panel for multi-change driver blocks — ONE entry per UNIQUE
sequence, with roles + natural identity as ANNOTATIONS (not in the name).

A phylo_mm "driver" can bundle several perfectly-correlated columns (e.g.
RE24/25/26): across the natural haplotypes those changes always co-occur, so the
data can't say which one matters. To break the correlation you test each change
in isolation on a fixed background.

We build, per multi-change driver:
  * background    — a representative carrier with ALL the driver's changes
                    reverted to the ancestral (non-carrier) allele (the negative
                    control / backbone the single edits sit on).
  * <one per change> — background plus exactly ONE change added back.
  * full          — the representative carrier itself (positive control).

But the SAME DNA can be several of these at once and can equal a real haplotype:
`full` is literally the carrier sequence; a `background` is often an existing
non-carrier; one driver's `background` could be another driver's something. So
this DEDUPLICATES by sequence and annotates each unique sequence with:
  - roles=  a |-list of {driver}:{role}[:{edit}] it plays (e.g. d3:background,
            d3:RE24:a>g, d3:full, d6:RE77:g>a)
  - natural= the re_variants seq_id it matches (if any), with that cluster's
            n_haps and MEAN % actuation — i.e. "this construct is really N
            existing haplotypes" — else natural=none (a novel synthetic seq).

Inputs are already-produced files (alignment + phylo_mm TSV + re_variants TSV);
no re-modeling. Run via preextracted_dir (flank=0) peaks or any per-peak run.

Usage:
    make_mutagenesis_panel.py --aligned aln.fa --phylo-mm peak.tsv \
        --re-variants peak.tsv [--driver d3] --out panel.fa \
        [--aligned-out panel.aln.fa] [--tsv panel.tsv]
"""
import argparse
from pathlib import Path


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in Path(path).read_text().splitlines():
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


def parse_phylo_mm(path):
    header, cols, drivers = {}, [], []
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
                    "cols": [int(c) for c in r["linked_sites"].split(",") if c],
                    "carriers": {d for d in r.get("descendants", "").split(",") if d},
                    "beta": r.get("beta", "NA"),
                })
    return header, drivers


def parse_re_variants(path):
    """re_variants TSV → {re_seq: (seq_id, n_haps, pct_acc)} + {sample: seq_id}."""
    by_seq, by_sample = {}, {}
    if not path or not Path(path).exists():
        return by_seq, by_sample
    lines = Path(path).read_text().splitlines()
    cols = lines[0].split("\t")
    for line in lines[1:]:
        r = dict(zip(cols, line.split("\t")))
        rs, sid = r.get("re_seq"), r.get("seq_id")
        if rs and sid:
            by_seq[rs] = (sid, r.get("n_haps", "?"), _pct(r.get("act_mean")))
            for s in (r.get("samples") or "").split(","):
                if s:
                    by_sample[s] = sid
    return by_seq, by_sample


def _num(x, default):
    try:
        return int(x)
    except (TypeError, ValueError):
        return default


def _pct(x):
    try:
        return f"{float(x) * 100:.1f}"
    except (TypeError, ValueError):
        return "NA"


def gc_frac(seq):
    bases = [b for b in seq.upper() if b in "ACGT"]
    return round(sum(b in "GC" for b in bases) / len(bases), 4) if bases else 0.0


def major(chars):
    counts = {}
    for c in chars:
        counts[c] = counts.get(c, 0) + 1
    return max(counts, key=lambda k: (counts[k], -list(counts).index(k)))


def merge_changes(cols, anc, car):
    cols = sorted(cols)
    is_indel = {c: ("-" in (anc[c], car[c])) for c in cols}
    groups, i = [], 0
    while i < len(cols):
        if is_indel[cols[i]]:
            j = i
            while j + 1 < len(cols) and cols[j + 1] == cols[j] + 1 and is_indel[cols[j + 1]]:
                j += 1
            groups.append(cols[i:j + 1])
            i = j + 1
        else:
            groups.append([cols[i]])
            i += 1
    return groups


MANIFEST_COLS = [
    "peak_id", "name", "re_len", "gc_frac", "n_roles", "roles",
    "natural_seq_id", "natural_n_haps", "natural_pct_acc", "is_synthetic", "seq",
]


def build(aligned, phylo_mm, re_variants=None, which=None):
    header, drivers = parse_phylo_mm(phylo_mm)
    natural, _by_sample = parse_re_variants(re_variants)
    peak = header.get("peak_id", Path(phylo_mm).stem)
    re_first = _num(header.get("re_first_col"), 1)
    mode = header.get("attribution_mode", "none")
    seqs = read_fasta(aligned)
    aln_len = len(next(iter(seqs.values()))) if seqs else 0
    re_last = _num(header.get("re_last_col"), aln_len)
    lo, hi = (1, aln_len) if mode == "region" or re_last <= 0 else (re_first, re_last)

    def relabel(c):
        return f"RE{c - re_first + 1}"

    # seq -> {aligned, roles[], anchor}  (deduplicated across drivers + roles)
    uniq, skipped = {}, []
    for drv in drivers:
        if which and drv["id"] != which:
            continue
        carriers = sorted(c for c in drv["carriers"] if c in seqs)   # deterministic rep
        noncarriers = [n for n in seqs if n not in set(drv["carriers"])]
        if not carriers or not noncarriers:
            continue
        cols = [c for c in drv["cols"] if 1 <= c <= aln_len]
        anc = {c: major([seqs[n][c - 1] for n in noncarriers]) for c in cols}
        car = {c: major([seqs[n][c - 1] for n in carriers]) for c in cols}
        changes = merge_changes(cols, anc, car)
        if len(changes) < 2:
            skipped.append((drv["id"], len(changes)))
            continue
        rep = next(c for c in carriers if seqs[c] == major([seqs[x] for x in carriers]))
        full = list(seqs[rep])
        backbone = full.copy()
        for c in cols:
            backbone[c - 1] = anc[c]

        def add(chars, role):
            sub = "".join(chars)[lo - 1:hi].replace("-", "")
            rec = uniq.setdefault(sub, {"aligned": "".join(chars), "roles": [], "anchor": rep})
            if role not in rec["roles"]:
                rec["roles"].append(role)

        add(backbone, f"{drv['id']}:background")
        for grp in changes:
            seq = backbone.copy()
            for c in grp:
                seq[c - 1] = car[c]
            if len(grp) == 1:
                c = grp[0]
                role = f"{drv['id']}:{relabel(c)}:{anc[c]}>{car[c]}"
            else:
                role = f"{drv['id']}:{relabel(grp[0])}-{relabel(grp[-1])}:indel"
            add(seq, role)
        add(full, f"{drv['id']}:full")

    # one record per UNIQUE sequence, roles + natural identity as annotations.
    records = []
    for i, (sub, rec) in enumerate(
        sorted(uniq.items(), key=lambda kv: (-len(kv[1]["roles"]), kv[0]))
    ):
        nat = natural.get(sub)            # (seq_id, n_haps, pct_acc) or None
        records.append({
            "name": f"{peak}_u{i:02d}", "seq": sub, "aligned": rec["aligned"],
            "anchor": rec["anchor"], "roles": rec["roles"], "natural": nat,
            "re_len": len(sub), "gc": gc_frac(sub), "peak": peak,
        })
    return records, skipped


VALIDATION_COLS = [
    "peak_id", "source", "id", "re_len", "gc_frac", "n_haps", "mean_pct_acc",
    "drivers", "roles", "clean_pair", "minimal_edit", "act_values", "samples", "re_seq",
]


def build_validation_table(records, re_variants_path, peak):
    """Unified per-peak "sequences to clone" table: every natural unique RE
    sequence (from re_variants) PLUS the synthetic add-back constructs, deduped
    on sequence. A construct that equals a real haplotype merges into that
    natural row (source=natural) with its role(s) attached; constructs with no
    natural match are their own source=synthetic rows (mean_pct_acc=NA — that's
    what you'd measure). Natural sequences that play no construct role still
    appear (roles empty)."""
    natural = {}   # seq_id -> row
    if re_variants_path and Path(re_variants_path).exists():
        lines = Path(re_variants_path).read_text().splitlines()
        cols = lines[0].split("\t")
        for line in lines[1:]:
            r = dict(zip(cols, line.split("\t")))
            sid = r.get("seq_id")
            if not sid:
                continue
            natural[sid] = {
                "peak_id": peak, "source": "natural", "id": sid,
                "re_len": r.get("re_len", ""), "gc_frac": r.get("gc_frac", ""),
                "n_haps": r.get("n_haps", ""), "mean_pct_acc": _pct(r.get("act_mean")),
                "drivers": r.get("driver_ids", ""), "roles": [],
                "clean_pair": r.get("clean_pair", ""),
                "minimal_edit": r.get("minimal_edit", ""),
                "act_values": r.get("act_values", ""), "samples": r.get("samples", ""),
                "re_seq": r.get("re_seq", ""),
            }
    synthetic = []
    for rec in records:
        if rec["natural"] and rec["natural"][0] in natural:
            natural[rec["natural"][0]]["roles"].extend(rec["roles"])
        elif rec["natural"] is None:
            synthetic.append({
                "peak_id": peak, "source": "synthetic", "id": rec["name"],
                "re_len": rec["re_len"], "gc_frac": rec["gc"], "n_haps": "0",
                "mean_pct_acc": "NA",
                "drivers": ",".join(sorted({role.split(":")[0] for role in rec["roles"]})),
                "roles": list(rec["roles"]), "clean_pair": "", "minimal_edit": "",
                "act_values": "", "samples": "", "re_seq": rec["seq"],
            })
    rows = list(natural.values()) + synthetic
    for r in rows:
        r["roles"] = "|".join(r["roles"])
    return rows


def _desc(r):
    nat = r["natural"]
    nat_s = (f"natural={nat[0]} natural_n_haps={nat[1]} natural_pct_acc={nat[2]}"
             if nat else "natural=none natural_n_haps=0 natural_pct_acc=NA")
    return (f"RE_ID={r['peak']} re_len={r['re_len']} gc={r['gc']} "
            f"{nat_s} n_roles={len(r['roles'])} roles={'|'.join(r['roles'])}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aligned", required=True)
    ap.add_argument("--phylo-mm", required=True)
    ap.add_argument("--re-variants", help="per-peak re_variants TSV for natural-identity annotation")
    ap.add_argument("--driver", help="restrict to one driver id, e.g. d3")
    ap.add_argument("--out", required=True, help="clonable (ungapped) FASTA, unique seqs")
    ap.add_argument("--aligned-out", help="alignment-frame FASTA for the validation plot")
    ap.add_argument("--tsv", help="optional construct manifest TSV (synthetic only)")
    ap.add_argument("--validation-table",
                    help="unified natural + synthetic 'sequences to clone' TSV "
                         "(needs --re-variants)")
    args = ap.parse_args()

    records, skipped = build(args.aligned, args.phylo_mm, args.re_variants, args.driver)

    # Uniqueness invariant: no two records share a sequence.
    assert len({r["seq"] for r in records}) == len(records), "duplicate sequences emitted"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        for r in records:
            fh.write(f">{r['name']} {_desc(r)}\n{r['seq']}\n")
    if args.aligned_out:
        Path(args.aligned_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.aligned_out, "w") as fh:
            for r in records:
                fh.write(f">{r['name']} anchor_leaf={r['anchor']} "
                         f"label={'|'.join(r['roles'])} "
                         f"natural={r['natural'][0] if r['natural'] else 'none'}\n{r['aligned']}\n")
    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("\t".join(MANIFEST_COLS) + "\n")
            for r in records:
                nat = r["natural"]
                row = {
                    "peak_id": r["peak"], "name": r["name"], "re_len": r["re_len"],
                    "gc_frac": r["gc"], "n_roles": len(r["roles"]),
                    "roles": ";".join(r["roles"]),
                    "natural_seq_id": nat[0] if nat else "",
                    "natural_n_haps": nat[1] if nat else "",
                    "natural_pct_acc": nat[2] if nat else "",
                    "is_synthetic": nat is None, "seq": r["seq"],
                }
                fh.write("\t".join(str(row[c]) for c in MANIFEST_COLS) + "\n")

    if args.validation_table:
        peak0 = records[0]["peak"] if records else Path(args.phylo_mm).stem
        vrows = build_validation_table(records, args.re_variants, peak0)
        Path(args.validation_table).parent.mkdir(parents=True, exist_ok=True)
        with open(args.validation_table, "w") as fh:
            fh.write("\t".join(VALIDATION_COLS) + "\n")
            for r in vrows:
                fh.write("\t".join(str(r[c]) for c in VALIDATION_COLS) + "\n")
        n_nat = sum(1 for r in vrows if r["source"] == "natural")
        print(f"[panel] validation table: {len(vrows)} rows "
              f"({n_nat} natural, {len(vrows) - n_nat} synthetic)")

    peak = records[0]["peak"] if records else "(none)"
    n_syn = sum(1 for r in records if r["natural"] is None)
    print(f"[panel] {peak}: {len(records)} unique sequences "
          f"({n_syn} synthetic, {len(records) - n_syn} = existing haplotypes)")
    for did, nch in skipped:
        print(f"[panel] skipped {did}: {nch} change(s), nothing to dissect")


def _self_test():
    import tempfile, os
    d = tempfile.mkdtemp()
    aln = os.path.join(d, "a.fa")
    with open(aln, "w") as f:
        f.write(">car1\nAGTTAC\n>car2\nAGTTGC\n>non1\nAACAAC\n>non2\nAACAAC\n")
    pm = os.path.join(d, "pm.tsv")
    with open(pm, "w") as f:
        f.write("# peak_id=chrT_1_6\n# attribution_mode=RE\n# re_first_col=1\n# re_last_col=6\n")
        f.write("\t".join(["branch_id", "is_split", "linked_sites", "descendants", "beta"]) + "\n")
        f.write("\t".join(["1", "TRUE", "2,3,4", "car1,car2", "-0.4"]) + "\n")
    rv = os.path.join(d, "rv.tsv")
    with open(rv, "w") as f:
        f.write("seq_id\tsamples\tn_haps\tact_mean\tre_seq\n")
        # the 'full' (carrier) AGTTAC and the 'background' AACAAC are both natural
        f.write("chrT_1_6_s00\tcar1\t1\t0.10\tAGTTAC\n")
        f.write("chrT_1_6_s01\tnon1,non2\t2\t0.80\tAACAAC\n")
    recs, skip = build(aln, pm, rv)
    seqs = [r["seq"] for r in recs]
    assert len(seqs) == len(set(seqs)), "not unique"
    byseq = {r["seq"]: r for r in recs}
    # full = carrier AGTTAC, annotated natural=s00 AND role d1:full
    full = byseq["AGTTAC"]
    assert "d1:full" in full["roles"] and full["natural"][0] == "chrT_1_6_s00", full
    # background = AACAAC, natural=s01 (2 haps) AND role d1:background
    bg = byseq["AACAAC"]
    assert "d1:background" in bg["roles"] and bg["natural"][0] == "chrT_1_6_s01", bg
    assert bg["natural"][1] == "2", bg["natural"]
    # the single edits are synthetic (no natural match)
    syn = [r for r in recs if r["natural"] is None]
    assert all(any(":RE" in role for role in r["roles"]) for r in syn), syn
    print("self-test OK:", len(recs), "unique;",
          "synthetic:", len(syn), "| full+bg are existing haplotypes")


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _self_test()
    else:
        main()
