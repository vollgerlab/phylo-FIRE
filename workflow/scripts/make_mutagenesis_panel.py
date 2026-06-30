#!/usr/bin/env python
"""Generate a single-change mutagenesis panel for multi-change driver blocks.

A phylo_mm "driver" can bundle several perfectly-correlated alignment columns
(e.g. RE24/25/26): across the natural haplotypes those changes ALWAYS co-occur,
so the data can't say which one drives actuation. To localize the effect you
have to break the correlation experimentally — test each change in isolation on
a constant background.

For each selected driver with >1 change this emits:
  * <driver>_background   — a representative carrier with ALL the driver's
                            changes reverted to the ancestral (non-carrier)
                            allele. The negative control; expect ~baseline.
  * <driver>_<edit>       — the background plus exactly ONE change added back.
  * <driver>_full         — the representative carrier itself. Positive control.

Every sequence differs from <driver>_background by exactly one change, on the
carrier's own background, so a reporter assay reads out each change's effect.
Contiguous gap columns (one indel) merge into a single change.

Headers point back to the re_variants table: the source carrier seq_id, its hap
support and observed % actuation, the baseline (nearest non-carrier) seq_id +
its support/actuation, the driver + block beta, and the specific edit — so each
clonable sequence is self-describing and traceable.

Usage:
    make_mutagenesis_panel.py --aligned aln.fa --phylo-mm peak.tsv \
        [--re-variants peak.tsv] [--driver d3] --out panel.fa [--tsv panel.tsv]
"""
import argparse
from pathlib import Path


def read_fasta(path):
    seqs, name, buf = {}, None, []
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


def parse_phylo_mm(path):
    header, cols, drivers = {}, [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                k, _, v = line[1:].strip().partition("=")
                header[k.strip()] = v.strip()
            elif not cols:
                cols = line.split("\t")
            else:
                row = dict(zip(cols, line.split("\t")))
                if row.get("is_split", "").upper() == "TRUE":
                    drivers.append({
                        "id": f"d{row.get('branch_id', '?')}",
                        "cols": [int(c) for c in row["linked_sites"].split(",") if c],
                        "carriers": [d for d in row.get("descendants", "").split(",") if d],
                        "beta": row.get("beta", "NA"),
                    })
    return header, drivers


def parse_re_variants(path):
    """re_variants TSV → ({seq_id: row}, {sample: seq_id}). Empty if no path."""
    by_seq, by_sample = {}, {}
    if not path or not Path(path).exists():
        return by_seq, by_sample
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(cols, line.rstrip("\n").split("\t")))
            sid = r.get("seq_id")
            if not sid:
                continue
            by_seq[sid] = r
            for s in (r.get("samples") or "").split(","):
                if s:
                    by_sample[s] = sid
    return by_seq, by_sample


def _num(x, default):
    try:
        return int(x)
    except (TypeError, ValueError):
        return default


def _pct(act):
    try:
        return f"{float(act) * 100:.1f}"
    except (TypeError, ValueError):
        return "NA"


def major(chars):
    counts = {}
    for c in chars:
        counts[c] = counts.get(c, 0) + 1
    return max(counts, key=lambda k: (counts[k], -list(counts).index(k)))


def merge_changes(cols, anc, car):
    """Group driver columns into changes: contiguous gap-edit columns (one indel)
    merge; substitutions stay individual. Returns list of column-lists."""
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
    "peak_id", "driver", "name", "type", "edit", "n_changes", "beta_block",
    "from_carrier_seq", "carrier_nhaps", "carrier_pct_acc",
    "baseline_seq", "baseline_nhaps", "baseline_pct_acc", "re_len", "seq",
]


def build(aligned, phylo_mm, re_variants=None, which=None):
    header, drivers = parse_phylo_mm(phylo_mm)
    by_seq, by_sample = parse_re_variants(re_variants)
    peak = header.get("peak_id", Path(phylo_mm).stem)
    re_first = _num(header.get("re_first_col"), 1)
    mode = header.get("attribution_mode", "none")
    seqs = read_fasta(aligned)
    aln_len = len(next(iter(seqs.values()))) if seqs else 0
    re_last = _num(header.get("re_last_col"), aln_len)
    lo, hi = (1, aln_len) if mode == "region" or re_last <= 0 else (re_first, re_last)

    def relabel(c):
        return f"RE{c - re_first + 1}"

    records, aligned_records, manifest, skipped = [], [], [], []
    for drv in drivers:
        if which and drv["id"] != which:
            continue
        carriers = [c for c in drv["carriers"] if c in seqs]
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
        # representative carrier background = modal carrier aligned sequence
        modal = major([seqs[c] for c in carriers])
        rep = next(c for c in carriers if seqs[c] == modal)
        full = list(seqs[rep])
        backbone = full.copy()
        for c in cols:
            backbone[c - 1] = anc[c]

        # ── tie back to the re_variants table ────────────────────────────────
        car_seq = by_sample.get(rep, "")
        car_row = by_seq.get(car_seq, {})
        base_seq = car_row.get("nearest_noncarrier_seq_id", "")
        base_row = by_seq.get(base_seq, {})
        ann = {
            "from_carrier_seq": car_seq or "NA",
            "carrier_nhaps": car_row.get("n_haps", "NA"),
            "carrier_pct_acc": _pct(car_row.get("act_mean")),
            "baseline_seq": base_seq or "NA",
            "baseline_nhaps": base_row.get("n_haps", "NA"),
            "baseline_pct_acc": _pct(base_row.get("act_mean")),
        }
        beta = drv["beta"]

        def emit(suffix, chars, edit, kind, npos):
            sub = "".join(chars)[lo - 1:hi].replace("-", "")
            name = f"{peak}_{drv['id']}_{suffix}"
            # alignment-frame copy (gapped, full length) for the validation plot:
            # anchored next to the representative carrier leaf it was built from.
            aligned_records.append({
                "name": name, "seq": "".join(chars), "anchor": rep,
                "label": f"{drv['id']}_{suffix}", "type": kind, "edit": edit,
            })
            desc = (
                f"RE_ID={peak} driver={drv['id']} beta_block={beta} "
                f"type={kind} edit={edit} n_changes={npos} "
                f"from_carrier={ann['from_carrier_seq']} "
                f"carrier_nhaps={ann['carrier_nhaps']} "
                f"carrier_pct_acc={ann['carrier_pct_acc']} "
                f"baseline={ann['baseline_seq']} "
                f"baseline_nhaps={ann['baseline_nhaps']} "
                f"baseline_pct_acc={ann['baseline_pct_acc']} re_len={len(sub)}"
            )
            records.append((name, desc, sub))
            manifest.append({
                "peak_id": peak, "driver": drv["id"], "name": name, "type": kind,
                "edit": edit, "n_changes": npos, "beta_block": beta,
                "from_carrier_seq": ann["from_carrier_seq"],
                "carrier_nhaps": ann["carrier_nhaps"],
                "carrier_pct_acc": ann["carrier_pct_acc"],
                "baseline_seq": ann["baseline_seq"],
                "baseline_nhaps": ann["baseline_nhaps"],
                "baseline_pct_acc": ann["baseline_pct_acc"],
                "re_len": len(sub), "seq": sub,
            })

        emit("background", backbone, "none(all_reverted)", "background", 0)
        for grp in changes:
            seq = backbone.copy()
            for c in grp:
                seq[c - 1] = car[c]
            if len(grp) == 1:
                c = grp[0]
                edit = f"{relabel(c)}:{anc[c]}>{car[c]}"
            else:
                edit = f"{relabel(grp[0])}-{relabel(grp[-1])}:indel"
            emit(edit.split(":")[0], seq, edit, "single_edit", 1)
        emit("full", full, "all(=carrier)", "full", len(changes))
    return records, aligned_records, manifest, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aligned", required=True)
    ap.add_argument("--phylo-mm", required=True)
    ap.add_argument("--re-variants", help="per-peak re_variants TSV for header annotation")
    ap.add_argument("--driver", help="restrict to one driver id, e.g. d3")
    ap.add_argument("--out", required=True, help="clonable (ungapped) FASTA")
    ap.add_argument("--aligned-out", help="alignment-frame FASTA for the validation plot")
    ap.add_argument("--tsv", help="optional manifest TSV")
    args = ap.parse_args()

    records, aligned_records, manifest, skipped = build(
        args.aligned, args.phylo_mm, args.re_variants, args.driver)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        for name, desc, seq in records:
            fh.write(f">{name} {desc}\n{seq}\n")
    if args.aligned_out:
        Path(args.aligned_out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.aligned_out, "w") as fh:
            for r in aligned_records:
                fh.write(
                    f">{r['name']} anchor_leaf={r['anchor']} type={r['type']} "
                    f"edit={r['edit']} label={r['label']}\n{r['seq']}\n"
                )
    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("\t".join(MANIFEST_COLS) + "\n")
            for m in manifest:
                fh.write("\t".join(str(m[c]) for c in MANIFEST_COLS) + "\n")
    peak = manifest[0]["peak_id"] if manifest else "(none)"
    n_drv = len({m["driver"] for m in manifest})
    print(f"[panel] {peak} → {len(records)} seqs over {n_drv} multi-change driver(s)")
    for did, n in skipped:
        print(f"[panel] skipped {did}: {n} change(s), nothing to dissect")


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
        f.write("seq_id\tsamples\tn_haps\tact_mean\tnearest_noncarrier_seq_id\n")
        f.write("chrT_1_6_s00\tcar1,car2\t2\t0.10\tchrT_1_6_s01\n")
        f.write("chrT_1_6_s01\tnon1,non2\t2\t0.80\t\n")
    recs, aligned, _man, _skip = build(aln, pm, rv)
    bg = next(s for n, _, s in recs if n.endswith("_background"))
    # aligned-frame records carry the anchor leaf + full-length (gapped) seq
    assert all(r["anchor"] in ("car1", "car2") for r in aligned), aligned
    assert {len(r["seq"]) for r in aligned} == {6}, aligned
    assert bg[:4] == "AACA", bg
    assert len(recs) == 5, [n for n, _, _ in recs]
    for n, _, s in recs:
        suf = n.split("_")[-1]
        if suf.startswith("RE"):
            assert sum(a != b for a, b in zip(s, bg)) == 1, (n, s, bg)
    # header carries table provenance
    desc0 = next(desc for n, desc, _ in recs if n.endswith("_background"))
    assert "from_carrier=chrT_1_6_s00" in desc0 and "carrier_pct_acc=10.0" in desc0, desc0
    assert "baseline=chrT_1_6_s01" in desc0 and "baseline_pct_acc=80.0" in desc0, desc0
    print("self-test OK:", len(recs), "seqs;", desc0)


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _self_test()
    else:
        main()
