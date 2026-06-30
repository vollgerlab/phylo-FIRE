#!/usr/bin/env python
"""Per-driver actuation violins: the Δ each contributing variant represents.

For every selected driver this draws a split violin at one x position — the
LEFT half is the % actuation of the haplotypes that CARRY the driver, the RIGHT
half is the haplotypes it's contrasted against (the rest). Individual haplotype
points are overlaid, the mean of each side is marked, and the marginal Δ
(mean_in − mean_out) + the joint β are annotated. The separation between the two
halves IS the effect we highlight on the tree/MSA.

Reads the per-peak phylo_mm TSV (drivers: carriers + mean_in/out + β) and the
actuation TSV (per-hap %). No re-modeling.

Usage:
    plot_driver_violins.py --phylo-mm peak.tsv --actuation peak.tsv \
        --peak-id chrX_.. --html out.html --pdf out.pdf
"""
import argparse
from pathlib import Path

import plotly.graph_objects as go


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def parse_phylo_mm(path):
    header, head, drivers = {}, [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            k, _, v = line[1:].partition("=")
            header[k.strip()] = v.strip()
        elif not head:
            head = line.split("\t")
        else:
            r = dict(zip(head, line.split("\t")))
            if r.get("is_split", "").upper() == "TRUE":
                drivers.append({
                    "id": f"d{r['branch_id']}",
                    "cols": [int(c) for c in r["linked_sites"].split(",") if c],
                    "carriers": {d for d in r.get("descendants", "").split(",") if d},
                    "mean_in": _f(r.get("mean_in")), "mean_out": _f(r.get("mean_out")),
                    "beta": _f(r.get("beta")), "p": _f(r.get("p")),
                })
    return header, drivers


def load_act(path):
    out = {}
    lines = Path(path).read_text().splitlines()
    cols = lines[0].split("\t")
    si, ai = cols.index("sample"), cols.index("actuation")
    for line in lines[1:]:
        p = line.split("\t")
        v = _f(p[ai]) if len(p) > ai else None
        if v is not None:
            out[p[si]] = v
    return out


def re_label(cols, re_first):
    labs = [f"RE{c - re_first + 1}" for c in sorted(cols)]
    return ", ".join(labs[:3]) + (f", +{len(labs) - 3}" if len(labs) > 3 else "")


def build(phylo_mm, actuation, peak_id):
    header, drivers = parse_phylo_mm(phylo_mm)
    act = load_act(actuation)
    rf = header.get("re_first_col", "NA")
    re_first = int(rf) if rf not in ("NA", "") else 1
    drivers = [d for d in drivers if d["carriers"]]

    fig = go.Figure()
    cx, cy, ch, rx, ry, rh, anns = [], [], [], [], [], [], []
    for d in drivers:
        lab = f"{d['id']}<br>{re_label(d['cols'], re_first)}"
        carriers = [(s, act[s]) for s in d["carriers"] if s in act]
        rest = [(s, v) for s, v in act.items() if s not in d["carriers"]]
        cx += [lab] * len(carriers); cy += [v * 100 for _, v in carriers]
        ch += [f"{s}: {v:.0%}" for s, v in carriers]
        rx += [lab] * len(rest); ry += [v * 100 for _, v in rest]
        rh += [f"{s}: {v:.0%}" for s, v in rest]
        delta = ((d["mean_in"] - d["mean_out"]) * 100
                 if d["mean_in"] is not None and d["mean_out"] is not None else None)
        bits = []
        if delta is not None:
            bits.append(f"Δ={delta:+.0f}%")
        if d["beta"] is not None:
            bits.append(f"β={d['beta']:+.2f}")
        bits.append(f"n={len(carriers)} vs {len(rest)}")
        anns.append(dict(x=lab, y=108, text="  ".join(bits), showarrow=False,
                         font=dict(size=11), xanchor="center"))

    fig.add_trace(go.Violin(
        x=cx, y=cy, name="carriers (with driver)", side="negative",
        line_color="#2c7fb8", fillcolor="rgba(44,127,184,0.3)",
        points="all", pointpos=-0.4, jitter=0.3, marker=dict(size=4),
        meanline_visible=True, hovertext=ch, hoverinfo="text+y", scalemode="width",
        spanmode="hard",
    ))
    fig.add_trace(go.Violin(
        x=rx, y=ry, name="compared against (rest)", side="positive",
        line_color="#999999", fillcolor="rgba(150,150,150,0.3)",
        points="all", pointpos=0.4, jitter=0.3, marker=dict(size=4),
        meanline_visible=True, hovertext=rh, hoverinfo="text+y", scalemode="width",
        spanmode="hard",
    ))
    fig.update_layout(
        violinmode="overlay", violingap=0,
        title=f"{peak_id} — driver actuation contrasts (carriers vs rest)",
        yaxis_title="% FIRE actuation", template="simple_white",
        annotations=anns, legend=dict(orientation="h", y=1.08),
        width=max(420, 300 * max(1, len(drivers))), height=560,
    )
    fig.update_yaxes(range=[-5, 115])
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phylo-mm", required=True)
    ap.add_argument("--actuation", required=True)
    ap.add_argument("--peak-id", required=True)
    ap.add_argument("--html")
    ap.add_argument("--pdf")
    a = ap.parse_args()
    fig = build(a.phylo_mm, a.actuation, a.peak_id)
    if a.html:
        Path(a.html).parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(a.html, include_plotlyjs="cdn")
        print(f"Saved {a.html}")
    if a.pdf:
        Path(a.pdf).parent.mkdir(parents=True, exist_ok=True)
        fig.write_image(a.pdf, scale=2)
        print(f"Saved {a.pdf}")


def _self_test():
    import tempfile, os
    d = tempfile.mkdtemp()
    pm = os.path.join(d, "pm.tsv")
    with open(pm, "w") as f:
        f.write("# peak_id=chrT_1_9\n# re_first_col=1\n")
        f.write("\t".join(["branch_id", "is_split", "linked_sites", "descendants",
                            "mean_in", "mean_out", "beta", "p"]) + "\n")
        f.write("\t".join(["1", "TRUE", "3,4", "a,b", "0.1", "0.8", "-0.6", "1e-4"]) + "\n")
    ac = os.path.join(d, "ac.tsv")
    with open(ac, "w") as f:
        f.write("sample\tactuation\tcoverage\n")
        for s, v in [("a", 0.1), ("b", 0.12), ("c", 0.8), ("e", 0.75)]:
            f.write(f"{s}\t{v}\t10\n")
    fig = build(pm, ac, "chrT_1_9")
    # two violin traces; carriers trace has a,b → 2 points
    assert len(fig.data) == 2, fig.data
    car = [t for t in fig.data if "carriers" in t.name][0]
    assert len(car.y) == 2, car.y
    rest = [t for t in fig.data if "rest" in t.name][0]
    assert len(rest.y) == 2, rest.y
    print("self-test OK: carriers", list(car.y), "rest", list(rest.y))


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        _self_test()
    else:
        main()
