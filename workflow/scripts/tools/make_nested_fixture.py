#!/usr/bin/env python
"""Fabricate a GENUINE nested-driver test peak by injecting a reversing variant
into a SUBSET of a real driver's carriers.

No real peak in the cohort has genuine nesting (an inner variant Z, carried by a
subset of an outer driver X's carriers, that REVERSES X's effect with opposite-
sign β). This builds one from a clean single-driver peak (chr1_157664372: a
strong repressive driver — 35 low-actuation carriers vs 42 high non-carriers):

  group        sequence              actuation   role
  non-X (~42)  baseline              ~0.55       X non-carrier (high)
  X-only       baseline + X SNPs     ~0.05       X represses (low)
  X+Z          baseline + X + Z SNP  ~0.35       Z reverses X (PARTWAY up)

Z's reversal is PARTIAL on purpose: a full reversal would make actuation stop
tracking the sequence tree (the X clade would be half-high/half-low), tank H²ₚ,
and the dual-H² gate would then attribute nothing. Partial keeps the whole X
clade lower than non-X (H² identified) while the inner Z still lifts X+Z above
X-only — so forward selection reports X (β≈−0.50) AND Z (β≈+0.30, opposite
sign), with Z's carriers a strict subset of X's. Run via preextracted_dir
(flank=0).

Usage:
    make_nested_fixture.py --fa <peak>.fa --tsv <peak>.tsv \
        --outdir resources/test/preextracted_nested --peak chrNest_100000_100189
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
            buf = []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fa", required=True)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--peak", default="chrNest_100000_100189")
    ap.add_argument("--x-threshold", type=float, default=0.15,
                    help="actuation below this = an X (outer driver) carrier")
    ap.add_argument("--n-inner", type=int, default=15, help="size of the Z subset")
    args = ap.parse_args()

    seqs = read_fasta(args.fa)
    rows = Path(args.tsv).read_text().splitlines()
    hdr = rows[0].split("\t")
    si, ai, ci = hdr.index("sample"), hdr.index("actuation"), hdr.index("coverage")
    act, cov = {}, {}
    for r in rows[1:]:
        p = r.split("\t")
        try:
            act[p[si]] = float(p[ai]); cov[p[si]] = int(float(p[ci]))
        except (ValueError, IndexError):
            continue

    samples = [s for s in seqs if s in act]
    x_carriers = sorted(s for s in samples if act[s] < args.x_threshold)
    non_x = sorted(s for s in samples if act[s] >= args.x_threshold)
    z_subset = set(x_carriers[:args.n_inner])          # inner Z ⊂ X (deterministic)
    x_only = [s for s in x_carriers if s not in z_subset]

    # Pick an invariant RE column to host the synthetic Z SNP, and a new base.
    L = len(next(iter(seqs.values())))
    inv_col = None
    for c in range(L):
        bases = {seqs[s][c] for s in samples}
        if len(bases) == 1 and next(iter(bases)) in "ACGT":
            inv_col = c
            break
    if inv_col is None:
        raise SystemExit("no invariant ACGT column to host the Z variant")
    base0 = seqs[samples[0]][inv_col]
    base1 = {"A": "T", "T": "A", "C": "G", "G": "C"}[base0]

    # Group target actuations (partial Z reversal). Small deterministic jitter so
    # clusters have within-cluster spread (H² needs replication variance).
    def jit(i):
        return ((i * 37) % 7 - 3) * 0.01

    def assign(s, i):
        if s in z_subset:
            return 0.35 + jit(i)          # X+Z: reversed partway back up
        if s in set(x_only):
            return 0.05 + jit(i)          # X-only: repressed
        return 0.55 + jit(i)              # non-X: high

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    fa_lines, tsv_lines = [], ["sample\tactuation\tcoverage\tfire_coverage"]
    for i, s in enumerate(samples):
        seq = list(seqs[s])
        if s in z_subset:
            seq[inv_col] = base1          # inject the Z SNP
        a = max(0.0, min(1.0, assign(s, i)))
        c = cov.get(s, 20)
        fa_lines += [f">{s}", "".join(seq)]
        tsv_lines.append(f"{s}\t{a:.4f}\t{c}\t{round(a * c)}")
    (out / f"{args.peak}.fa").write_text("\n".join(fa_lines) + "\n")
    (out / f"{args.peak}.tsv").write_text("\n".join(tsv_lines) + "\n")

    re_pos = inv_col + 1
    print(f"[nested] {args.peak}: {len(non_x)} non-X (~0.55), {len(x_only)} X-only "
          f"(~0.05), {len(z_subset)} X+Z (~0.35); Z SNP at RE col {re_pos} "
          f"({base0}->{base1}). X carriers={len(x_carriers)} ⊃ Z carriers={len(z_subset)}")


if __name__ == "__main__":
    main()
