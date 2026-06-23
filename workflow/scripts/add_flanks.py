#!/usr/bin/env python3
"""Add assembly flanks to a peak's RE-only FASTA by pulling from an AGC archive.

Reads a per-peak RE FASTA (``{peak}.fa``, headers = sample names) plus its
actuation TSV (``{peak}.tsv``, carrying ``asm_chr / asm_start / asm_end`` per
sample) and writes a FASTA where each sequence is the assembly substring with
``--flank`` bp of context on each side.

The AGC region is ``{asm_chr}:{asm_start-flank}-{asm_end-1+flank}`` — asm_start
is 1-based and asm_end is exclusive (verified against the union-BED `seq`
column), and asm_chr is the PanSN contig name, which is exactly the AGC contig
key. The RE then sits at ``[flank, flank+RE_len)`` of each pulled sequence.

HG002 uses T2T ``chrN_MATERNAL`` / ``chrN_PATERNAL`` contig names; these are
mapped to the archive's PanSN form ``HG002#2#chrN`` (maternal) / ``HG002#1#chrN``
(paternal) by :func:`agc_contig`. A sample is written with its RE-only sequence
(unchanged) only when its AGC region truly can't be resolved: missing/NA asm
coords, an unrecognised non-PanSN name, or a contig absent from the archive. So
a flank run never drops a haplotype.

With ``--flank 0`` this is a passthrough copy (no AGC needed), which keeps the
Snakemake DAG uniform.

Usage:
    python add_flanks.py --fa peak.fa --tsv peak.tsv --agc all.agc \
        --flank 250 --out flanked.fa
"""
import argparse
import subprocess
import sys
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Return [(header, sequence), ...] preserving input order."""
    out: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(parts)))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        out.append((name, "".join(parts)))
    return out


def load_asm_coords(tsv: Path) -> dict[str, tuple[str, int, int]]:
    """Return {sample: (asm_chr, asm_start, asm_end)} from the actuation TSV.

    Samples whose row lacks usable asm columns are simply absent from the map.
    """
    import csv

    coords: dict[str, tuple[str, int, int]] = {}
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asm_chr = row.get("asm_chr")
            if not asm_chr or asm_chr == "NA":
                continue
            try:
                coords[row["sample"]] = (
                    asm_chr, int(row["asm_start"]), int(row["asm_end"])
                )
            except (KeyError, ValueError):
                continue
    return coords


def parse_agc_fasta(text: str) -> dict[str, str]:
    """Parse ``agc getctg`` stdout into {contig: sequence}.

    Headers are ``{contig}:{from}-{to}``; PanSN contig names carry no ':', so
    the contig is everything before the last ':'.
    """
    out: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = "".join(parts)
            name = line[1:].rsplit(":", 1)[0]
            parts = []
        else:
            parts.append(line.strip())
    if name is not None:
        out[name] = "".join(parts)
    return out


def agc_contig(sample: str, asm_chr: str) -> str | None:
    """Map a BED `asm_chr` to an AGC-queryable contig name.

    Most samples are already PanSN (`HG02886#1#CMxxxx.1`) and pass through.
    HG002 is the Q100 T2T assembly with `chrN_MATERNAL` / `chrN_PATERNAL`
    contigs; the AGC stores these PanSN-style as `HG002#2#chrN` (maternal) and
    `HG002#1#chrN` (paternal), i.e. hap1=paternal / hap2=maternal (verified
    against the archive). We key on the asm_chr's MATERNAL/PATERNAL suffix —
    the actual parental assembly — which is correct even though the BED's
    `Haplotype` 1/2 column for HG002 is flipped. Returns None for an
    unrecognised non-PanSN name.
    """
    if "#" in asm_chr:
        return asm_chr
    ind = sample.split("__")[-1].rsplit("_", 1)[0]   # 'HG002_1' -> 'HG002'
    for suffix, hap in (("_MATERNAL", "2"), ("_PATERNAL", "1")):
        if asm_chr.endswith(suffix):
            return f"{ind}#{hap}#{asm_chr[: -len(suffix)]}"
    return None


def pull_regions(agc: str, regions: list[str], threads: int = 4) -> dict[str, str]:
    """Batch ``agc getctg`` over regions; return {contig: sequence}.

    All of a peak's haplotype regions go out in one (chunked) invocation —
    getctg accepts any number of specs and returns them in request order. ``-p``
    disables prefetch (these are scattered short pulls); ``-t`` pins agc's
    thread count to the rule's allocation instead of agc's default (half the
    machine's cores), which would oversubscribe a shared node.
    """
    pulled: dict[str, str] = {}
    for k in range(0, len(regions), 400):   # chunk so argv never nears ARG_MAX
        chunk = regions[k : k + 400]
        try:
            res = subprocess.run(
                ["agc", "getctg", agc, *chunk, "-p", "-t", str(threads)],
                capture_output=True, text=True, check=True,
            )
        except (subprocess.CalledProcessError, FileNotFoundError) as ex:
            print(f"[add_flanks] WARNING: agc getctg failed for {len(chunk)} "
                  f"regions ({ex}); those stay RE-only", file=sys.stderr)
            continue
        pulled.update(parse_agc_fasta(res.stdout))
    return pulled


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fa", required=True, type=Path, help="RE-only peak FASTA")
    ap.add_argument("--tsv", required=True, type=Path, help="actuation TSV with asm coords")
    ap.add_argument("--agc", type=str, default=None, help="AGC archive (needed when --flank>0)")
    ap.add_argument("--flank", type=int, default=0, help="bp of context each side")
    ap.add_argument("--threads", type=int, default=4, help="agc -t thread count")
    ap.add_argument("--out", required=True, type=Path, help="output flanked FASTA")
    args = ap.parse_args()

    records = read_fasta(args.fa)

    if args.flank <= 0:
        # Passthrough — write the RE sequences unchanged.
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with args.out.open("w") as fh:
            for name, seq in records:
                fh.write(f">{name}\n{seq}\n")
        print(f"[add_flanks] flank=0; copied {len(records)} seqs to {args.out}",
              file=sys.stderr)
        return

    if not args.agc:
        sys.exit("[add_flanks] --flank > 0 requires --agc")

    coords = load_asm_coords(args.tsv)
    regions: list[str] = []
    region_chrom: dict[str, str] = {}   # sample -> contig key for lookup
    for name, _seq in records:
        pos = coords.get(name)
        if pos is None:
            continue
        asm_chr, s, e = pos
        contig = agc_contig(name, asm_chr)
        if contig is None:                       # unrecognised non-PanSN naming
            continue
        frm, to = max(1, s - args.flank), e - 1 + args.flank
        if to < frm:
            continue
        regions.append(f"{contig}:{frm}-{to}")
        region_chrom[name] = contig

    pulled = pull_regions(args.agc, regions, args.threads) if regions else {}

    n_flanked = 0
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        for name, seq in records:
            chrom = region_chrom.get(name)
            flanked = pulled.get(chrom) if chrom else None
            if flanked:
                fh.write(f">{name}\n{flanked.upper()}\n")
                n_flanked += 1
            else:
                fh.write(f">{name}\n{seq}\n")   # fall back to RE-only
    print(f"[add_flanks] flank={args.flank}: {n_flanked}/{len(records)} seqs "
          f"flanked from AGC ({len(records) - n_flanked} kept RE-only)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
