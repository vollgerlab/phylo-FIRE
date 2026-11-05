#!/usr/bin/env python3
"""Single-pass hg38/CHM13 projection for a WHOLE run (Stage 2, no snakemake DAG).

Reads the consolidated haplotypes.parquet (RE seq + asm coords per haplotype),
projects every peak's variable SNV sites to CHM13 + hg38 + assembly via
`project_variants.project_peak` (AGC reference pulls), flags drivers from
phylo_mm, then reuses `aggregate_hg38.py` to emit the variant table + sorted,
bgzipped, tabix'd VCF + BEDs. Sequences come from the parquet + the AGC only --
no fastas/ dir. mafft/agc run per peak, so we parallelize with a thread pool
(each worker just waits on its subprocess).

    pixi run python workflow/scripts/project_variants_all.py \
        --haplotypes results/<run>/tables/haplotypes.parquet \
        --phylo-mm-dir results/<run>/phylo_mm \
        --coords results/<run>/tables/peak_coords.tsv \
        --peaks results/<run>/fastas/peak_ids.txt \
        --agc resources/HPRC_r2_assemblies_0.6.1.agc \
        --outdir results/<run>/tables --threads 16
"""
import argparse
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage1_tables import Stage1Store  # noqa: E402
from common_io import read_fasta  # noqa: E402
from project_variants import (  # noqa: E402
    project_peak, _COLS, _header, _unlifted_header, _UNLIFTED_COLS,
)


def _load_coords(path: Path) -> dict[str, dict]:
    lines = Path(path).read_text().splitlines()
    hdr = lines[0].split("\t")
    out = {}
    for ln in lines[1:]:
        f = ln.split("\t")
        if f:
            out[f[0]] = dict(zip(hdr, f))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--haplotypes", required=True, type=Path)
    ap.add_argument("--tables-dir", required=True, type=Path,
                    help="results/<run>/tables (Stage-1 summary parquets)")
    ap.add_argument("--coords", required=True, type=Path, help="peak_coords.tsv")
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--agc", required=True)
    ap.add_argument("--flank", type=int, default=0)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--pad", type=int, default=25)
    ap.add_argument("--min-carriers", type=int, default=2)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    coords = _load_coords(args.coords)
    parts = pl.read_parquet(args.haplotypes).partition_by("peak_id", as_dict=True)

    with tempfile.TemporaryDirectory(prefix="proj_all_") as tmpdir:
        tmp = Path(tmpdir)
        hg38_dir = tmp / "hg38"
        hg38_dir.mkdir()
        store = Stage1Store(args.tables_dir)
        missing_coords: list[str] = []
        missing_aligned: list[str] = []
        unlifted_all: list[dict] = []

        def do_peak(peak: str) -> None:
            sub = parts.get((peak,))
            crow = coords.get(peak)
            aln_fa = store.aligned(peak)
            out = hg38_dir / f"{peak}.tsv"
            if sub is None or crow is None or aln_fa is None:
                if crow is None:
                    missing_coords.append(peak)
                if aln_fa is None:
                    missing_aligned.append(peak)
                out.write_text(_header())
                return
            re_seqs = dict(zip(sub["sample"].to_list(), sub["seq"].to_list()))
            aln_haps = read_fasta(aln_fa, strip_r=False)
            asm = {
                s: {"chr": c, "start": st, "end": en}
                for s, c, st, en in zip(
                    sub["sample"].to_list(), sub["asm_chr"].to_list(),
                    sub["asm_start"].to_list(), sub["asm_end"].to_list(),
                )
            }
            hg38 = ("NA", "NA", "NA")
            if str(crow.get("exists_in_hg38", "")).upper() in ("TRUE", "1"):
                hg38 = (crow["hg38_chr"], crow["hg38_start"], crow["hg38_end"])
            rows, unlifted = project_peak(
                peak, aln_haps, re_seqs, asm, store.phylo_mm(peak),
                crow["cons_chr"], int(crow["cons_start"]), int(crow["cons_end"]),
                hg38[0], hg38[1], hg38[2], args.agc,
                flank=args.flank, threads=1, pad=args.pad, min_carriers=args.min_carriers,
            )
            with out.open("w") as fh:
                fh.write(_header())
                for r in rows:
                    fh.write("\t".join(str(r[c]) for c in _COLS) + "\n")
            unlifted_all.extend(unlifted)

        print(f"[project_all] projecting {len(peaks)} peaks with {args.threads} "
              f"workers ...", file=sys.stderr)
        done = 0
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            for _ in as_completed([ex.submit(do_peak, p) for p in peaks]):
                done += 1
                if done % 1000 == 0 or done == len(peaks):
                    print(f"[project_all]   {done}/{len(peaks)} peaks projected",
                          file=sys.stderr)
        if missing_coords:
            print(f"[project_all] {len(missing_coords)} peaks not in coords table "
                  f"(skipped): {', '.join(missing_coords[:10])}", file=sys.stderr)
        if missing_aligned:
            print(f"[project_all] {len(missing_aligned)} peaks missing aligned/*.fa "
                  f"(skipped): {', '.join(missing_aligned[:10])}", file=sys.stderr)

        # Reuse aggregate_hg38 (variant table + plain VCF + BEDs), then bgzip + tabix the VCF.
        args.outdir.mkdir(parents=True, exist_ok=True)
        vcf_plain = args.outdir / "variants.hg38.vcf"
        agg = Path(__file__).resolve().parent / "aggregate_hg38.py"
        subprocess.run(
            ["python", str(agg), "--indir", str(hg38_dir), "--peaks", str(args.peaks),
             "--coords", str(args.coords),
             "--out-variants", str(args.outdir / "variants.hg38.tsv.gz"),
             "--out-vcf", str(vcf_plain),
             "--out-variants-bed", str(args.outdir / "variants.hg38.bed"),
             "--out-variants-chm13-bed", str(args.outdir / "variants.chm13.bed"),
             "--out-elements-bed", str(args.outdir / "elements.hg38.bed")],
            check=True,
        )
    subprocess.run(["bgzip", "-f", "-@", str(args.threads), str(vcf_plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_plain) + ".gz"], check=True)
    # BEDs: sort (tabix needs contiguous chroms, increasing starts), bgzip, index
    for name in ("variants.hg38.bed", "variants.chm13.bed", "elements.hg38.bed"):
        bed = args.outdir / name
        subprocess.run(
            f"sort -k1,1 -k2,2n {bed} | bgzip -@ {args.threads} > {bed}.gz && "
            f"tabix -f -p bed {bed}.gz && rm {bed}",
            shell=True, check=True,
        )

    # Drivers that are only an indel (a gap run, no single-base substitution
    # column) can't be a bed/VCF SNV; list them so they aren't silently dropped.
    unl_out = args.outdir / "drivers.unlifted.tsv"
    with unl_out.open("w") as fh:
        fh.write(_unlifted_header())
        for r in unlifted_all:
            fh.write("\t".join(str(r[c]) for c in _UNLIFTED_COLS) + "\n")
    print(f"[project_all] done -> {args.outdir}/variants.{{hg38,chm13}}.* "
          f"({len(unlifted_all)} indel-only drivers -> {unl_out.name})", file=sys.stderr)


if __name__ == "__main__":
    main()
