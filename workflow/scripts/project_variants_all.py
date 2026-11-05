#!/usr/bin/env python3
"""Single-pass hg38/CHM13 projection for a WHOLE run (Stage 2, no snakemake DAG).

Reads the consolidated haplotypes.parquet (RE seq + asm coords per haplotype),
projects every peak's variable SNV sites to CHM13 + hg38 + assembly via
`project_variants.project_peak` (AGC reference pulls), flags drivers from
phylo_mm, then runs `aggregate_hg38.aggregate` to emit the variant table +
sorted, bgzipped, tabix'd VCF + BEDs. Sequences come from the parquet + the AGC
only -- no fastas/ dir. mafft/agc run per peak, so we parallelize with a thread
pool (each worker just waits on its subprocess).

    pixi run python workflow/scripts/project_variants_all.py \
        --haplotypes results/<run>/tables/haplotypes.parquet \
        --tables-dir results/<run>/tables \
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
from aggregate_hg38 import aggregate  # noqa: E402
from common_io import read_fasta  # noqa: E402
from stage1_tables import Stage1Store  # noqa: E402
from project_variants import (  # noqa: E402
    project_peak, _COLS, _header, _unlifted_header, _UNLIFTED_COLS,
)

BEDS = ("variants.hg38.bed", "variants.chm13.bed", "elements.hg38.bed")


def _load_coords(path: Path) -> dict[str, dict]:
    lines = Path(path).read_text().splitlines()
    hdr = lines[0].split("\t")
    out = {}
    for ln in lines[1:]:
        f = ln.split("\t")
        if f:
            out[f[0]] = dict(zip(hdr, f))
    return out


def _asm_coords(sub: pl.DataFrame) -> dict[str, dict]:
    """sample -> {chr, start, end} assembly coordinates for one peak."""
    return {
        s: {"chr": c, "start": st, "end": en}
        for s, c, st, en in zip(
            sub["sample"].to_list(), sub["asm_chr"].to_list(),
            sub["asm_start"].to_list(), sub["asm_end"].to_list(),
        )
    }


def _hg38_coords(crow: dict) -> tuple[str, str, str]:
    """(chr, start, end) of the peak in hg38, or NA x3 when it does not lift."""
    if str(crow.get("exists_in_hg38", "")).upper() in ("TRUE", "1"):
        return crow["hg38_chr"], crow["hg38_start"], crow["hg38_end"]
    return "NA", "NA", "NA"


def project_peaks(peaks: list[str], parts: dict, coords: dict, store: Stage1Store,
                  outdir: Path, args) -> list[dict]:
    """Write ``{outdir}/{peak}.tsv`` (project_variants rows) for every peak and
    return the unlifted (indel-only) driver rows. A peak missing coords or an
    alignment gets a header-only TSV. A peak whose projection raises gets no TSV
    (aggregate_hg38 reports it as missing)."""
    missing_coords: list[str] = []
    missing_aligned: list[str] = []
    unlifted_all: list[dict] = []

    def do_peak(peak: str) -> None:
        sub = parts.get((peak,))
        crow = coords.get(peak)
        aln_fa = store.aligned(peak)
        out = outdir / f"{peak}.tsv"
        if sub is None or crow is None or aln_fa is None:
            if crow is None:
                missing_coords.append(peak)
            if aln_fa is None:
                missing_aligned.append(peak)
            out.write_text(_header())
            return
        re_seqs = dict(zip(sub["sample"].to_list(), sub["seq"].to_list()))
        aln_haps = read_fasta(aln_fa, strip_r=False)
        hg38 = _hg38_coords(crow)
        rows, unlifted = project_peak(
            peak, aln_haps, re_seqs, _asm_coords(sub), store.phylo_mm(peak),
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
    return unlifted_all


def compress_and_index(outdir: Path, threads: int) -> None:
    """bgzip + tabix the plain VCF and BEDs aggregate_hg38 wrote. BEDs are sorted
    first (tabix needs contiguous chroms and increasing starts)."""
    vcf_plain = outdir / "variants.hg38.vcf"
    subprocess.run(["bgzip", "-f", "-@", str(threads), str(vcf_plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_plain) + ".gz"], check=True)
    for name in BEDS:
        bed = outdir / name
        subprocess.run(
            f"sort -k1,1 -k2,2n {bed} | bgzip -@ {threads} > {bed}.gz && "
            f"tabix -f -p bed {bed}.gz && rm {bed}",
            shell=True, check=True,
        )


def write_unlifted(path: Path, rows: list[dict]) -> None:
    """Drivers that are only an indel (a gap run, no single-base substitution
    column) can't be a bed/VCF SNV; list them so they aren't silently dropped."""
    with path.open("w") as fh:
        fh.write(_unlifted_header())
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in _UNLIFTED_COLS) + "\n")


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
    store = Stage1Store(args.tables_dir)
    out = args.outdir
    with tempfile.TemporaryDirectory(prefix="proj_all_") as tmpdir:
        hg38_dir = Path(tmpdir) / "hg38"
        hg38_dir.mkdir()
        unlifted = project_peaks(peaks, parts, coords, store, hg38_dir, args)
        out.mkdir(parents=True, exist_ok=True)
        aggregate(hg38_dir, peaks, args.coords,
                  out_variants=out / "variants.hg38.tsv.gz",
                  out_vcf=out / "variants.hg38.vcf",
                  out_variants_bed=out / "variants.hg38.bed",
                  out_elements_bed=out / "elements.hg38.bed",
                  out_variants_chm13_bed=out / "variants.chm13.bed")
    compress_and_index(out, args.threads)
    unl_out = out / "drivers.unlifted.tsv"
    write_unlifted(unl_out, unlifted)
    print(f"[project_all] done -> {out}/variants.{{hg38,chm13}}.* "
          f"({len(unlifted)} indel-only drivers -> {unl_out.name})", file=sys.stderr)


if __name__ == "__main__":
    main()
