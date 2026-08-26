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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from project_variants import project_peak, _COLS, _header  # noqa: E402


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
    ap.add_argument("--phylo-mm-dir", required=True, type=Path)
    ap.add_argument("--coords", required=True, type=Path, help="peak_coords.tsv")
    ap.add_argument("--peaks", required=True, type=Path, help="run manifest")
    ap.add_argument("--agc", required=True)
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--pad", type=int, default=25)
    ap.add_argument("--min-carriers", type=int, default=2)
    args = ap.parse_args()

    peaks = [p.strip() for p in args.peaks.read_text().splitlines() if p.strip()]
    coords = _load_coords(args.coords)
    parts = pl.read_parquet(args.haplotypes).partition_by("peak_id", as_dict=True)

    tmp = Path(tempfile.mkdtemp(prefix="proj_all_"))
    hg38_dir = tmp / "hg38"
    hg38_dir.mkdir()
    missing_coords: list[str] = []

    def do_peak(peak: str) -> None:
        sub = parts.get((peak,))
        crow = coords.get(peak)
        out = hg38_dir / f"{peak}.tsv"
        if sub is None or crow is None:
            if crow is None:
                missing_coords.append(peak)
            out.write_text(_header())
            return
        haps = dict(zip(sub["sample"].to_list(), sub["seq"].to_list()))
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
        pmm = args.phylo_mm_dir / f"{peak}.tsv"
        rows = project_peak(
            peak, haps, asm, pmm if pmm.exists() else None,
            crow["cons_chr"], int(crow["cons_start"]), int(crow["cons_end"]),
            hg38[0], hg38[1], hg38[2], args.agc,
            threads=1, pad=args.pad, min_carriers=args.min_carriers,
        )
        with out.open("w") as fh:
            fh.write(_header())
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in _COLS) + "\n")

    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        list(ex.map(do_peak, peaks))
    if missing_coords:
        print(f"[project_all] {len(missing_coords)} peaks not in coords table "
              f"(skipped): {', '.join(missing_coords[:10])}", file=sys.stderr)

    # Reuse aggregate_hg38 (writes the variant table + plain VCF + BEDs), then
    # bgzip + tabix the VCF (matching the old aggregate_hg38 rule shell).
    args.outdir.mkdir(parents=True, exist_ok=True)
    vcf_plain = args.outdir / "variants.hg38.vcf"
    agg = Path(__file__).resolve().parent / "aggregate_hg38.py"
    subprocess.run(
        ["python", str(agg), "--indir", str(hg38_dir), "--peaks", str(args.peaks),
         "--coords", str(args.coords),
         "--out-variants", str(args.outdir / "variants.hg38.tsv.gz"),
         "--out-vcf", str(vcf_plain),
         "--out-variants-bed", str(args.outdir / "variants.hg38.bed"),
         "--out-elements-bed", str(args.outdir / "elements.hg38.bed")],
        check=True,
    )
    subprocess.run(["bgzip", "-f", "-@", str(args.threads), str(vcf_plain)], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_plain) + ".gz"], check=True)
    print(f"[project_all] done -> {args.outdir}/variants.hg38.*", file=sys.stderr)


if __name__ == "__main__":
    main()
