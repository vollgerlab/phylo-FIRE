# phylo-FIRE

Snakemake pipeline for building per-peak phylogenetic trees across HPRC fiber-seq samples.
For each consensus FIRE peak, sequences are extracted from all 90 haplotypes (45 individuals),
aligned with MAFFT, and used to infer a maximum-likelihood tree with IQ-TREE. Each tree is
exported alongside a FIRE actuation heatmap showing per-sample accessibility.

## Pipeline overview

```
union-peaks BED (bgzipped, ~3 GB)
        |
        v
[checkpoint] extract_fastas.py   -- one-pass stream; writes fastas/{peak_id}.fa
                                     and fastas/{peak_id}.tsv (per-sample actuation)
        |
        v (per peak, parallel)
rule mafft        aligned/{peak_id}.fa     -- MAFFT --auto --maxiterate 1000
rule iqtree       trees/{peak_id}.treefile -- IQ-TREE ModelFinder + bootstrap (if >=4 unique seqs)
rule plot         plots/{peak_id}.pdf      -- tree + FIRE actuation heatmap
```

Intermediate aligned FASTAs are marked `temp()` and deleted after IQ-TREE completes.
IQ-TREE auxiliary files (`.log`, `.iqtree`, `.bionj`, etc.) are isolated in a shadow
directory and discarded; only the `.treefile` is kept.

## Quick start

```bash
pixi install              # set up the environment
pixi run dry-run          # preview the DAG (no execution)
pixi run test             # clean build + run on 10 test peaks
pixi run unit             # pytest unit tests for workflow/scripts
```

## Running on the full dataset

```bash
# All haplotype-selective peaks
snakemake --profile profiles/slurm \
    --config haplotype_distinct="Haplotype-selective" -j 200

# A specific list of peaks
snakemake --profile profiles/slurm \
    --config peak_ids=my_peaks.txt -j 200

# Both filters combined
snakemake --profile profiles/slurm \
    --config peak_ids=my_peaks.txt haplotype_distinct="Haplotype-selective" -j 200
```

## Outputs

| Path | Description |
|------|-------------|
| `fastas/{peak_id}.fa` | Unaligned FASTA (one entry per primary sample/haplotype) |
| `fastas/{peak_id}.tsv` | Per-sample FIRE actuation fraction (`fire_cov_OV_cov`, 0-1) |
| `trees/{peak_id}.treefile` | Newick ML tree from IQ-TREE |
| `plots/{peak_id}.pdf` | Tree + actuation heatmap (white-to-dark-red, 10% bins) |

## Configuration (`config.yaml`)

| Key | Default | Description |
|-----|---------|-------------|
| `input_table` | `Tables/union-peaks-cons.bed.gz` | bgzipped union-peaks BED |
| `peak_ids` | `null` | Path to a file of consensus_peak_ids (one per line) |
| `haplotype_distinct` | `null` | Filter peaks by classification (e.g. `Haplotype-selective`) |
| `primary_only` | `true` | Restrict to primary samples |
| `min_seqs` | `2` | Minimum unique sequences to build a tree |
| `mafft_maxiterate` | `1000` | MAFFT iterative refinement cycles |
| `bgzip_threads` | `16` | Threads for bgzip decompression |

## Input table

The input is the bgzipped BED produced by the `union-peaks.Rmd` pipeline. Key columns used:

| Column | Index | Description |
|--------|-------|-------------|
| `consensus_peak_id` | 4 | Unique peak identifier (e.g. `chr8_8453862_8454062`) |
| `fire_cov_OV_cov` | 9 | FIRE actuation fraction (fire_coverage / coverage, 0-1) |
| `seq` | 23 | DNA sequence at this peak for this haplotype (NA if unmapped) |
| `Individual_ID` | 25 | Sample individual (e.g. `HG00438`) |
| `Haplotype` | 26 | `1` or `2` |
| `is_primary_sample` | 36 | `TRUE` for the preferred sample per individual |
| `haplotype_distinct` | 37 | Peak classification: Shared, Haplotype-selective, etc. |

Peak IDs containing `#` (assembly-coordinate peaks) are sanitized to `_` for filenames.

## Test dataset

`test/test-union-peaks.bed.gz` contains 810 rows covering 9 peaks (90 haplotypes each):

- 4 Haplotype-selective (genetic variants between haplotypes)
- 3 Epigenetically haplotype-selective
- 3 Shared (expect similar sequences; useful as controls)
- `chr8_70445067_70445271`: two-driver exemplar (main-figure example). RE112 (β≈−0.33)
  closes the element on two haplotype backgrounds; RE63 (β≈+0.17) opens it on a separate clade.

`test_peak_ids.txt` lists all 9 peak IDs. All 9 produce trees; peaks with fewer than 4
unique sequences run without ultrafast bootstrap.

## Data locations (relative to `/mmfs1/gscratch/stergachislab/HPRC/mrv-figures/`)

| File | Size | Description |
|------|------|-------------|
| `Tables/union-peaks-cons.bed.gz` | ~3 GB | Full union-peaks table (~65 M rows) |
| `Tables/peaks-cons.bed.gz` | ~700 MB | Peaks-only subset (~7.7 M rows) |
| `phylo-trees/test/test-union-peaks.bed.gz` | small | 900-row test subset |

## Dependencies

Managed via `pixi.toml` (conda-forge/bioconda):

- `bgzip` / `htslib` (decompression)
- `mafft` (multiple sequence alignment)
- `iqtree` (ML tree inference; installs as `iqtree`, version 3+)
- `biopython` (tree parsing and drawing)
- `matplotlib` (visualization)
- `snakemake >= 8.0` with `snakemake-executor-plugin-slurm`
