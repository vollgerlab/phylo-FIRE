# Phylogenetic Tree Pipeline for Consensus Peaks

## Overview
Build phylogenetic trees for consensus peaks across all HPRC fiber-seq samples.
For each consensus peak, collect the DNA sequence from every sample/haplotype,
align with MAFFT, and infer a tree with IQ-TREE.

## Pipeline Steps

### Step 0: Extract FASTAs (one-pass)
- Read the union peaks table once with `bgzip -cd -@ 16`
- Filter to primary samples with non-NA sequences
- For each consensus_peak_id, write a FASTA: `>Individual_ID_Haplotype\nsequence`
- Sanitize peak IDs for filenames (`#` → `_`)
- Skip peaks with <3 unique sequences (no meaningful tree)
- Output: `fastas/{peak_id}.fa`

### Step 1: MAFFT Alignment
- `mafft --maxiterate 1000 fastas/{peak_id}.fa > aligned/{peak_id}.fa`
- For simple/many sequences: `--retree 2` (faster progressive method)
- Output: `aligned/{peak_id}.fa`

### Step 2: IQ-TREE
- `iqtree -s aligned/{peak_id}.fa`
- Uses ModelFinder automatically for best-fit model
- Output: `trees/{peak_id}.fa.treefile` (NEWICK format)

## Directory Structure
```
phylo-trees/
├── Snakefile
├── README.md
├── config.yaml
├── test_peak_ids.txt          # 10 test peaks
├── test/
│   └── test-union-peaks.bed.gz  # subset for testing
├── scripts/
│   └── extract_fastas.py      # one-pass FASTA extraction
├── fastas/                    # per-peak FASTA files
├── aligned/                   # MAFFT aligned FASTAs
└── trees/                     # IQ-TREE output
```

## Test Dataset
`test_peak_ids.txt` contains 10 consensus peaks:
- 4 Haplotype-selective (genetic variants between haplotypes)
- 3 Epigenetically haplotype-selective
- 3 Shared (control, expect similar sequences)

## Usage
```bash
# Test run (10 peaks)
snakemake --config peak_ids=test_peak_ids.txt input_table=test/test-union-peaks.bed.gz -j 4

# Full run (all haplotype-selective peaks)
snakemake --config haplotype_distinct=Haplotype-selective -j 50 --profile profiles/slurm

# Specific peaks from a file
snakemake --config peak_ids=my_peaks.txt -j 50 --profile profiles/slurm
```

## Input Table Structure
The input is a bgzipped BED from the `union-peaks.Rmd` pipeline (see `Tables/union-peaks-README.md`
for full column docs). Key columns for this pipeline:

| Column | Description |
|--------|-------------|
| consensus_peak_id | Unique ID per consensus peak (e.g. `chr8_8453862_8454062` or `HG00438#1#CM089167.1_120346427_120346657`) |
| Individual_ID | Individual (e.g. `HG00438`) |
| Haplotype | `1` or `2` |
| seq | DNA sequence at this peak for this sample/haplotype (can be NA) |
| haplotype_distinct | Classification: Shared, Haplotype-selective, Epigenetically haplotype-selective, Genetically haplotype-specific, Complex, chrX, chrY, Untestable |
| is_primary_sample | TRUE for the preferred sample per individual (PacBio > ONT; deep HG002 > shallow) |
| is_peak | Whether this site passes peak thresholds in this sample |

- 90 haplotypes from 45 individuals (39 primary, 6 non-primary duplicates)
- Each consensus_peak_id has one row per sample-haplotype (~90 rows)
- The `seq` column can be NA when the peak has no assembly mapping for that sample

## File Locations (relative to repo root `/mmfs1/gscratch/stergachislab/HPRC/mrv-figures/`)
- Full union table: `Tables/union-peaks-cons.bed.gz` (~3GB, ~65M rows)
- Peak-only table: `Tables/peaks-cons.bed.gz` (~700MB, ~7.7M rows)
- This pipeline: `phylo-trees/`
- Test subset: `phylo-trees/test/test-union-peaks.bed.gz` (900 rows, 10 peaks)
- Extraction script (partially written): `phylo-trees/scripts/extract_fastas.py`

## Conventions
- Use `bgzip -cd -@ 16` for decompression (fast, multithreaded)
- Use pixi for environment management (see `pixi.toml` at repo root)
- Snakemake: use `rule all` with explicit outputs, `temp()` for intermediates, parameterize threads via `resources:`
- Filename sanitization: replace `#` with `_` in consensus_peak_ids for filesystem safety

## Dependencies
- bgzip/tabix (htslib)
- MAFFT (alignment)
- IQ-TREE (tree inference)
- Python 3 (FASTA extraction)
