# CLAUDE.md

## Project

Snakemake pipeline (`phylo-FIRE`) for per-peak phylogenetic tree inference across HPRC
fiber-seq samples. Builds ML trees (IQ-TREE) from per-haplotype DNA sequences at consensus
FIRE peaks, then renders each tree alongside a FIRE actuation heatmap.

## Quick commands

```bash
pixi install              # first-time setup
pixi run test             # clean run on 10 test peaks (deletes results/ first)
pixi run dry-run          # DAG preview
pixi run fmt              # snakefmt workflow/Snakefile
```

## Pipeline architecture

Three Snakemake steps driven by a checkpoint:

1. `checkpoint extract_fastas` -- one-pass stream of the bgzipped BED via `bgzip -cd | python`.
   Writes `results/fastas/{peak_id}.fa` (sequences) and `results/fastas/{peak_id}.tsv` (actuation values).
   The checkpoint sentinel is `results/fastas/.done`; downstream rules are discovered via
   `glob_wildcards("results/fastas/{peak_id}.fa")` after it completes.

2. `rule mafft` -- per-peak alignment (`--auto --maxiterate 1000 --reorder --adjustdirection`).
   Output is `temp()` (deleted after plot finishes). `--adjustdirection` adds `_R_` prefix to
   reversed sequences; the plot script strips this prefix from both alignment and tree leaf names.

3. `rule iqtree` -- per-peak ML tree. Uses `shadow: "minimal"` so IQ-TREE's auxiliary
   output files (`.log`, `.iqtree`, `.bionj`, `.mldist`, `.model.gz`) are discarded;
   only the `.treefile` survives. Bootstrap (`-B 1000`) is conditional: enabled only when
   the aligned FASTA has >=4 distinct sequences (IQ-TREE refuses bootstrap otherwise).

4. `rule plot` -- produces `results/plots/{peak_id}.html` and `results/plots/{peak_id}.pdf`
   via `workflow/scripts/plot_tree_actuation_html.py`: midpoint-rooted, ladderized tree on the
   left; FIRE actuation heatmap in the center; MSA with per-cell consensus coloring on the right.

## Key implementation decisions

- **`iqtree` not `iqtree2`**: the conda-forge package installs as `iqtree` (version 3+).
- **Bootstrap gated on unique seq count**: `count_unique_seqs()` parses the multi-line FASTA
  properly. Counting `>` headers is wrong because IQ-TREE deduplicates internally.
- **`min_seqs: 2`** (config default): allows peaks with only 2 distinct haplotype sequences.
  With `min_seqs: 3` (old default), some valid peaks were silently skipped.
- **Actuation = `fire_cov_OV_cov`** (column 9, 0-indexed): `fire_coverage / coverage`,
  a 0-1 fraction. Stored in `fastas/{peak_id}.tsv` alongside the FASTA.
- **`haplotype_distinct` is a peak-level filter**: a peak is included if any primary-sample
  row has the target classification. All qualifying sequences go into the tree regardless.
- **Streaming extraction**: `extract_fastas.py` processes one peak at a time (constant memory),
  relying on the BED being sorted by genomic coordinate so rows for a peak are consecutive.

## File layout

```
pixi.toml
config/
    config.yaml               pipeline settings
    test_peak_ids.txt         10 peak IDs for test runs
workflow/
    Snakefile                 main pipeline
    scripts/
        extract_fastas.py         one-pass FASTA + TSV extraction
        plot_tree_actuation_html.py  interactive HTML + PDF tree/MSA/actuation plot
    rules/                    (reserved for future rule modules)
resources/
    test/
        test-union-peaks.bed.gz   900-row subset (10 peaks)
profiles/
    slurm/config.yaml         SLURM executor profile (CHPC Redwood)
results/                      all pipeline outputs (gitignored)
    fastas/  {peak_id}.fa  {peak_id}.tsv
    aligned/ {peak_id}.fa  (temp, deleted after plot)
    trees/   {peak_id}.treefile
    plots/   {peak_id}.html  {peak_id}.pdf
```

## Input table columns referenced

| Column | Index | Used for |
|--------|-------|----------|
| `consensus_peak_id` | 4 | grouping key, output filename |
| `fire_cov_OV_cov` | 9 | actuation fraction for heatmap |
| `seq` | 23 | DNA sequence (skip if NA) |
| `Individual_ID` | 25 | FASTA header / tree leaf prefix |
| `Haplotype` | 26 | FASTA header / tree leaf suffix |
| `is_primary_sample` | 36 | filter when `primary_only: true` |
| `haplotype_distinct` | 37 | optional peak-level filter |

## SLURM notes

Profile at `profiles/slurm/config.yaml` targets CHPC Redwood (`redwood15.chpc.utah.edu`).
Verify `slurm_partition` and `slurm_account` values with `sinfo` before the first cluster run.
