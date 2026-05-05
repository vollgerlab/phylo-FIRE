# CLAUDE.md

## Project

Snakemake pipeline (`phylo-FIRE`) for per-peak phylogenetic tree inference across HPRC
fiber-seq samples. Builds ML trees (IQ-TREE) from per-haplotype DNA sequences at consensus
FIRE peaks, then renders each tree alongside a FIRE actuation heatmap.

## Quick commands

```bash
pixi install              # first-time setup
pixi run test             # 10 test peaks via test-union-peaks.bed.gz
pixi run test-paralog     # NOTCH2NL paralog group via PARE-style preextracted dir
pixi run dry-run          # DAG preview
pixi run replot           # force re-run of the plot rule on existing trees
pixi run fmt              # snakefmt workflow/Snakefile
```

There is **no `clean` task**: snakemake handles its own state via
`--forcerun <rule>`, `--rerun-incomplete`, or `--forceall`. Don't manually
`rm -rf results/`.

## Pipeline architecture

Three Snakemake steps driven by a checkpoint:

1. `checkpoint extract_fastas` -- runs in one of two modes:

   - **Default**: one-pass stream of the bgzipped BED via `bgzip -cd | python`.
     Writes `results/fastas/{peak_id}.fa` (sequences) and `results/fastas/{peak_id}.tsv`
     (actuation values).
   - **`preextracted_dir`**: when `config["preextracted_dir"]` is set, the
     checkpoint becomes a noop wrapper that symlinks pre-extracted `.fa`+`.tsv`
     pairs from there into `results/fastas/`. This is the contract for
     ingesting PARE's per-paralog-group output (see "Cross-repo: PARE" below).
     `config["input_table"]` is ignored when `preextracted_dir` is set.

   In both cases the checkpoint sentinel is `results/fastas/.done`; downstream
   rules are discovered via `glob_wildcards("results/fastas/{peak_id}.fa")`
   after it completes.

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
        preextracted_pg_9/        PARE-emitted NOTCH2 paralog group (pg_9.fa + pg_9.tsv)
profiles/
    slurm/config.yaml         SLURM executor profile (CHPC Redwood)
results/                      all pipeline outputs (gitignored)
    fastas/  {peak_id}.fa  {peak_id}.tsv
    aligned/ {peak_id}.fa  (temp, deleted after plot)
    trees/   {peak_id}.treefile
    plots/   {peak_id}.html  {peak_id}.pdf
```

## Cross-repo: PARE

The PARE pipeline (https://github.com/vollgerlab/PARE) lives at
`/scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/PARE/` and emits per-
paralog-group FASTA + actuation TSV under `results/paralog_msas/` in
exactly the on-disk format `extract_fastas` would produce. To run a tree
per paralog group instead of per peak:

```bash
pixi run -- snakemake \
    --config preextracted_dir=/scratch/.../PARE/results/paralog_msas \
    -j 32
```

The `pg_9.fa` test fixture (`resources/test/preextracted_pg_9/`) is
PARE's chr1 NOTCH2 paralog group: NOTCH2 + NOTCH2NLA + NOTCH2NLB +
NOTCH2NLC, 276 sequences, 71 unique. Headers carry `{peak_id}__{gene}__{ind}_{hap}`
(e.g. `chr1_145345837_145346088__NOTCH2NLA__HG01891_1`) so leaf labels
in the rendered tree are paralog-readable. PARE prefers the
protein-coding gene name over alternate-TSS LOC* lncRNAs in its
gene-annotation pipeline; see PARE's `CLAUDE.md` for the full story.

`pixi run test-paralog` exercises the full pipeline on this fixture
(~90 s wallclock), producing `results/{trees,plots}/pg_9.{treefile,html,pdf}`.

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
