# CLAUDE.md

## Project

Snakemake pipeline (`phylo-FIRE`) for per-peak phylogenetic tree inference across HPRC
fiber-seq samples. Builds ML trees (IQ-TREE) from per-haplotype DNA sequences at consensus
FIRE peaks, then renders each tree alongside a FIRE actuation heatmap.

## Quick commands

```bash
pixi install              # first-time setup
pixi run setup-pclai      # one-time: bulk-fetch HPRC R2 PCLAI BEDs → resources/pclai/hprc_r2_asm.bed.gz
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
   via `workflow/scripts/plot_tree_actuation_html.py`. Layout left → right:
   midpoint-rooted ladderized tree, leaf labels, **PCLAI ancestry stripe**, FIRE
   actuation heatmap, MSA with per-cell consensus coloring (toggleable).
   The PCLAI column is shown only when `config["pclai_bed"]` points at a valid
   bgzipped + tabix-indexed BED; otherwise it collapses to zero width.

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
- **PCLAI lookup is CHM13-coord** (not asm). The plot script parses `peak_id` as
  `{chrom}_{start}_{end}` and runs **one tabix query per peak**; each returned row's
  source sample/hap is encoded in BED col 4 (`HG00097/h1/...`). This wins over
  per-leaf asm-coord lookup because: (1) one query vs. N queries, (2) works for
  `preextracted_dir` mode (PARE) without any upstream changes, (3) no PanSN/accession
  heterogeneity. Tradeoff: HG002-only consensus peaks like `chr8_MATERNAL_8587817_8588089`
  silently fall through (chrom doesn't exist in CHM13 PCLAI) — small fraction in practice.
  `extract_fastas.py` still emits `asm_chr`/`asm_start`/`asm_end` in the per-peak `.tsv`
  for forward-compat / debugging; the plot ignores them.

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
    pclai/                        gitignored: built by `pixi run setup-pclai`
        hprc_r2_chm13.bed.gz      merged HPRC R2 PCLAI chm13-coord BED (~126 MB, 11.7M rows)
        hprc_r2_chm13.bed.gz.tbi
        raw_chm13/                per-sample-hap source BEDs (cache; safe to delete)
    test/
        test-union-peaks.bed.gz             900-row subset (10 peaks)
        preextracted_pg_31/                 small NOTCH2NL promoter fixture (24 seqs, 6 haps)
        preextracted_pg_31_all_haps/        full-cohort version (319 seqs, 78 haps)
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

The `pg_31.fa` test fixture (`resources/test/preextracted_pg_31/`) is
PARE's canonical chr1 NOTCH2NL promoter family: NOTCH2 (1p11) +
NOTCH2NLA / NOTCH2NLB / NOTCH2NLC promoters (1q21) + a 5th paralog at
chr1:120,736,956-120,737,187. The 5th cp's nearest_pc gene name is
"NBPF26" because its actual gene (NOTCH2NLR) is annotated as a
pseudogene/lncRNA in the RefSeq Liftoff v5.2 track and gets skipped by
the protein-coding filter; sequence identity to the other four is
99.5-99.7%, the same level the four "real" paralogs share with each
other (NBPF genes carry DUF1220/Olduvai domains and are NOT
sequence-similar to NOTCH2 family). Treat it as a NOTCH2NLR proxy.

5 consensus peaks x 6 primary haplotype-samples = 24 sequences in the
small fixture (NBPF26/NOTCH2NLR coverage is sparser than the four
named paralogs). Headers carry `{peak_id}__{gene}__{ind}_{hap}`
(e.g. `chr1_145345837_145346088__NOTCH2NLA__HG01891_1`) so leaf labels
in the rendered tree are paralog-readable.

The fixture is *small by design* -- one homologous regulatory element
across the paralog loci, not a tile of every regulatory element in
the NOTCH2NL landscape. PARE's `target_peak_coverage` filter (80%
threshold on the lifted projection's overlap with the target peak's
interval) prevents N x M edge proliferation between paralog loci that
each contain multiple tiling peaks; without it, a previous incarnation
of this fixture (pg_9) had 46 cps and 276 sequences from collapsing
13 distinct regulatory elements across 4 paralogs into one giant
group. The fixture was renamed pg_9 -> pg_64 -> pg_31 over successive
PARE runs because PARE's union-find paralog-group IDs aren't stable
between runs (any change to the upstream filter cascade or graph
ordering reshuffles them); pg_31 is the chr1 SD belt rerun that
picked up the 5th NOTCH2NLR-like paralog.

`pixi run test-paralog` exercises the full pipeline on this fixture
(~30 s wallclock), producing `results/{trees,plots}/pg_31.{treefile,html,pdf}`.

A second fixture `resources/test/preextracted_pg_31_all_haps/` carries
the *same* 5 NOTCH2NL peaks but pulled from the full HPRCv2 union BED
(45 individuals; 78 primary haps after `is_primary_sample` filter;
319 sequences total once NA-seq haps are dropped). The fixture file
basenames are deliberately `pg_31_all_haps.{fa,tsv}` (not `pg_31`) so
the all-haps build artifacts at `results/{trees,plots}/pg_31_all_haps.{...}`
do not overwrite the small `pg_31.{...}` outputs from `test-paralog`.
Run with `pixi run test-paralog-all-haps`. Useful for validating that
the tree and heatmap renderers scale to the production-shape input.

## Input table columns referenced

| Column | Index | Used for |
|--------|-------|----------|
| `consensus_peak_id` | 4 | grouping key, output filename |
| `fire_cov_OV_cov` | 9 | actuation fraction for heatmap |
| `asm_chr` | 11 | PanSN per-sample asm contig (e.g. `HG02886#1#CMxxxx.1`); used for PCLAI lookup |
| `asm_start` | 12 | peak start in asm coords (PCLAI tabix range) |
| `asm_end` | 13 | peak end in asm coords (PCLAI tabix range) |
| `seq` | 23 | DNA sequence (skip if NA) |
| `Individual_ID` | 25 | FASTA header / tree leaf prefix |
| `Haplotype` | 26 | FASTA header / tree leaf suffix; matches PCLAI's `#N#` directly |
| `is_primary_sample` | 36 | filter when `primary_only: true` |
| `haplotype_distinct` | 37 | optional peak-level filter |

## PCLAI ancestry overlay

`pixi run setup-pclai` (one-time) bulk-downloads the HPRC R2 PCLAI v1.1
**chm13-coord** BEDs (463 sample-haps including CHM13 itself, source index at
`hprc_intermediate_assembly/data_tables/annotation/pclai/`), concatenates them,
sorts, bgzips, and tabix-indexes the merged file. Output:
`resources/pclai/hprc_r2_chm13.bed.gz` (~126 MB, 11.7M rows, gitignored).

Each leaf's stripe color is HPRC's pre-baked `itemRgb` for the PCLAI window
with the largest overlap of the peak interval — same palette HPRC uses, so
visualizations stay consistent across tools. Lookup is one tabix query per
peak; rows are then matched to leaves via the `{sample}/h{hap}` prefix in
BED col 4.

The merged BED column layout (BED9+1):
- col 1: CHM13 reference chrom (`chr1` … `chrY`, `chrM`)
- col 4: `{sample}/h{hap}/{chrom}_w{NNNN}_(pc1,pc2)` — sample/hap identity
  and per-window PC coordinates embedded here
- col 5: confidence score (higher = better)
- col 9: itemRgb for the PCLAI palette
- col 10: centroid PC `(pc1,pc2)` — one of 4 buckets

Why CHM13 over asm coords: phylo-FIRE peak_ids are CHM13-coord-native, so a
single tabix query per peak returns one row per sample-haplotype that overlaps
the locus — no per-leaf assembly lookups, no PanSN/accession naming
heterogeneity, no PARE changes needed for paralog mode.

CHM13-only peaks work cleanly. HG002-only consensus peaks
(e.g. `chr8_MATERNAL_8587817_8588089`) parse to a chrom that doesn't exist in
CHM13 PCLAI and silently fall through to no stripe. If you ever need PCLAI
on those, the asm-coord variant of the BED handles them — point `pclai_bed`
at a different file and update the plot script's lookup accordingly.

## SLURM notes

Profile at `profiles/slurm/config.yaml` targets CHPC Redwood (`redwood15.chpc.utah.edu`).
Verify `slurm_partition` and `slurm_account` values with `sinfo` before the first cluster run.
