# CLAUDE.md

## Project

Snakemake pipeline (`phylo-FIRE`) for per-peak phylogenetic tree inference across HPRC
fiber-seq samples. Builds ML trees (IQ-TREE) from per-haplotype DNA sequences at consensus
FIRE peaks, then renders each tree alongside a FIRE actuation heatmap.

## Quick commands

```bash
pixi install                       # first-time setup
pixi run setup-pclai               # one-time: bulk-fetch HPRC R2 PCLAI BEDs → resources/pclai/hprc_r2_asm.bed.gz
pixi run setup-selectiveinference  # one-time: install selectiveInference R pkg from CRAN (not on conda-forge)
pixi run test                      # 8 test peaks via test-union-peaks.bed.gz
pixi run test-paralog              # NOTCH2NL paralog group via PARE-style preextracted dir
pixi run hap-selective-5           # 5 hap-selective peaks pulled from hyak via tabix-over-ssh
pixi run dry-run                   # DAG preview
pixi run replot                    # force re-run of the plot rule on existing trees
pixi run fmt                       # snakefmt workflow/Snakefile

# Cluster (CHPC) production run — uses the asm-coord union BED + AGC paths:
pixi run -- snakemake --configfile config/cluster.yaml \
    --config peak_ids=<peaks.txt> --profile profiles/slurm -j <N>
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

1b. `rule add_flanks` (only when `config["flank"] > 0`) -- pulls each haplotype's
   sequence with `flank` bp of assembly context each side from the AGC archive
   (`config["agc"]`, e.g. `resources/HPRC_r2_assemblies_0.6.1.agc`), keyed by the
   PanSN `asm_chr` in the per-peak TSV. Region is
   `{asm_chr}:{asm_start-flank}-{asm_end-1+flank}` (asm_start 1-based, asm_end
   exclusive; verified against the BED `seq` column). One batched `agc getctg`
   call per peak over all its haplotype regions (`-p`, `-t {threads}`). Rows
   whose contig isn't PanSN/in the archive (e.g. HG002 `chrN_PATERNAL`) fall
   back to RE-only, so a flank run never drops a haplotype. The RE then sits at
   `[flank, flank+RE_len)` of each pulled sequence. `_mafft_input` routes mafft
   to `results/flanked/{peak}.fa` when flanking is on, RE-only otherwise, so the
   checkpoint's peak manifest stays independent of flank. Flanks add variable
   sites that resolve the per-peak tree and reduce H²ₚ prior-dependence;
   `preextracted_dir`/PARE controls its own flanks (this rule doesn't apply).

2. `rule mafft` -- per-peak alignment (`--auto --maxiterate 1000 --reorder --adjustdirection`).
   Output is `temp()` (deleted after plot finishes). `--adjustdirection` adds `_R_` prefix to
   reversed sequences; the plot script strips this prefix from both alignment and tree leaf names.

3. `rule encode_indels` -- Simmons-Ochoterena (2000) simple indel coding via
   `workflow/scripts/encode_indels.py`. Each distinct *internal* contiguous gap
   (count >= 2) becomes a binary character; leading/trailing gaps are treated as
   MISSING DATA (`?`), not coded, because a gap at an alignment end is
   indistinguishable from the extraction window clipping that haplotype.
   Invariant columns are dropped so the indel matrix is variable-only. Emits
   `{peak_id}_nucleic.phy` + `{peak_id}_indel.phy` + a partition `.nex`
   (`GTR+I+G` on DNA, `GTR2+ASC` on the indel partition; `+ASC` corrects the
   ascertainment bias from recording only variable indel characters).

4. `rule iqtree` -- per-peak ML tree from the partitioned (DNA + indel) matrix.
   Cannot use `shadow: "minimal"` (conflicts with the absolute `--prefix`), so a
   post-run `find` deletes every aux file sharing the prefix except `.treefile`.
   Models are pinned by the `.nex` charpartition, so IQ-TREE estimates parameters
   but does NOT run model selection. Bootstrap (`-B 1000`) is conditional:
   enabled only when the nucleic PHYLIP has >=4 distinct sequences (IQ-TREE
   refuses bootstrap otherwise).

4b. `rule root_tree` (`workflow/scripts/root_tree.R`) -- roots the ML
   `.treefile` ONCE into `results/trees/{peak_id}.rooted.nwk`, which BOTH
   `phylo_mm` and `plot` read (neither re-roots), so the clade each branch
   represents is identical in the stats and the figure. Root = the **ancestral
   proxy**: the haplotype(s) closest to the MSA consensus (the modal allele ≈
   ancestral in young intraspecific variation), via `ape::root(outgroup=…)`;
   falls back to midpoint when that set isn't monophyletic, then to the first
   tip. This polarizes the tree ancestral→derived without midpoint's
   molecular-clock assumption and puts the reference allele at the base, which
   the sister-clade contrast (below) needs to read cleanly. Needs the alignment
   (to compute the consensus) as an input alongside the treefile.

5. `rule phylo_mm` -- per-peak phylogenetic mixed-model statistics via
   `workflow/scripts/phylo_mm.R`. Outputs `results/phylo_mm/{peak_id}.tsv` +
   `.pdf`. Three analyses:

   - **H²ₚ multi-prior sensitivity.** Intercept-only MCMCglmm fit 3 times with
     `alpha.V ∈ {1000, 1, 0.25}` on the phylogenetic random effect (PX prior;
     marginal on σ_phylo is half-Cauchy with scale √alpha.V). Trait is
     logit-transformed (with a coverage-derived offset for 0/1 saturation),
     residual likelihood uses `mev = 1/(coverage · p · (1-p))` (binomial SE
     on logit scale, the Hadfield 2010 meta-analysis trick), tree is
     ultrametricized via `chronos()` before `inverseA`. Posterior summary
     is median + 95% HPD. **Plot headline is the `alpha.V=1` fit** (Gelman 2006
     default half-Cauchy(0,1) on σ_phylo); when the `alpha.V=0.25`
     conservative-prior fit differs by ≥0.05 the plot shows
     `H²ₚ=X.XX (α=0.25→Y.YY) ⚠` to flag prior dependence. `alpha.V=1000` is
     retained in the TSV as a stress test but **not** in the headline:
     on weakly-identified peaks (small n_unique, big within-cluster spread)
     the wide half-Cauchy(0,~32) systematically inflates H²ₚ toward 1
     without adding identification information. R prior stays at
     `V=1, nu=0.002` (data-dominated); a `tightR` (V=1, nu=1) sensitivity
     is computed in the diagnostic stderr but not surfaced — it introduces
     bias on well-identified peaks where σ²_R is genuinely small.

   - **Per-split sister-clade contrast (headline).** Each internal node splits
     the tree into two daughter subtrees. For the SMALLER daughter (the split's
     representative, flagged `is_split`) a coverage-weighted WLS contrasts its
     actuation against its SISTER subtree (the other daughter), **not** against
     the rest of the tree. This is the local phylogenetic contrast — "does this
     lineage differ from the one it just split from?" — so a clade that looks
     significant only because the global mean lumps in a distant high/low clade
     is correctly NOT significant. BH FDR is over the split representatives (one
     test per node). A split is tested when both daughters have ≥2 covered
     leaves and the contrast spans ≥2 distinct sequence clusters (this skips
     arbitrary within-identical-cluster resolutions without dropping a real
     split where one identical leaf happens to sit on the sister side, which
     the LASSO's `clean_partition` would wrongly reject — e.g. an
     ancestral/derived root split). `Δ = mean_in − mean_out` (clade − sister)
     is the headline effect, `p_fdr` the headline p.

   - **LASSO joint branch selection (secondary).** Over `clean_partition`
     bipartitions (no identical-cluster straddle), QR-reduced to a linearly
     independent set (tree clades are collinear: parent = child₁ + child₂),
     `cv.glmnet` (alpha=1, standardize=FALSE, lambda.min) +
     `selectiveInference::fixedLassoInf()` give a joint post-selection β + CI.
     Retained in the TSV / plot hover as a secondary effect-size view; it NO
     LONGER drives the comparison or which branches are annotated, because the
     joint β over correlated nested clades is unstable (sign-flips vs the
     marginal Δ). OLS fallback (`inference_method=ols_fallback`) when the
     polyhedral check fails.

   The TSV header carries `h2_alpha{1000,1,025}_{mean,lower,upper}`,
   `n_branches_tested` (= split reps tested), `n_splits_sig` (BH<0.05),
   `n_branches_selected` (LASSO), `lasso_lambda`, `inference_method`.
   Per-branch rows carry `branch_id`, `descendants`, sister-contrast stats
   (`mean_in`, `mean_out`, `t_stat`, `p_raw`, `p_fdr`), `is_split`,
   `linked_sites`, plus secondary LASSO columns (`selected_by_lasso`,
   `lasso_coef`, `beta`, `se`, `ci_lower`, `ci_upper`, `p`).

5. `rule plot` -- produces `results/plots/{peak_id}.html` and `results/plots/{peak_id}.pdf`
   via `workflow/scripts/plot_tree_actuation_html.py`. Layout left → right:
   ancestral-proxy-rooted ladderized tree, leaf labels, **PCLAI ancestry stripe**,
   FIRE actuation heatmap (light-grey gaps with dark box outlines around
   contiguous `-` runs, N/unmapped stay dark grey), MSA with per-cell
   consensus coloring. The header surfaces H²ₚ (alpha.V=1 point estimate **with
   its 95% HPD**, alpha=0.25 flank when sensitivity ≥0.05 with `⚠`, and
   `⚠ unidentified` when the HPD spans >0.5 of [0,1] — i.e. the posterior is
   nearly flat and the point estimate is meaningless, distinct from
   prior-sensitivity), the significant-split count, and the standard peak stats.
   Each significant split (`is_split` rep with `p_fdr < 0.05`) is annotated on
   the tree with a single horizontal segment at the split's branch and a
   `Δ=… p_fdr=…` label (Δ = clade − sister, the local contrast; the joint LASSO
   β + CI live in hover only). One annotation per node, each in its own color,
   labels anchored to avoid the left edge near the root and staggered
   above/below to avoid overlap. Linked variants (when a split's bipartition is
   mirrored by an MSA SNP) are dashed-highlighted in the MSA in the split's
   color. Light dotted leader lines run from each tip to the sample-name
   column. The PCLAI column is shown only when `config["pclai_bed"]` points at a
   valid bgzipped + tabix-indexed BED; otherwise it collapses to zero width.

## Key implementation decisions

- **`iqtree` not `iqtree2`**: the conda-forge package installs as `iqtree` (version 3+).
- **Bootstrap gated on unique seq count**: `count_unique_seqs_phylip()` parses the nucleic
  PHYLIP. Counting `>` headers is wrong because IQ-TREE deduplicates internally. (This
  counts DNA-only uniqueness, so two tips identical in DNA but distinct in indel coding
  count as one — slightly conservative for the >=4 bootstrap gate.)
- **Indel-partition `+ASC`**: the binary indel matrix is variable-only by construction
  (every retained column has both a `0` taxon and a `1` taxon; `encode_indels` drops any
  invariant column), so ascertainment-bias correction is required or the indel-partition
  branch lengths inflate. Terminal gaps are coded `?` (missing), not as indel state — see
  rule `encode_indels` above and the "regulatory-element lengths" note below.
- **What H²<sub>p</sub> means here**: the tree is the regulatory element's *own* sequence
  genealogy and the trait is actuation of that same element, so H²<sub>p</sub> is "how much
  actuation variance tracks sequence-relatedness at this locus", NOT transmissible
  heritability. It is a within-locus descriptive statistic; the comparative-method
  machinery (independent lineages, a species tree) does not hold for one short intraspecific
  element. Treat it as descriptive, especially on prior-sensitive (`⚠`) or non-converged
  (`[ESS low]`) peaks.
- **Paralog mode + gene conversion (caveat)**: for SD paralog families (NOTCH2NL), a strictly
  bifurcating ML tree is the wrong model — interlocus gene conversion makes the promoters
  evolve reticulately, so haplotypes cluster by their most-recently-converted tract rather
  than by orthology. Read the paralog-mode tree (and its H²<sub>p</sub>) with that in mind.
- **Regulatory-element lengths**: peaks are short (~150-300 bp) with real per-haplotype
  length variation (the fixture spans 238-282 bp). Length differences at the *alignment ends*
  are partly an extraction-window artifact (the per-hap sequence is the assembly substring at
  the projected consensus-peak coords, not the true element boundary), which is why
  `encode_indels` masks terminal gaps as `?`. The principled fix is to extract with flanking
  sequence and track the RE sub-span in the MSA (planned, upstream change).
- **`min_seqs: 2`** (config default): allows peaks with only 2 distinct haplotype sequences.
  With `min_seqs: 3` (old default), some valid peaks were silently skipped.
- **Actuation = `fire_cov_OV_cov`** (column 9, 0-indexed): `fire_coverage / coverage`,
  a 0-1 fraction. Stored in `fastas/{peak_id}.tsv` alongside the FASTA.
- **`haplotype_distinct` is a peak-level filter**: a peak is included if any primary-sample
  row has the target classification. All qualifying sequences go into the tree regardless.
- **Streaming extraction**: `extract_fastas.py` processes one peak at a time (constant memory),
  relying on the BED being sorted by `consensus_peak_id` so rows for a peak are consecutive.
- **Columns resolved by header NAME, not fixed index** (`extract_fastas.resolve_columns`):
  the union-peaks BED carries a `#chrom ...` header, and the column order has changed
  across dataset builds (e.g. `fire_cov_OV_cov` moved, and `asm_chr/asm_start/asm_end`
  were collapsed into a single `asm_peak_id = '{asm_chr}_{start}_{end}'`). Resolving by
  name survives both. asm coords come from explicit columns if present, else from
  splitting `asm_peak_id`. Falls back to legacy fixed indices only for a headerless file.
- **HG002 maternal/paternal is FLIPPED in the peaks BED** (confirmed against the AGC:
  `HG002#1#` = PATERNAL, `HG002#2#` = MATERNAL, i.e. hap1=paternal/hap2=maternal). The BED
  has hap1→`chrN_MATERNAL`, hap2→`chrN_PATERNAL` — backwards. The *sequences* are correct
  (they match their labeled parental assembly exactly), only the `Haplotype` 1/2 label is
  swapped, so HG002 leaves carry the wrong hap number vs. the rest of the cohort. Not fixed
  in the asm-coord rebuild as of 2026-06; needs an upstream fix (swap HG002 Haplotype 1↔2).
  HG002 uses non-PanSN T2T contig names (`chr12_MATERNAL`); `add_flanks.agc_contig` maps
  these to the archive's `HG002#2#chrN` (maternal) / `HG002#1#chrN` (paternal), keyed on the
  asm_chr suffix (the true assembly), so HG002 now flanks correctly *despite* the label flip.
- **Variant anchoring is RE-scoped; the split test and H²ₚ are not.** When `flank > 0`,
  `phylo_mm.R` restricts the variable sites that anchor branches (`var_cols`) — hence the
  LASSO candidate set and the linked-site annotation — to the regulatory-element columns
  (`re_band`, reported as `re_first_col`/`re_last_col` in the TSV header); neutral flank
  variation only resolves the tree. The **headline sister-clade split test is NOT
  column-scoped** — it contrasts daughter-subtree actuation regardless of where the
  defining variants sit, so a split driven by a flank variant still gets a Δ + p_fdr (it
  just carries no linked RE site). H²ₚ uses the full flanked tree. The plot labels linked
  sites RE-relative (`RE13` = 13 bp into the element) via the same band.
- **`phylo_mm.R` uses `standardize=FALSE` in cv.glmnet**: indicator columns are already
  on the same 0/1 scale, and `selectiveInference::fixedLassoInf` requires the same
  parameterisation as the LASSO that produced the coefficients (otherwise the polyhedral
  check fails). With our 0/1 indicators, standardisation also implicitly down-weights
  small-clade branches — we want all branches penalised equally.
- **LASSO uses `lambda.min`, not `lambda.1se`**: small effective N (peaks with only 2
  unique phenotype levels) makes lambda.1se overly conservative; it routinely picks 0
  branches even where the per-branch FDR scan flags real splits.
- **H²ₚ identifiability depends on `n_unique` and within-cluster spread, not tree depth
  per se**. Peaks with many haps per unique sequence (78 haps / 15 unique) have σ²_R
  pinned by within-cluster replication; H²ₚ is data-driven and stable across alpha.V
  priors. Peaks with few unique sequences AND large within-cluster spread on the logit
  scale (e.g. chr15/chr20/chr22 in the test fixture: 2-3 unique, wvar(y) ≈ 5-7 vs
  mev ≈ 0.5) have σ²_phylo unidentified by the data; the prior dictates the result.
  The 3-alpha sensitivity sweep (`alpha.V ∈ {1000, 1, 0.25}`) surfaces this in the
  diagnostic stderr (`[phylo_mm] <peak> H²: alpha=...`). Headline H²ₚ uses alpha=1;
  alpha=0.25 conservative-prior flank with `⚠` is shown when they differ by ≥0.05.
  alpha=1000 is *not* used in the headline because its scale-32 half-Cauchy on σ_phylo
  systematically pushes σ²_phylo (and H²ₚ) toward huge values on weakly-identified peaks
  without adding information beyond the prior itself.
- **Headline test is the per-split sister-clade contrast, `p_fdr` (BH over splits).**
  Each internal node is tested ONCE: its smaller daughter (`is_split`) vs its sister
  daughter, coverage-weighted WLS, then BH FDR over the split reps. The comparison is
  the SISTER, not the rest of the tree — so a sub-clade that differs from the global
  mean only because the tree contains a distant high/low clade is correctly *not*
  significant against its own sibling. This collapses the headline from "every clade
  that beats the global mean" (which on a high-H²ₚ peak is 10+ correlated, sign-flipping
  branches) to the few real shifts, and makes Δ a coherent local effect (Δ ≈ the WLS
  coefficient, no joint-β instability). The plot annotation reads `Δ=… p_fdr=…`.
- **The sister contrast replaced two earlier designs** that both produced confusing
  output on high-H²ₚ peaks: (a) clade-vs-*rest-of-tree* (every derived sub-clade looked
  significant against the high-actuation reference background), and (b) clade-vs-LASSO-
  *container* (Δ conditioned on the LASSO selection hierarchy, but lambda.min over the
  collinear clade set over-selected and the joint β sign-flipped vs the marginal Δ).
  The sister (rest-of-parent) comparison is the local phylogenetic contrast and needs
  no LASSO hierarchy to be coherent.
- **LASSO + `selectiveInference::fixedLassoInf` are now SECONDARY** (TSV + plot hover
  only). They still run on the QR-reduced `clean_partition` set and supply a joint
  post-selection β + CI, but no longer choose what to draw or the comparison group.
  Kept because the joint multi-branch effect is occasionally a useful cross-check;
  demoted because (a) the polyhedral check fails on collinear tree clades (→ OLS
  fallback with invalid CIs), and (b) the joint β can disagree in sign with the local
  Δ, which is exactly the confusion the sister contrast removes.
- **Identical-haplotype collapse** (Hadfield 2010 / Halliwell BRV 2025): tips with
  EXACT same aligned sequence (gaps counted; only truly-identical strings collapse)
  share a single `animal` factor level in the MCMCglmm random effect. The collapsed
  tree (via `ape::keep.tip`) gives a non-singular `inverseA`; within-cluster
  observation variance lands in σ²_R rather than inflating σ²_phylo. H²ₚ is the
  honest between-sequence-cluster heritability. `n_unique_seqs` is reported alongside
  `n_tips` in the TSV header. Peaks with `n_unique_seqs < 2` (monomorphic
  alignment) are written with empty TSV + note and skipped from any modeling.
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
        phylo_mm.R                MCMCglmm H²ₚ + LASSO branch selection + selectiveInference
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
    fastas/    {peak_id}.fa  {peak_id}.tsv  peak_ids.txt
    aligned/   {peak_id}.fa  (temp, deleted after plot)
    trees/     {peak_id}.treefile
    phylo_mm/  {peak_id}.tsv  {peak_id}.pdf
    plots/     {peak_id}.html  {peak_id}.pdf
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

Resolved by header NAME (see `extract_fastas.resolve_columns`), not fixed index —
the column order differs between dataset builds.

| Column | Used for |
|--------|----------|
| `consensus_peak_id` | grouping key, output filename |
| `fire_cov_OV_cov` | actuation fraction for heatmap |
| `asm_peak_id` (or `asm_chr`/`asm_start`/`asm_end`) | per-sample asm contig + coords; flanks + PCLAI. New builds carry one `asm_peak_id = '{asm_chr}_{start}_{end}'`; old builds had the three split columns. |
| `seq` | DNA sequence (skip if NA) |
| `Individual_ID` | FASTA header / tree leaf prefix |
| `Haplotype` | FASTA header / tree leaf suffix; matches PCLAI's `#N#` directly (but see HG002 flip above) |
| `is_primary_sample` | filter when `primary_only: true` |
| `haplotype_distinct` | optional peak-level filter |

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
