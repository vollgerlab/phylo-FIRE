# Developer notes

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

# Cluster (CHPC) production run — TWO stages (see "Two-stage split" below):
pixi run hsca-hapsel          # Stage 1: per-peak trees + phylo_mm (expensive)
pixi run hsca-hapsel-stage2   # Stage 2: re_variants + hg38 + plasmid + plots (minutes)
# full 180K: `pixi run -- snakemake --configfile config/cluster.yaml --rerun-triggers mtime -j 400`
# then `pixi run hsca-stage2`.
```

There is **no `clean` task**: snakemake handles its own state via
`--forcerun <rule>`, `--rerun-incomplete`, or `--forceall`. Don't manually
`rm -rf results/`.

## Two-stage split

The pipeline is split into two Snakemake entry points sharing one
`results/{run}/` layout (a single 180K-peak DAG took a day+ just to *build*):

- **Stage 1 — `workflow/Snakefile`** (expensive, run once): the per-peak
  sequence → alignment → tree → `phylo_mm` chain. Persists
  `results/{run}/aligned/` (the flanked MSA, now a real output — was `temp`),
  `trees/{peak}.rooted.nwk`, `phylo_mm/{peak}.tsv`, and one consolidated
  `tables/haplotypes.parquet` (per-hap RE seq + actuation + asm coords, replacing
  the 180K `fastas/*` which are now `temp`). `rule all` targets exactly these.
- **Stage 2 — `workflow/analysis.smk`** (cheap, rerunnable, ~6-node DAG that
  builds in seconds): `re_variants_all → (aggregate) → project_all (hg38/CHM13) →
  annotate_rsids → plasmid_all → plot_all`. Each is ONE single-pass script
  looping peaks internally (the `build_plasmid_table.py` pattern), reading only
  the Stage-1 dirs + `haplotypes.parquet` + the AGC — **never recomputes a
  tree**. Run: `snakemake -s workflow/analysis.smk --config run_name=<run>
  --rerun-triggers mtime -j <N>` (or the `*-stage2` pixi tasks). Add
  `--config dbsnp=<dbsnp.vcf.gz>` to also emit `variants.hg38.rsid.vcf.gz`.

**Why:** iterating on any downstream analysis (the hg38 projection fix, rs IDs,
plasmid, plots) reruns only Stage 2 — no trees, no MCMCglmm. Stage 2 is
pure Python/bcftools, so it also runs **locally** after downloading
`results/{run}/{aligned,trees,phylo_mm}` + `tables/haplotypes.parquet` and using
the local AGC (`resources/HPRC_r2_assemblies_0.6.1.agc`) — no genome fasta, no
`fastas/` dir. Sequences come from the parquet + AGC; the AGC lacks actuation
(the phenotype lives only in the union BED), which is why Stage 1 emits it once.

**AGC is 0-based.** `agc getctg contig:S-E` returns the standard 1-based base
`S+1`; `project_variants.py` passes 0-based starts directly (a historical `+1`
had made every reported hg38/CHM13 position 1 too low — fixed, verified
52764/52764 REF bases match GRCh38).

**Driver projection uses the Stage-1 alignment, flagged by `linked_sites`
column.** `project_variants.py` projects onto the pipeline's own
`aligned/{peak}.fa` (adding the CHM13/hg38 ref pulls with `mafft --add`, which
only inserts all-hap-gap columns so hap columns — and thus phylo_mm's — are
preserved; a `_col_bridge` maps them to the ref-added alignment for coords). A
driver is flagged at exactly its phylo_mm `linked_sites` column(s), the SAME
anchor the tree/MSA plot uses, so the bed can't disagree with the plot. This
replaced an earlier design that RE-only **re-aligned** and matched drivers by
**carrier-set**, which silently dropped (a) flank drivers (not in the RE-only
input) and (b) **homoplastic** drivers whose clade (`descendants`) isn't any
single column's carriers — the MHC/HLA Group-A peaks Annelise flagged (e.g.
`chr6_32416483` went from 0 projected rows to 17 driver SNVs). Carrier-set match
is now only a fallback (RE-band-scoped) when a driver's linked columns aren't
substitutions. Drivers defined only by an **indel** (a contiguous gap run, no
single-base column — e.g. `chr6_31400514`) can't be a bed/VCF SNV, so they are
written to `tables/drivers.unlifted.tsv` with a reason instead of a fabricated
position. `re_first_col`/`re_last_col`/`attribution_mode` come from the phylo_mm
header; `linked_sites` are 1-based alignment columns.

**Divergent driver sets are flagged, never filtered.** A driver set whose
`linked_sites` cover more than `DIVERGENT_FRAC` (5%, `build_re_variants.py`) of the
attributed band (RE band, or whole alignment in `region` mode) is a divergent
haplotype (paralog / mismapped copy / misalignment) rather than a clonable variant:
e.g. `chr12_131726083` has 111 linked SNVs in a 226 bp RE. The flag is
`driver_divergent` in `re_variants` (comma-sep, aligned to `driver_ids`),
`divergent` in `variants.hg38.tsv.gz` / `drivers.unlifted.tsv` / the variant BED,
and `DIVERGENT` in the VCF INFO. In `hsca_hapsel` these few sets carry ~1/4 of all
driver SNVs, so count drivers per set (not per SNV) in enrichment tests.

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

   In both cases the peak set comes from the static manifest
   `results/<run>/fastas/peak_ids.txt`, read at parse time.

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

   - **H²ₚ — RE-only first, expand to whole region, and gate on it.** The
     intercept-only MCMCglmm (`compute_h2`) is fit on TWO scopes: the **RE-only**
     tree/alignment (`rule re_trim`→`re_iqtree`→`re_root_tree`: the same mafft
     alignment sliced to the RE band, run through the identical encode→iqtree→
     root pipeline) and the **whole flanked** tree. Each scope fits the
     `alpha.V ∈ {1000, 1, 0.25}` (+ `tightR`) sensitivity (PX prior; logit trait
     with coverage offset; `mev = 1/(cov·p·(1-p))`; `chronos()` then `inverseA`,
     node labels stripped first). The clusters for σ²_R are by **that scope's**
     sequence identity, so high-diversity peaks with big within-cluster
     actuation spread come out unidentified (wide HPD). A fit is **heritable**
     when its alpha.V=1 HPD width ≤ `h2_hpd_width_max` (0.5, identified) AND its
     lower bound ≥ `h2_hpd_lower_min` (0.1, confidently > 0). **Attribution
     mode**: RE if the RE-only H² is heritable → attribute **RE-band** variants;
     else if the whole-region H² is heritable → mode `region`, attribute the
     **whole-alignment** (RE + flank) variants; else mode `none` → no drivers
     (e.g. chr22, both unidentified — identical RE sequences scatter 0.2–0.7 in
     actuation, so it isn't sequence-driven). The headline `h2_alpha1_*` is the
     winning scope's; `h2_re_*`, `h2_region_*`, `attribution_mode` carry both.
     Thresholds are config params (`h2_hpd_width_max` / `h2_hpd_lower_min`). The
     `.pdf` H² panel shows the winning scope. (At flank=0 the RE *is* the whole
     alignment, so `_re_tree`/`_re_alignment` fall back to the main tree.)

   - **Variant attribution — forward-selection joint regression (headline).**
     The question is "given this heritable trait, which variants contribute?",
     so the candidates are the VARIABLE SITES themselves (RE-scoped via
     `var_cols`), not tree splits. Each site is a 0/1 minor-allele indicator;
     sites carrying the SAME partition (perfectly correlated / complementary)
     merge into one block reporting all positions ("multiple variants on one
     branch"). Forward selection adds the candidate that explains the most
     ADDITIONAL weighted-R² of actuation (coverage-weighted WLS), stopping when
     the best gain falls below `R2_MIN` (0.05 = 5% of variance) or stops being
     significant. We select on **variance explained, not p-value**: at N≈78 a
     β≈0.08 variant adding <2% variance still clears p<0.05, so a significance
     gate over-selects weak passengers (chr19: pos 77 explains 85%, two others
     ~1% each yet "significant"). The reported `beta` is the JOINT (conditional)
     contribution — so a passenger that merely tracks a stronger variant adds ~0
     R² and is dropped, while a NESTED variant that switches actuation back off
     among a stronger site's carriers explains real variance and is kept with a
     negative β (verified: Y β=+0.50 "creates", W⊂Y β=−0.48 "removes"). Few
     candidates (2–10 blocks), so the joint fit is exact and stable — no LASSO /
     selectiveInference. `mean_in`/`mean_out` are the marginal carriers-vs-rest
     means; `beta`/`se`/`ci`/`p` the joint fit; `is_split` flags the selected
     (contributing) variants.

   The TSV header carries `h2_alpha{1000,1,025}_{mean,lower,upper}`,
   `n_variants_tested`, `n_variants_selected`, `n_branches_selected` (alias),
   `inference_method=forward_wls`. Per-variant rows carry `branch_id` (block
   index), `n_below` (n carriers), `descendants` (the carrier leaves, so the
   plot can highlight them), marginal `mean_in`/`mean_out`/`t_stat`/`p_raw`/
   `p_fdr`, `is_split` (= `selected_by_lasso`, both flag a contributing
   variant), `linked_sites` (the RE position(s)), and the joint `beta`/`se`/
   `ci_lower`/`ci_upper`/`p`. Epigenetic peaks (flanks resolve the tree but the
   RE is monomorphic) get a real H² with an empty variant table.

5. `rule plot` -- produces `results/plots/{category}/{peak_id}.html` and `.pdf`
   via `workflow/scripts/plot_tree_actuation_html.py`. `{category}` is one of
   `drivers` / `heritable_no_driver` / `random` (see `aggregate_results` below) —
   plots are NOT rendered for every peak, only the selected subset. Layout left → right:
   ancestral-proxy-rooted ladderized tree, leaf labels, **PCLAI ancestry stripe**,
   FIRE actuation heatmap (light-grey gaps with dark box outlines around
   contiguous `-` runs, N/unmapped stay dark grey), MSA with per-cell
   consensus coloring. The header surfaces H²ₚ (alpha.V=1 point estimate **with
   its 95% HPD**, alpha=0.25 flank when sensitivity ≥0.05 with `⚠`, and
   `⚠ unidentified` when the HPD spans >0.5 of [0,1] — i.e. the posterior is
   nearly flat and the point estimate is meaningless, distinct from
   prior-sensitivity), the contributing-variant count, and the standard peak
   stats. Each contributing variant (`is_split`) is annotated at its **MSA
   column** (top of the alignment panel) with a `β=… p=…` label + its
   RE-relative position(s) — NOT on a tree branch, because a variant's carriers
   are often not monophyletic (homoplasy / recurrent allele), so the column is
   the reliable anchor. β is the JOINT contribution (a nested "removal" variant
   shows negative β). The carriers' MSA column is dashed-highlighted and their
   leader lines coloured to match; when the carriers happen to form a clade its
   branch is also drawn in the same colour. Within each contiguous block of
   IDENTICAL aligned sequences the rows are ordered by % actuation (high at top)
   so the gradient inside a haplotype cluster is visible (leaf order is reordered
   but y stays monotonic so the heatmaps render correctly). Light dotted leader
   lines run from each tip to the sample-name column. The PCLAI column is shown
   only when `config["pclai_bed"]` points at a valid bgzipped + tabix-indexed
   BED; otherwise it collapses to zero width.

6. `rule re_variants` (`workflow/scripts/build_re_variants.py`) -- the
   **validation deliverable**: one row per *unique RE sequence* per peak, built
   from already-produced files (the mafft alignment + actuation TSV + the
   `phylo_mm` TSV), no re-modeling. Haplotypes identical over the attributed band
   collapse to one row. Grain columns: `re_seq` (literal clonable DNA), `re_len`,
   `gc_frac`, cov-weighted `act_mean`+spread, `act_values` (comma-sep raw per-hap
   actuations), `samples` (provenance). Driver columns (inline): `has_driver`,
   `driver_ids`/`driver_positions`/`driver_betas`/`driver_p` (comma-sep when a
   sequence carries >1), and `driver_delta` = sum of the joint β of the driver
   alleles it carries = the additive model's predicted actuation shift vs the
   non-carrier baseline (the number to validate). **Plasmid-design helpers**
   (within-peak Hamming): `nearest_noncarrier_seq_id`, `n_other_diffs` (diffs
   OUTSIDE the driver columns; `0` ⇒ a natural single-variable pair exists),
   `clean_pair`, `minimal_edit` (`RE80 c->g`, the exact mutagenesis). The
   attributed band = RE band normally, whole alignment in `region` mode (so the
   clonable unit always contains the drivers). RE-relative positions are
   alignment-column-based (`col - re_first_col + 1`, matching the plot).

7. `checkpoint aggregate_results`
   (`workflow/scripts/aggregate_re_variants.py`) -- concatenates every peak's
   `re_variants` TSV into `results/tables/re_variants.parquet` (+ `.tsv.gz`), the
   unified table, and writes three plot manifests that drive data-driven
   plotting: `plot_drivers.txt` (peaks with ≥1 selected driver), 
   `plot_heritable_no_driver.txt` (H² identified — `h2_re_ok` or `h2_region_ok` —
   but 0 drivers; **unlimited**), and `plot_random.txt` (`config["random_n"]`=100
   seeded-random **background** peaks — sampled from those that are NEITHER
   drivers NOR heritable-no-driver, i.e. the epigenetic / not-sequence-driven
   majority — as the negative-control contrast to the drivers dir; it does NOT
   overlap the other two). It is a **checkpoint** because
   `rule all`'s plot targets depend on which peaks land in each manifest.
   `rule all` = the parquet table + the categorized plots; requiring the table
   forces the full per-peak stats (tree → phylo_mm → re_variants) for the WHOLE
   peak set, so stats are computed for all peaks and only the selected subset is
   plotted.

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
results/{run_name}/           all pipeline outputs (gitignored)
    fastas/      {peak_id}.fa  {peak_id}.tsv  peak_ids.txt
    trees/       {peak_id}.treefile  {peak_id}.rooted.nwk
    phylo_mm/    {peak_id}.tsv          (no PDF; --pdf is ad-hoc only)
    re_variants/ {peak_id}.tsv          per-peak unique-RE-sequence rows
    tables/      re_variants.parquet  re_variants.tsv.gz  plot_*.txt  (aggregate)
    plots/       {category}/{peak_id}.html  .pdf   (category ∈ drivers/heritable_no_driver/random)
temp/{run_name}/              all temp() intermediates (gitignored; kept by --notemp)
    flanked/     {peak_id}.fa
    aligned/     {peak_id}.fa
    re/aligned/  {peak_id}.fa
```

**Run namespacing**: every output is under `results/{run_name}/` + `temp/{run_name}/`
(`run_name` config, default `default`; each pixi task sets its own). Separate
runs never clobber or orphan each other and can run concurrently (separate
manifests). All `temp()` outputs live under `temp/{run_name}/`, so the pixi test
tasks (which pass `--notemp`) leave every intermediate in one inspectable dir;
production runs without `--notemp` delete them.

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
