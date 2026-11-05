# Developer notes

## Overview

phylo-FIRE has two Snakemake entry points that share `config/config.yaml` and the
`results/{run_name}/` layout. Stage 1 (`workflow/Snakefile`) is the expensive
per-peak chain and runs once. Stage 2 (`workflow/analysis.smk`) reads only the
Stage-1 summary tables and reruns in minutes. See `README.md` for outputs and
`docs/data-setup.md` for input data.

## Commands

| Task | Purpose |
|------|---------|
| `test`, `hap-selective-30`, ... | smoke tests; see "Test fixtures" |
| `replot` | Stage 2 on the `test` run, forcing `plot_all` |
| `unit` | pytest over `workflow/scripts/tests/` (AGC-backed tests skip without the archive) |
| `fmt` | `snakefmt workflow/` |
| `setup-pclai` | one-time: build `resources/pclai/hprc_r2.bed.gz` |
| `setup-hsca`, `hsca-peak-ids` | one-time: HSCA table and the 180K peak list |
| `hsca-hapsel-ids`, `hsca-hapsel` | Stage 1 on the Haplotype-selective class (~23K peaks) |
| `hsca-hapsel-tables` | rebuild the Stage-1 summary tables outside the DAG |
| `hsca-hapsel-stage2`, `hsca-stage2` | Stage 2 for `hsca_hapsel` and the full `hsca` run |

Every task accepts extra snakemake flags, e.g. `pixi run test -n` for a dry run.
There is no `clean` task. Use `--forcerun <rule>`, `-R`, `--rerun-incomplete`
or `--forceall` instead of deleting `results/`.

Full 180K run: `pixi run -- snakemake --configfile config/cluster.yaml
--rerun-triggers mtime -j 400`, then `pixi run hsca-stage2`.

## Stage 1 rules (`workflow/rules/`)

- `extract_fastas` (extract.smk): streams the union BED into temp
  `fastas/{peak}.{fa,tsv}`, `fastas/peak_ids.txt` and `tables/haplotypes.parquet`.
  With `preextracted_dir` it links ready-made `.fa`/`.tsv` pairs (e.g. PARE) instead.
- `add_flanks` (align.smk): adds `flank` bp of assembly context per haplotype from
  the AGC. Only runs when `flank > 0`.
- `mafft` (align.smk): `aligned/{peak}.fa`, persisted.
- `encode_indels`, `iqtree`, `root_tree` (trees.smk): indel coding, IQ-TREE ML
  tree, ancestral-proxy rooting to `trees/{peak}.rooted.nwk`.
- `re_trim`, `re_encode_indels`, `re_iqtree`, `re_root_tree` (trees.smk): the same
  chain on the MSA sliced to the RE, for the RE-only H²ₚ. Flanked runs only.
- `phylo_mm` (phylo_mm.smk): `phylo_mm.R` writes `phylo_mm/{peak}.tsv` (H²ₚ in
  the header, one row per candidate variant).
- `stage1_tables` (phylo_mm.smk): combines every peak's MSA, tree and phylo_mm
  TSV into `tables/{alignments,trees,phylo_mm}.parquet`.

## Stage 2 rules (`workflow/analysis.smk`)

Each rule is one job that loops over all peaks.

- `re_variants_all`: `tables/re_variants.{parquet,tsv.gz}` plus the plot manifests
  `plot_{drivers,heritable_no_driver,random}.txt`.
- `peak_coords`: CHM13 and hg38 coordinates per peak, from the union BED.
- `project_all`: every variable site to CHM13, hg38 and the carrier assemblies:
  `variants.hg38.{tsv.gz,vcf.gz}`, `variants.{hg38,chm13}.bed.gz`,
  `elements.hg38.bed.gz`, `drivers.unlifted.tsv`. Needs the AGC; skipped for
  `preextracted_dir` runs.
- `annotate_rsids`: `variants.hg38.rsid.vcf.gz`. Opt-in via `dbsnp`.
- `significant_driver_rsids`: `phewas/significant_driver_rsids.{tsv,txt}`.
  Opt-in via `significant_loci` plus `dbsnp`.
- `plasmid_all`: `plasmid_fiberseq/plasmid_constructs.{tsv,fa}`. Opt-in via `plasmid`.
- `plot_all`: `plots/{category}/{peak}.pdf` (and `.html` with `plot_html`) for
  every driver peak, every heritable-no-driver peak and `random_n` background peaks.

## Design decisions

**Two stages joined by summary tables.** One 180K-peak DAG took over a day just
to build. Stage 1 now ends in four parquet tables (`haplotypes`, `alignments`,
`trees`, `phylo_mm`), and Stage 2 reads only those plus `fastas/peak_ids.txt`, the
union BED (peak coordinates) and the AGC. Iterating on projection, rs IDs,
plasmids or plots never touches a tree. Stage 2 also runs locally after copying
`results/<run>/tables/` (~40 MB) and `fastas/peak_ids.txt`:

```bash
pixi run -- snakemake -s workflow/analysis.smk --config run_name=<run> \
    input_table=<union.bed.gz> peak_ids=<peaks.txt> plasmid=true \
    --rerun-triggers mtime -j 8
```

**Scripts are params, and runs use `--rerun-triggers mtime`.** Rules pass script
paths as `params`, never as `input`, so editing a script leaves outputs current.
The cluster and Stage-2 tasks also pass `--rerun-triggers mtime`, so a params or
rule-code edit cannot rerun 180K trees either. After a change that should alter
outputs, force it with `-R <rule>` or `--forcerun <rule>`.

**Static peak set, no checkpoint.** `PEAK_IDS` is resolved at parse time from
`preextracted_dir`, then `fastas/peak_ids.txt`, then config `peak_ids`.
Re-evaluating a checkpoint hangs at 100K+ peaks. A fresh run therefore needs
`peak_ids`; `haplotype_distinct` alone resolves zero peaks.

**Flanks resolve the tree; attribution stays RE-scoped.** With `flank > 0` each
haplotype gets AGC context, which adds variable sites and reduces the prior
dependence of H²ₚ. Haplotypes whose contig is not in the AGC fall back to RE-only,
so flanking never drops a haplotype. Variant candidates are the RE columns
(`re_first_col`/`re_last_col` in the phylo_mm header) unless the region scope wins.

**Indel coding with `+ASC`.** `encode_indels.py` codes each internal gap as a
binary character (Simmons-Ochoterena). Terminal gaps become `?` because they are
indistinguishable from extraction-window clipping. The indel partition is
variable-only, so it uses `GTR2+ASC`; DNA uses `GTR+I+G`. Models are fixed by the
`.nex`, so IQ-TREE skips model selection. Bootstrap (`-B 1000`) runs only with
4 or more unique sequences, because IQ-TREE refuses fewer.

**Ancestral-proxy rooting, done once.** `root_tree.R` roots on the haplotype(s)
closest to the MSA consensus, falling back to midpoint and then the first tip.
`phylo_mm.R` and the plot read the same rooted tree, so clades match.

**H²ₚ gate.** `phylo_mm.R` fits an intercept-only MCMCglmm on the RE-only and the
flanked tree, each under `alpha.V` = 1000, 1 and 0.25. A scope is heritable when
the `alpha.V=1` 95% HPD width is at most `h2_hpd_width_max` (0.5) and its lower
bound is at least `h2_hpd_lower_min` (0.1). The RE scope wins first, then
`region`, else `attribution_mode=none` and no drivers. `alpha.V=1000` is not the
headline because its prior inflates H²ₚ on weakly identified peaks. Haplotypes
with identical aligned sequences share one random-effect level, so
within-cluster spread lands in the residual. Peaks with few unique sequences and
large within-cluster spread come out unidentified (wide HPD).

**Forward-selection joint WLS, not LASSO.** Candidates are variable sites; sites
with the same carrier partition merge into one block. Selection adds the block
with the largest coverage-weighted R² gain and stops when the gain is below 5%
(`R2_MIN`) or p ≥ 0.05. It selects on variance explained because at N≈78 weak
passengers still reach p < 0.05. Betas are joint, so a nested variant that
reverses a stronger one is kept with the opposite sign. `is_split` marks drivers.
The columns `selected_by_lasso` and `lasso_coef` keep their legacy names for
compatibility; selection is forward WLS (`inference_method=forward_wls`).

**Driver projection uses the Stage-1 alignment.** `project_variants.py` adds the
CHM13/hg38 sequences to `aligned/{peak}.fa` with `mafft --add`, which leaves hap
columns intact, and flags a driver at its phylo_mm `linked_sites` columns. The
BED therefore matches the plot, including flank and homoplastic drivers.
Drivers defined only by an indel cannot be a single-base site and go to
`drivers.unlifted.tsv` with a reason.

**Divergent driver sets are flagged, not filtered.** A set whose `linked_sites`
cover more than 5% (`DIVERGENT_FRAC`, `build_re_variants.py`) of the attributed
band is a paralog, mismapped copy or misalignment rather than a clonable variant.
The flag is `driver_divergent` in `re_variants`, `divergent` in the hg38 TSV,
`drivers.unlifted.tsv` and the BED, and `DIVERGENT` in the VCF INFO. These few
sets carry about a quarter of driver SNVs in `hsca_hapsel`, so count drivers per
set in enrichment tests.

**What H²ₚ means.** The tree is the element's own genealogy and the trait is that
element's actuation, so H²ₚ measures how much actuation variance tracks local
sequence relatedness. It is descriptive, not transmissible heritability. For SD
paralog families (NOTCH2NL), gene conversion makes evolution reticulate, so read
paralog-mode trees with care.

## Data gotchas

- AGC is 0-based: `agc getctg contig:S-E` returns standard base `S+1`, so
  callers pass 0-based starts.
- Columns are resolved by header name (`extract_fastas.resolve_columns`).
  Column order changes between union-BED builds, and asm coords come either as
  `asm_chr`/`asm_start`/`asm_end` or as one `asm_peak_id`.
- HG002 haplotype labels are swapped in the union BED: hap1 is labelled
  MATERNAL, but the AGC has `HG002#1#` = paternal. Sequences are correct.
  `add_flanks.agc_contig` keys on the contig's MATERNAL/PATERNAL suffix, so flanks
  are right; the hap number on HG002 leaves is not. Needs an upstream fix.
- MAFFT `--adjustdirection` prefixes reversed sequences with `_R_`. Code that
  matches MSA names to samples must strip it.
- PCLAI (`pclai_bed`, built by `setup-pclai`): the merged BED holds both CHM13
  and assembly coordinates. The plot runs one tabix query per peak at the CHM13
  coords parsed from `peak_id` and matches rows to leaves by the
  `{sample}/h{hap}` prefix in column 4. HG002-only peaks (`chr8_MATERNAL_...`) get
  no stripe.

## Output layout

```
results/{run_name}/
    fastas/      peak_ids.txt (+ temp {peak}.fa/.tsv)
    aligned/     {peak}.fa
    trees/       {peak}.rooted.nwk
    phylo_mm/    {peak}.tsv
    tables/      Stage-1 parquet tables + every Stage-2 table
    plots/       {drivers,heritable_no_driver,random}/{peak}.pdf
    plasmid_fiberseq/  phewas/      (opt-in)
temp/{run_name}/ flanked/  coded/  trees/*.treefile  re/
```

Smoke tests pass `--notemp`, so temp files stay inspectable; production runs
delete them.

## Test fixtures

| Task | Fixture | What it exercises |
|------|---------|-------------------|
| `test` | `resources/test/test-union-peaks.bed.gz`, 9 peaks | mixed classes; `chr8_70445067_70445271` is the two-driver figure example |
| `hap-selective-30` | `resources/hap-selective-30/`, 30 peaks | both stages; 20 across classes, 10 with many unique RE seqs |
| `genetic-candidates` | `resources/genetic-candidates/`, 12 peaks | 6 low-spread (drivers) vs 6 high-spread (mode none) |
| `hap-selective-5` | `resources/hap-selective-5/`, 5 peaks | chm13-coord tabix pull |
| `hap-selective-5-samples` | same 5 peaks | only hap-selective individuals, `highlight_file` |
| `test-nested` | `resources/test/preextracted_nested/` | synthetic opposite-sign nested driver |
| `test-paralog` | `resources/test/preextracted_pg_31/`, 24 seqs | `preextracted_dir`, NOTCH2NL paralog group |
| `test-paralog-all-haps` | `resources/test/preextracted_pg_31_all_haps/`, 319 seqs | same group at full-cohort scale |

Only `hap-selective-30` runs Stage 2; the others run Stage 1 only.
`chr14_100272296` in genetic-candidates has two same-sign co-occurring drivers,
not a nested pair; `test-nested` is the only nested case, built by
`workflow/scripts/tools/make_nested_fixture.py`.

The paralog fixtures come from PARE (https://github.com/vollgerlab/PARE), whose
`paralog_msas/` output matches the `extract_fastas` format, so any PARE run can
be passed as `preextracted_dir`. `pg_31` is the chr1 NOTCH2NL promoter family
(NOTCH2, NOTCH2NLA/B/C and a NOTCH2NLR-like copy annotated as NBPF26), with headers
`{peak}__{gene}__{sample}_{hap}`. PARE group ids change between runs.

New fixtures: see `docs/data-setup.md`.

## SLURM notes

`profiles/slurm/config.yaml` targets CHPC Redwood. It batches the per-peak chain
(`group: "peak_chain"`, `group-components: peak_chain=100`) so 180K peaks do not
become ~2M jobs. Check `slurm_partition` and `slurm_account` with `sinfo` before
the first cluster run.
