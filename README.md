# phylo-FIRE

Snakemake pipeline that asks, for each consensus FIRE peak, whether chromatin
actuation tracks the local sequence genealogy across HPRC haplotypes, and if so,
which variants drive it. Each peak gets a flanked MSA, an ML tree, a phylogenetic
mixed-model H²ₚ estimate and a forward-selection driver model. The deliverable is
a per-RE-sequence variant table with hg38/CHM13 coordinates and plasmid designs.

## Quick start

```bash
pixi install
pixi run test               # Stage 1 on 9 test peaks (results/test/)
pixi run replot             # Stage 2 on the test run, re-rendering plots
pixi run hap-selective-30   # both stages on 30 HSCA peaks
pixi run test -n            # dry run; tasks pass extra flags to snakemake
pixi run unit               # unit tests for the scripts
```

## Two stages

Both stages write to `results/{run_name}/` and share `config/config.yaml`.

**Stage 1, `workflow/Snakefile`** (rules in `workflow/rules/`). The expensive
per-peak chain: extract sequences from the union BED, add AGC flanks, align with
MAFFT, code indels, build an IQ-TREE tree, root it, then run `phylo_mm.R`
(H²ₚ gate, then joint WLS driver selection). Key outputs:

| Path | Contents |
|------|----------|
| `aligned/`, `trees/`, `phylo_mm/` | per-peak MSA, rooted tree, phylo_mm TSV |
| `tables/{alignments,trees,phylo_mm}.parquet` | the per-peak files combined |
| `tables/haplotypes.parquet` | per-haplotype sequence, actuation, assembly coords |
| `fastas/peak_ids.txt` | peaks that survived extraction |

**Stage 2, `workflow/analysis.smk`**. Single-pass jobs that read only the
Stage-1 tables, so they finish in minutes and never recompute a tree:

| Path | Contents |
|------|----------|
| `tables/re_variants.{parquet,tsv.gz}` | one row per unique RE sequence per peak, with drivers |
| `tables/variants.hg38.{tsv.gz,vcf.gz}` | every variable site projected to hg38 |
| `tables/variants.{hg38,chm13}.bed.gz`, `elements.hg38.bed.gz` | tabix-indexed BEDs |
| `tables/drivers.unlifted.tsv` | indel-only drivers that cannot be a single-base site |
| `plasmid_fiberseq/` | plasmid construct panel (`plasmid=true`) |
| `plots/{drivers,heritable_no_driver,random}/` | tree + MSA + actuation PDFs for a subset |

Optional: `dbsnp=<vcf.gz>` adds `variants.hg38.rsid.vcf.gz`, and
`significant_loci=<tsv>` (with `dbsnp`) adds `phewas/significant_driver_rsids.*`.

Stage 2 runs anywhere once you copy `results/<run>/tables/` (~40 MB) and
`fastas/peak_ids.txt`. The hg38/CHM13 projection also needs the AGC archive and
the union BED.

## Cluster runs (CHPC Redwood)

```bash
pixi run hsca-hapsel-ids && pixi run hsca-hapsel   # Stage 1, ~23K peaks
pixi run hsca-hapsel-stage2                        # Stage 2
```

The full 180K-peak run uses `config/cluster.yaml`; see `docs/data-setup.md`.

## More detail

- `config/config.yaml`: every config key, with comments.
- `docs/data-setup.md`: input data, the HSCA peak set and new test fixtures.
- `docs/DEVELOPMENT.md`: rule list, design decisions and data gotchas.
