# Data setup: union peaks → HSCA + test fixtures

How to get the input data, carve the HSCA production set, and mint test
fixtures. All paths are relative to the repo root.

## 1. The union-peaks BED (the input)

The pipeline ingests a consensus-FIRE-peaks **union BED** with one row per
peak × haplotype-sample (sequence, actuation `fire_cov_OV_cov`, asm coords,
`haplotype_distinct` class, `is_primary_sample`, …). Columns are resolved by
**header name** (`extract_fastas.resolve_columns`), so the column order may
differ across builds — don't hardcode indices.

Production copy (CHPC Redwood, `config/cluster.yaml`):

```
input_table: /scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/fiber-seq-pilot/union-peaks-asm.bed.gz   # ~5.6 GB, asm-coord-sorted
agc:         /scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/impg/HPRC_r2_assemblies_0.6.1.agc          # flank sequence (flank=250)
```

(A chm13-coord variant on hyak is used for tabix region pulls — see
`resources/hap-selective-5/`.)

## 2. Get it locally (for testing / building fixtures)

`config/config.yaml` defaults to `input_table: resources/union-peaks-asm.bed.gz`,
so copy the union BED (and the AGC, for flanks) into `resources/`:

```bash
rsync -avP <user>@redwood15.chpc.utah.edu:/scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/fiber-seq-pilot/union-peaks-asm.bed.gz resources/
rsync -avP <user>@redwood15.chpc.utah.edu:/scratch/ucgd/lustre-labs/vollger/projects/HPRCv2/impg/HPRC_r2_assemblies_0.6.1.agc        resources/
```

The union BED and AGC are large and gitignored; never commit them.

## 3. HSCA reference table + the production peak-id list

```bash
pixi run setup-hsca      # carves resources/hsca/{hsca-peaks.bed.gz, hsca-index.tsv} from the union BED
pixi run hsca-peak-ids   # writes config/hsca_peak_ids.txt (180,684 ids) from the index
```

- `hsca-peaks.bed.gz` — all rows for every HSCA peak (a peak is HSCA if any
  primary-sample row carries a haplotype-selective `haplotype_distinct` class).
- `hsca-index.tsv` — one row per peak: `haplotype_distinct`, `n_primary_haps`,
  `min_act`, `max_act`, `act_spread`, `consensus_peak_id`. This is the **menu**
  for picking examples without rescanning the BED.
- `config/hsca_peak_ids.txt` — the 180,684 HSCA peak ids; `config/cluster.yaml`
  pins the production run to this (gitignored; regenerate with the task).

## 4. Mint a test fixture (a small peak subset)

Pick peak ids — e.g. sort `hsca-index.tsv` by class / `act_spread` / hap count —
into a `config/<name>_peak_ids.txt`, then carve a fixture BED:

```bash
bash workflow/scripts/make_peak_subset.sh \
    config/genetic_candidates_peak_ids.txt \
    resources/union-peaks-asm.bed.gz \
    resources/genetic-candidates/genetic-candidates.bed.gz
```

`make_peak_subset.sh` resolves `consensus_peak_id` by header name and emits a
bgzipped subset (small fixtures are committed; the full union BED is not). Add a
pixi task that points `peak_ids` + `input_table` + `run_name` at the fixture
(see existing `test`, `genetic-candidates`, `hap-selective-30`, … in
`pixi.toml`).

How the current fixtures were screened (for reference):
- **genetic-candidates** — 12 Genetically-haplotype-specific peaks chosen by
  within-cluster actuation spread: 6 low-spread + high-unique-seq (identified
  RE-H² → drivers, incl. the nested chr14_100272296) + 6 high-spread (RE-H²
  fails → mode `none`).
- **hap-selective-30** — hap-selective-20 + 10 high-unique-seq-count peaks
  (likely multi-variant drivers).

Selection is eyeballing `hsca-index.tsv` (`haplotype_distinct`, `act_spread`)
plus an occasional quick unique-sequence count.

## 5. Run HSCA on the cluster

```bash
# one-time, on the cluster (after the union BED + AGC are in place)
pixi install
pixi run setup-hsca && pixi run hsca-peak-ids

# verify the partition/account first
sinfo

# stage on a slice (override run_name + a smaller peak list) before the full run
head -5000 config/hsca_peak_ids.txt > config/hsca_stage.txt
pixi run -- snakemake --configfile config/cluster.yaml \
    --config run_name=hsca_stage peak_ids=config/hsca_stage.txt \
    --profile profiles/slurm -j 200 -n          # dry-run first; watch DAG build + head RAM

# full run — outputs under results/hsca/, temp under temp/hsca/
pixi run -- snakemake --configfile config/cluster.yaml --profile profiles/slurm -j 200
```

The SLURM profile batches the per-peak chain via `group: "peak_chain"` +
`group-components: peak_chain=100` so ~180K peaks don't become ~2M sbatch jobs;
calibrate the component count and the per-job allocation on the staging slice.
Do **not** pass `--notemp` in production (it would keep `temp/hsca/`). The
deliverable is `results/hsca/tables/re_variants.parquet`; plots render only for
the drivers / heritable-no-driver / random subsets.
