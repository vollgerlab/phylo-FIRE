# Driver-variant schematic: figure spec

Five panels (a to e), left to right or 2+3 grid, ~180 mm wide (double column).
Toy data are **real** values from `chr8_70445067_70445271` (test run), subsampled to 14 tips.
Full numbers: `toy_data.json` (same directory). Do not invent values outside that file.

## Corrections to the plan

1. **Haplotype count.** The pipeline uses 39 individuals x 2 = **78 haplotypes** (primary samples;
   the exemplar and the production `haplotypes.parquet` both have 78; HG002 is included). The plan's
   37 x 2 = 74 does not match the code or data. Panel a below uses 78; confirm with the authors
   and swap both numbers together if the paper cohort really is 37.
2. **Two trees, not one.** The RE-only H² uses its own tree: the MSA is sliced to the RE columns and
   rerun through indel coding, IQ-TREE and rooting (`re_trim`, `re_iqtree`, `re_root_tree`).
   Panel c must show an RE-only tree and a flanked-region tree.
3. **Gate outcome sets the candidate sites.** RE-only H² passes: test RE sites only. RE fails but
   region passes: test all sites (RE + flanks). Both fail: no drivers. The plan omits this.
4. **Candidates are site blocks.** Sites with identical carrier sets are merged into one candidate
   (exemplar: RE74 + RE148). Gap is a valid allele; columns >50% gap/N are skipped. Indicator =
   second most common allele.
5. **Regression details.** Forward selection is weighted least squares on the raw 0 to 1 actuation
   (weights = fibers per haplotype), using haplotypes with >=4 fibers. H² instead uses logit
   actuation with a coverage-based measurement-error term on an ultrametric tree, identical
   sequences pooled (prior alpha.V = 1 is the headline). Keep these out of the labels; they are context.
6. **Root fallback.** Ancestral-proxy rooting falls back to midpoint when the closest-to-consensus
   set is not monophyletic. Panel b label stays "rooted on ancestral proxy".
7. **Stopping rule.** A step is accepted only if the best new site adds **>=5% R² AND P < 0.05**
   (its P in the model at that step). The plan's wording is right; the arrow on step 3 should cite both.
8. **Panel e statistic.** Label with the joint-model β and P (as in methods), not the marginal
   carrier-vs-rest difference (RE112 marginal Δ = -40% vs joint β = -0.33). β is in actuation
   units (fraction of fibers), so β = -0.33 means 33 percentage points less FIRE.
9. IQ-TREE models: GTR+I+G on DNA plus GTR2+ASC on the binary indel partition (Simmons-Ochoterena
   simple coding); 1000 ultrafast bootstraps when >=4 unique sequences. MAFFT runs with
   `--auto --maxiterate 1000 --adjustdirection`. Mention only "DNA + indel partitions".

## Palette (from `plot_tree_actuation_html.py` / `plot_tree_compact.py`)

| Use | Hex |
|---|---|
| Actuation bins 0-10% ... 90-100% (10 bins, `ACT_COLORS`) | `#ffffff #fff5f0 #fee0d2 #fcbba1 #fc9272 #fb6a4a #ef3b2c #cb181d #a50f15 #67000d` |
| Actuation missing | `#d0d0d0` |
| Driver 1 (RE112, closes; first by FDR) | `#d62728` |
| Driver 2 (RE63, opens) | `#1f77b4` |
| Linked, not selected (RE74+RE148, RE159) | `#888888` (suggested; not in code) |
| Rest / non-carriers: fill (alpha 0.4), point edge, mean line | `#bbbbbb`, `#888888`, `#666666` dashed |
| Tree branches / axis lines | `#444444` |
| Text / secondary text | `#222222` / `#444444` |
| DNA bases A / T / G / C (variable columns) | `#009600 / #FF4444 / #FFB300 / #4444FF` |
| MSA invariant / gap | `#EEEEEE` / `#555555` |
| Gate pass / fail (suggested; not in code) | `#2ca02c` / `#7f7f7f` |

Font: Arial (Helvetica fallback), editable text. Body 6-7 pt, panel letters 9 pt bold.
Caveat: driver red `#d62728` sits near the top of the actuation scale. Keep it for consistency
with the data figures, but never place a red driver mark directly on a red actuation cell.

## Panel a: Haplotype-resolved input

- Left: one individual icon split into two haplotype bars, `HG02257` -> `hap 1`, `hap 2`.
- Each bar: grey flank, RE block, grey flank. Labels above: `250 bp`, `RE`, `250 bp`.
- Under each RE: 4 to 6 short fiber lines; FIRE elements as filled red ticks over the RE.
  hap 1: 1 of 15 fibers FIRE (`7%`). hap 2: 7 of 14 fibers FIRE (`50%`). Draw ~1/6 and ~3/6.
- Actuation chip right of each bar, filled by bin color (`#fff5f0` for 7%, `#fb6a4a` for 50%).
- Bottom caption: `78 haplotypes (39 individuals x 2)`.
- Definition text: `FIRE actuation = fraction of fibers with a FIRE element`.

## Panel b: Align and build tree

- Stacked sequence bars (5 to 6 rows) -> arrow `MAFFT` -> small MSA block with flank shading and
  the RE band bracketed (`RE`). Show a few colored variable columns; one gap run.
- Below the MSA, a 2-row binary strip labeled `indels coded 0/1`.
- Arrow `IQ-TREE (DNA + indel partitions)` -> the 14-tip tree from `toy_data.json`
  (`tree_newick`), drawn rectangular, root at left.
- Mark the root group with a small open circle and label `rooted on ancestral proxy
  (closest to consensus)`.
- Tip colored chips = actuation (ACT_COLORS). No driver colors yet in this panel.

## Panel c: Heritability gate

- Two mini trees side by side: `RE only` and `RE ± 250 bp`, both with actuation chips.
- Model line: `Phylogenetic mixed model (MCMCglmm)` and `H² = phylogenetic share of variance`.
- Three H² posterior curves on x = 0 to 1, 95% HPD shaded under each:

| Curve | Label | Shape hint | Verdict |
|---|---|---|---|
| Pass (exemplar) | `H² = 1.00 [0.99, 1.00]` | spike at 1, e.g. Beta(80, 0.6) | `heritable` (pass color) |
| Too wide (chr22) | `H² = 0.49 [0.00, 0.93]` | broad hump, e.g. Beta(1.4, 1.4) | `width > 0.5` (fail) |
| Near zero (chr7) | `H² = 0.05 [0.00, 0.45]` | mass at 0, e.g. Beta(0.7, 9) | `lower < 0.1` (fail) |

- Gate box text: `Heritable if 95% HPD width ≤ 0.5 and lower bound ≥ 0.1`.
- Flow: `RE only` -> pass -> `test RE sites`; fail -> `RE ± 250 bp` -> pass -> `test all sites`;
  fail -> `no drivers`. The exemplar passes at the RE step (region H² 1.00 [0.96, 1.00]).

## Panel d: Forward selection

- Top: candidate table (one row per RE variable site, from `re_variable_sites`):
  `RE63 A>G (28)`, `RE74+RE148 (8)`, `RE112 T>C (10)`, `RE159 (7)`, and `RE170 (1)` struck through
  with note `<2 haplotypes per allele`. Bracket RE74 and RE148: `same carriers, merged`.
- Main: R² staircase. x = step (0 to 3), y = weighted R² (0 to 0.7, ticks 0, 0.2, 0.4, 0.6).

| Step | Label | R² | Bar color |
|---|---|---|---|
| 0 | `none` | 0.00 | none |
| 1 | `+ RE112` | 0.445 (ΔR² 44.5%) | `#d62728` |
| 2 | `+ RE63` | 0.604 (ΔR² 15.9%) | `#1f77b4` |
| 3 | `+ RE74+RE148` | 0.604 (ΔR² 0.05%, P = 0.77) | `#888888`, marked `✗ stop` |

- Dashed reference showing the 5% increment needed at step 3.
- Rule text: `Add the site that explains the most extra variance while ΔR² ≥ 5% and P < 0.05`.
- Subtitle: `coverage-weighted regression of actuation`.

## Panel e: Driver variants

- Left: the panel b tree repeated, now with RE112 stem and carrier tips in `#d62728`, RE63 stem
  and tips in `#1f77b4`. Mark RE112's two backgrounds: `bg1` (RE112 alone: HG01891_1, HG02257_1)
  and `bg2` (with RE74, RE148). Tag both: `closed on 2 backgrounds`.
- Right: two half-violin or strip rows (carriers below, rest above, as in `plot_tree_compact.py`),
  x = % FIRE 0 to 100. Real distributions: RE112 carriers 9/10 at 0%, one at 7%; rest mean 41%.
  RE63 carriers mean 51%; rest mean 27%. Dashed mean lines.
- Row labels (driver color, bold): `RE112 T>C` and `RE63 A>G`.
- Stat labels: `β = -0.33, P = 2e-10` and `β = +0.17, P = 5e-7`.
- Footer: `β and P from joint model of all drivers`. Use a true minus sign (U+2212) for -0.33.

## Global label rules

- Use `H²`, `ΔR²`, `β`, `P` (italic P optional); `±` not `+/-`; `≥`/`≤` glyphs.
- No em dashes anywhere. No GRCh38 coordinates or rsIDs (projection is out of scope).
- Haplotype names optional (≤5 pt) in b and e; drop them if crowded.
