#!/usr/bin/env bash
# Build a compressed HSCA-only peak table from the full union-peaks BED, for
# quickly minting new test cases (we never run all of it — it's a reference).
#
# A peak is HSCA if any primary-sample row carries a haplotype-selective
# `haplotype_distinct` class (Genetically haplotype-specific / Epigenetically
# haplotype-selective / Haplotype-selective / Complex). We emit:
#   resources/hsca/hsca-peaks.bed.gz   all rows for every HSCA peak (bgzip)
#   resources/hsca/hsca-index.tsv      per-peak class + actuation spread, for
#                                      picking examples without scanning the BED
#
# Usage: bash workflow/scripts/build_hsca_table.sh [full_union_bed.gz]
set -euo pipefail

FULL="${1:-resources/union-peaks-asm.bed.gz}"
OUT_DIR="resources/hsca"
IDS="$(mktemp)"
trap 'rm -f "$IDS"' EXIT
mkdir -p "$OUT_DIR"

is_hsca='($38=="Genetically haplotype-specific" || $38=="Epigenetically haplotype-selective" || $38=="Haplotype-selective" || $38=="Complex") && $37=="TRUE"'

echo "[hsca] pass 1: indexing HSCA peaks (class + actuation spread) ..." >&2
bgzip -dc "$FULL" | awk -F'\t' "
NR==1{next}
$is_hsca {
  pk=\$5; a=\$13+0; cat[pk]=\$38; cnt[pk]++;
  if(!(pk in mn)||a<mn[pk]) mn[pk]=a;
  if(!(pk in mx)||a>mx[pk]) mx[pk]=a;
}
END{ for(pk in cnt) printf \"%s\t%d\t%.3f\t%.3f\t%.3f\t%s\n\", cat[pk], cnt[pk], mn[pk], mx[pk], mx[pk]-mn[pk], pk }
" | sort -t$'\t' -k5,5gr > "$OUT_DIR/hsca-index.body.tsv"

{ printf 'haplotype_distinct\tn_primary_haps\tmin_act\tmax_act\tact_spread\tconsensus_peak_id\n'; \
  cat "$OUT_DIR/hsca-index.body.tsv"; } > "$OUT_DIR/hsca-index.tsv"
rm -f "$OUT_DIR/hsca-index.body.tsv"
cut -f6 "$OUT_DIR/hsca-index.tsv" | tail -n +2 > "$IDS"
echo "[hsca] $(wc -l < "$IDS") HSCA peaks" >&2

echo "[hsca] pass 2: extracting all rows for those peaks ..." >&2
bgzip -dc "$FULL" | \
  awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($5 in keep){print}' "$IDS" - | \
  bgzip -@ 4 > "$OUT_DIR/hsca-peaks.bed.gz"

echo "[hsca] wrote $OUT_DIR/hsca-peaks.bed.gz + hsca-index.tsv" >&2
