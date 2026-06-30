#!/usr/bin/env bash
# Write the consensus_peak_ids for ONE haplotype_distinct class out of the HSCA
# index (resources/hsca/hsca-index.tsv: col 1 = class, col 6 = peak id). Used to
# scope a cluster run to a single class.
#
# Usage: hsca_class_ids.sh "<class>" <out.txt>
#   classes: "Genetically haplotype-specific" | "Haplotype-selective"
#            | "Epigenetically haplotype-selective" | "Complex"
set -euo pipefail

CLASS="${1:?usage: hsca_class_ids.sh \"<class>\" <out.txt>}"
OUT="${2:?missing output path}"
IDX="resources/hsca/hsca-index.tsv"
[[ -f "$IDX" ]] || { echo "ERROR: $IDX missing — run 'pixi run setup-hsca' first" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")"
tail -n +2 "$IDX" | awk -F'\t' -v c="$CLASS" '$1==c {print $6}' > "$OUT"
n=$(wc -l < "$OUT" | tr -d ' ')
[[ "$n" -gt 0 ]] || { echo "ERROR: 0 peaks for class '$CLASS' (typo? check col 1 of $IDX)" >&2; exit 1; }
echo "[hsca_class_ids] wrote $n '$CLASS' peak ids -> $OUT" >&2
