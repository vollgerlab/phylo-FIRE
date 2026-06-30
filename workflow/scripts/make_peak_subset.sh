#!/usr/bin/env bash
# Carve a small peak-subset BED out of the full union-peaks BED, given a list of
# consensus_peak_ids (one per line). This is how the test fixtures (e.g.
# genetic-candidates, hap-selective-N) are minted, and how you scope any run to
# a chosen set of peaks. The consensus_peak_id column is resolved by HEADER NAME
# (the union BED's column order varies across builds), so this survives rebuilds.
#
# Usage:
#   make_peak_subset.sh <peak_ids.txt> <full_union.bed.gz> <out.bed.gz>
#
# Example (12 genetic candidates):
#   bash workflow/scripts/make_peak_subset.sh config/genetic_candidates_peak_ids.txt \
#        resources/union-peaks-asm.bed.gz resources/genetic-candidates/genetic-candidates.bed.gz
set -euo pipefail

IDS="${1:?usage: make_peak_subset.sh <peak_ids.txt> <full_union.bed.gz> <out.bed.gz>}"
FULL="${2:?missing full union BED}"
OUT="${3:?missing output bed.gz}"
mkdir -p "$(dirname "$OUT")"

bgzip -dc "$FULL" | awk -F'\t' -v ids="$IDS" '
  BEGIN { while ((getline line < ids) > 0) if (line != "") keep[line] = 1 }
  NR == 1 {
    for (i = 1; i <= NF; i++) if ($i == "consensus_peak_id") pc = i
    if (!pc) { print "ERROR: no consensus_peak_id column in header" > "/dev/stderr"; exit 1 }
    print; next
  }
  ($pc in keep) { print }
' | bgzip -@ 4 > "$OUT"

n_peaks=$(bgzip -dc "$OUT" | tail -n +2 | cut -f1 | wc -l)
echo "[make_peak_subset] wrote $OUT ($(wc -l < "$IDS") ids requested, $n_peaks rows)" >&2
