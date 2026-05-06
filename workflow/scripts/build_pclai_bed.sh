#!/usr/bin/env bash
# Bulk-download HPRC R2 PCLAI per-haplotype BEDs in BOTH coordinate spaces
# (CHM13 + per-sample assembly), concatenate everything into one bgzip+tabix
# file. Run once; output is reused thereafter.
#
# Why both spaces? Most phylo-FIRE peaks are CHM13-native, so a single tabix
# query at the peak's CHM13 coords returns one row per sample-haplotype. But
# some peaks/loci have no CHM13 representation (HG002-only ``chr*_MATERNAL_*``
# consensus peaks; paralogs in segmental duplications that don't land on
# chm13). For those the input table provides per-sample asm-coord coordinates,
# and the asm-coord BEDs (PanSN-keyed: ``{sample}#{hap}#{contig}``) carry
# the PCLAI window for that haplotype's own assembly. Combining both into one
# tabix file lets the plot script try CHM13 first, fall back to asm-coord per
# leaf when CHM13 misses — without juggling two index files.
#
# Source index:
#   https://github.com/human-pangenomics/hprc_intermediate_assembly/tree/main/data_tables/annotation/pclai
#
# BED format (BED9+1, identical for both coord spaces):
#   chrom = CHM13 ref ("chr1" ... "chrM") OR PanSN ("{sample}#{hap}#{contig}")
#   col 4 = "{sample}/h{hap}/{chrom}_w{window}_(pc1,pc2)"
#   col 5 = score (confidence, higher = better)
#   col 9 = itemRgb (HPRC's pre-baked PCLAI palette)
#   col 10 = centroid PC coordinates "(pc1,pc2)" — one of 4 buckets
#
# Output:
#   resources/pclai/hprc_r2.bed.gz (+ .tbi)
set -euo pipefail

INDEX_BASE="https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/annotation/pclai"
S3_PREFIX="s3://human-pangenomics/"
HTTPS_PREFIX="https://human-pangenomics.s3.amazonaws.com/"

OUT_DIR="resources/pclai"
OUT_BED="${OUT_DIR}/hprc_r2.bed.gz"

PARALLEL_DOWNLOADS=${PARALLEL_DOWNLOADS:-16}

mkdir -p "${OUT_DIR}"

# fetch_space <coord_label> <raw_subdir>
#   - downloads the per-haplotype BEDs for one coordinate space,
#   - fails loudly if any download is missing,
#   - leaves the raw BEDs in raw_subdir for later concat.
fetch_space() {
    local coord="$1"      # "chm13" or "asm"
    local subdir="$2"
    local index_url="${INDEX_BASE}/pclai_v1.1_${coord}_coord_local_hprc_r2.index.csv"
    local index_local="${OUT_DIR}/pclai_v1.1_${coord}_coord.index.csv"
    local url_list="${OUT_DIR}/.${coord}_urls.txt"
    local raw_dir="${OUT_DIR}/${subdir}"

    mkdir -p "${raw_dir}"

    echo "[${coord}] fetching index..."
    curl -fsSL "${index_url}" -o "${index_local}"
    local n_haps
    n_haps=$(($(wc -l < "${index_local}") - 1))
    echo "       ${n_haps} sample-haplotypes listed."

    # GitHub-hosted CSVs are CRLF — strip the trailing CR before passing to curl.
    awk -F, -v s3="${S3_PREFIX}" -v https="${HTTPS_PREFIX}" '
        NR > 1 { sub(/\r$/, "", $4); sub(s3, https, $4); print $4 }
    ' "${index_local}" > "${url_list}"

    local n_urls
    n_urls=$(wc -l < "${url_list}" | tr -d ' ')
    [ "${n_urls}" -gt 0 ] || { echo "ERROR: empty URL list for ${coord}" >&2; exit 1; }

    echo "[${coord}] downloading ${n_urls} BEDs (parallel=${PARALLEL_DOWNLOADS})..."
    # No -I {} here: without it xargs appends the input line as the trailing
    # positional arg, which lets `$1` inside `bash -c` be the URL.
    RAW_DIR="${raw_dir}" xargs -n 1 -P "${PARALLEL_DOWNLOADS}" bash -c '
        url="$1"
        out="${RAW_DIR}/$(basename "$url")"
        if [ -s "$out" ]; then exit 0; fi
        curl -fsSL --retry 3 --retry-delay 2 "$url" -o "$out.part" \
            && mv "$out.part" "$out" \
            || { echo "FAILED: $url" >&2; exit 1; }
    ' _ < "${url_list}"

    local n_have
    n_have=$(find "${raw_dir}" -maxdepth 1 -name "*.pclai_v1.1.${coord}_coord.bed" | wc -l | tr -d ' ')
    if [ "${n_have}" -ne "${n_haps}" ]; then
        echo "ERROR [${coord}]: have ${n_have} BEDs but expected ${n_haps}" >&2
        exit 1
    fi
}

fetch_space "chm13" "raw_chm13"
fetch_space "asm"   "raw"

echo "[merge] concatenating both coord spaces, sorting, bgzipping..."
# CHM13 chroms ("chr1"..."chrM") and PanSN chroms ("HG00097#1#...") never
# collide, so a single LC_ALL=C sort on chrom,start gives a valid tabix layout.
# Tabix queries against either chrom space hit only rows in that space.
{
    find "${OUT_DIR}/raw_chm13" -maxdepth 1 -name '*.pclai_v1.1.chm13_coord.bed' -print0 \
        | xargs -0 cat
    find "${OUT_DIR}/raw"       -maxdepth 1 -name '*.pclai_v1.1.asm_coord.bed'   -print0 \
        | xargs -0 cat
} | LC_ALL=C sort -k1,1 -k2,2n -S 2G \
  | bgzip -@ 8 > "${OUT_BED}"

echo "[merge] tabix indexing..."
tabix -f -p bed "${OUT_BED}"

size=$(du -h "${OUT_BED}" | cut -f1)
nlines=$(bgzip -dc "${OUT_BED}" | wc -l | tr -d ' ')
echo "Done: ${OUT_BED} (${size}, ${nlines} rows)"
