#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Select hyper- and hypomethylated regions for Healthy samples
#
# Purpose:
#   For each QC-passed Healthy sample, select:
#     - 52 hypermethylated regions: beta > 0.89, highest beta values
#     - 52 hypomethylated regions: beta < 0.25, lowest beta values
#
# Input:
#   BED directory containing files like:
#     01_hyper.bed
#     01_hypo.bed
#     02_hyper.bed
#     02_hypo.bed
#     ...
#
# Expected BED columns:
#   chr start end beta
#
# Usage:
#   bash scripts/06_hyper_hypo_regions/01_select_hyper_hypo_regions.sh \
#     /path/to/input_bed_directory \
#     /path/to/output_selected_bed_directory \
#     results/tables/qc/final_sample_metadata.csv
#
# Output:
#   selected BED files:
#     01_hyper_beta_gt_0.89_top52.bed
#     01_hypo_beta_lt_0.25_top52.bed
#
# Notes:
#   This script uses only QC-passed Healthy samples.
# ============================================================

BED_DIR="${1:-/path/to/input_hyper_hypo_beds}"
OUT_DIR="${2:-results/intermediate/hyper_hypo_regions/selected_beta_regions}"
QC_METADATA="${3:-results/tables/qc/final_sample_metadata.csv}"

HYPER_CUTOFF=0.89
HYPO_CUTOFF=0.25
TOP_N=52

mkdir -p "$OUT_DIR"

echo "Input BED directory: $BED_DIR"
echo "Output directory: $OUT_DIR"
echo "QC metadata: $QC_METADATA"
echo "Hyper cutoff: beta > $HYPER_CUTOFF"
echo "Hypo cutoff: beta < $HYPO_CUTOFF"
echo "Top N per class: $TOP_N"
echo

if [ ! -d "$BED_DIR" ]; then
    echo "ERROR: BED directory not found: $BED_DIR"
    exit 1
fi

if [ ! -f "$QC_METADATA" ]; then
    echo "ERROR: QC metadata file not found: $QC_METADATA"
    exit 1
fi

# ------------------------------------------------------------
# Extract QC-passed Healthy sample numbers from metadata.
#
# Expected metadata columns include:
#   sample, group, include_analysis
#
# sample should look like:
#   Cap01, Cap02, ...
# ------------------------------------------------------------

healthy_samples=$(awk -F',' '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "sample") sample_col = i
        if ($i == "group") group_col = i
        if ($i == "include_analysis") include_col = i
    }

    if (!sample_col || !group_col || !include_col) {
        print "ERROR: metadata must contain sample, group, include_analysis columns" > "/dev/stderr"
        exit 1
    }

    next
}

$group_col == "Healthy" && $include_col == "yes" {
    sample = $sample_col
    gsub("Cap", "", sample)
    print sample
}
' "$QC_METADATA")

if [ -z "$healthy_samples" ]; then
    echo "ERROR: no QC-passed Healthy samples found in metadata."
    exit 1
fi

echo "QC-passed Healthy sample numbers:"
echo "$healthy_samples"
echo

n_processed=0
n_skipped=0

for i in $healthy_samples; do

    hyper_bed="$BED_DIR/${i}_hyper.bed"
    hypo_bed="$BED_DIR/${i}_hypo.bed"

    hyper_out="$OUT_DIR/${i}_hyper_beta_gt_0.89_top52.bed"
    hypo_out="$OUT_DIR/${i}_hypo_beta_lt_0.25_top52.bed"

    echo "Processing sample $i ..."

    if [ ! -f "$hyper_bed" ]; then
        echo "WARNING: missing hyper BED: $hyper_bed"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    if [ ! -f "$hypo_bed" ]; then
        echo "WARNING: missing hypo BED: $hypo_bed"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    awk -v cutoff="$HYPER_CUTOFF" 'BEGIN{OFS="\t"} $4 > cutoff {print $0}' "$hyper_bed" \
        | sort -k4,4nr \
        | head -n "$TOP_N" \
        > "$hyper_out"

    awk -v cutoff="$HYPO_CUTOFF" 'BEGIN{OFS="\t"} $4 < cutoff {print $0}' "$hypo_bed" \
        | sort -k4,4n \
        | head -n "$TOP_N" \
        > "$hypo_out"

    n_hyper=$(wc -l < "$hyper_out" | tr -d ' ')
    n_hypo=$(wc -l < "$hypo_out" | tr -d ' ')

    echo "  hyper selected: $n_hyper"
    echo "  hypo selected:  $n_hypo"

    if [ "$n_hyper" -ne "$TOP_N" ]; then
        echo "WARNING: sample $i has fewer than $TOP_N hyper regions after filtering."
    fi

    if [ "$n_hypo" -ne "$TOP_N" ]; then
        echo "WARNING: sample $i has fewer than $TOP_N hypo regions after filtering."
    fi

    n_processed=$((n_processed + 1))
    echo

done

echo "============================================================"
echo "Hyper/hypo region selection completed."
echo "Processed samples: $n_processed"
echo "Skipped samples: $n_skipped"
echo "Output directory:"
echo "$OUT_DIR"
echo "============================================================"
