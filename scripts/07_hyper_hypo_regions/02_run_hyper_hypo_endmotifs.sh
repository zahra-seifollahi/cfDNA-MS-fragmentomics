#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Run FinaleToolkit interval-end-motifs for selected
# hypermethylated and hypomethylated regions
#
# Purpose:
#   For each QC-passed Healthy sample, run end-motif extraction
#   separately on:
#     - selected hypermethylated regions
#     - selected hypomethylated regions
#
# Input:
#   1. Fragment files directory
#   2. Selected BED directory from:
#        01_select_hyper_hypo_regions.sh
#   3. hg19 reference genome in .2bit format
#   4. Output directory
#   5. QC metadata
#
# Expected selected BED files:
#   01_hyper_beta_gt_0.89_top52.bed
#   01_hypo_beta_lt_0.25_top52.bed
#
# Expected fragment files:
#   Cap01.dedup.frag.gz
#   Cap02.dedup.frag.gz
#   ...
#
# Usage:
#   bash scripts/06_hyper_hypo_regions/02_run_hyper_hypo_endmotifs.sh \
#     /path/to/Fragment_files \
#     /path/to/selected_beta_regions \
#     /path/to/hg19.2bit \
#     /path/to/output_dir \
#     results/tables/qc/final_sample_metadata.csv
#
# Output:
#   interval_end_motifs/
#     Cap01_hyper_endmotif.tsv
#     Cap01_hypo_endmotif.tsv
# ============================================================

FRAG_DIR="${1:-/path/to/fragment_files}"
BED_DIR="${2:-results/intermediate/hyper_hypo_regions/selected_beta_regions}"
REF_2BIT="${3:-/path/to/hg19.2bit}"
OUT_DIR="${4:-results/intermediate/hyper_hypo_regions}"
QC_METADATA="${5:-results/tables/qc/final_sample_metadata.csv}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MAPQ=30
KMER=4
MIN_MOTIF_LEN=50

MOTIF_OUT_DIR="$OUT_DIR/interval_end_motifs"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$MOTIF_OUT_DIR"
mkdir -p "$LOG_DIR"

echo "============================================================"
echo "Hyper/hypo interval end-motif extraction"
echo "============================================================"
echo "Fragment directory:       $FRAG_DIR"
echo "Selected BED directory:   $BED_DIR"
echo "Reference genome:         $REF_2BIT"
echo "Output directory:         $OUT_DIR"
echo "Motif output directory:   $MOTIF_OUT_DIR"
echo "QC metadata:              $QC_METADATA"
echo "MAPQ:                     $MAPQ"
echo "K-mer:                    $KMER"
echo "Minimum motif length:     $MIN_MOTIF_LEN"
echo "Threads:                  $THREADS"
echo "============================================================"
echo

if ! command -v finaletoolkit >/dev/null 2>&1; then
    echo "ERROR: finaletoolkit not found. Activate your FinaleToolkit environment first."
    exit 1
fi

if [ ! -d "$FRAG_DIR" ]; then
    echo "ERROR: fragment directory not found: $FRAG_DIR"
    exit 1
fi

if [ ! -d "$BED_DIR" ]; then
    echo "ERROR: selected BED directory not found: $BED_DIR"
    exit 1
fi

if [ ! -f "$REF_2BIT" ]; then
    echo "ERROR: reference .2bit file not found: $REF_2BIT"
    exit 1
fi

if [ ! -f "$QC_METADATA" ]; then
    echo "ERROR: QC metadata file not found: $QC_METADATA"
    exit 1
fi

echo "FinaleToolkit:"
which finaletoolkit
finaletoolkit --version
echo

# ------------------------------------------------------------
# Extract QC-passed Healthy sample numbers from metadata.
# Expected metadata columns:
#   sample, group, include_analysis
# sample values:
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

    sample="Cap${i}"

    frag="$FRAG_DIR/${sample}.dedup.frag.gz"
    hyper_bed="$BED_DIR/${i}_hyper_beta_gt_0.89_top52.bed"
    hypo_bed="$BED_DIR/${i}_hypo_beta_lt_0.25_top52.bed"

    hyper_out="$MOTIF_OUT_DIR/${sample}_hyper_endmotif.tsv"
    hypo_out="$MOTIF_OUT_DIR/${sample}_hypo_endmotif.tsv"

    hyper_log="$LOG_DIR/${sample}_hyper_endmotif.log"
    hypo_log="$LOG_DIR/${sample}_hypo_endmotif.log"

    echo "------------------------------------------------------------"
    echo "Processing $sample"
    echo "Fragment file: $frag"
    echo "Hyper BED:     $hyper_bed"
    echo "Hypo BED:      $hypo_bed"
    echo "------------------------------------------------------------"

    if [ ! -f "$frag" ]; then
        echo "WARNING: missing fragment file: $frag"
        n_skipped=$((n_skipped + 1))
        continue
    fi

    if [ ! -f "${frag}.tbi" ]; then
        echo "WARNING: missing fragment index: ${frag}.tbi"
        n_skipped=$((n_skipped + 1))
        continue
    fi

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

    echo "Running interval end motifs: $sample hyper"

    finaletoolkit interval-end-motifs \
        "$frag" \
        "$REF_2BIT" \
        "$hyper_bed" \
        -q "$MAPQ" \
        -k "$KMER" \
        -min "$MIN_MOTIF_LEN" \
        -w "$THREADS" \
        -o "$hyper_out" \
        -v \
        > "$hyper_log" 2>&1

    echo "Finished end motifs: $sample hyper"

    echo "Running interval end motifs: $sample hypo"

    finaletoolkit interval-end-motifs \
        "$frag" \
        "$REF_2BIT" \
        "$hypo_bed" \
        -q "$MAPQ" \
        -k "$KMER" \
        -min "$MIN_MOTIF_LEN" \
        -w "$THREADS" \
        -o "$hypo_out" \
        -v \
        > "$hypo_log" 2>&1

    echo "Finished end motifs: $sample hypo"

    n_processed=$((n_processed + 1))
    echo

done

echo "============================================================"
echo "Hyper/hypo end-motif extraction completed."
echo "Processed samples: $n_processed"
echo "Skipped samples: $n_skipped"
echo "Output files:"
echo "$MOTIF_OUT_DIR"
echo "Logs:"
echo "$LOG_DIR"
echo "============================================================"
