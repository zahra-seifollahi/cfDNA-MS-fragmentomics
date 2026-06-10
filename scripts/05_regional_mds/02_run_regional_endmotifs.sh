#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Run FinaleToolkit interval-end-motifs for regional MDS
#
# Purpose:
#   Extract 4-mer end-motif frequencies for predefined marker
#   regions using FinaleToolkit interval-end-motifs.
#
# Input:
#   1. Fragment files directory containing *.frag.gz
#   2. BED file of 167 marker regions
#   3. hg19 reference genome in .2bit format
#   4. Output directory
#
# Usage:
#   bash scripts/05_regional_mds/02_run_regional_endmotifs.sh \
#     /path/to/Fragment_files \
#     /path/to/ms_panel_regions_unique.bed \
#     /path/to/hg19.2bit \
#     /path/to/output_interval_endmotifs
#
# Output:
#   *.interval_endmotifs.tsv
# ============================================================

FRAG_DIR="${1:-/path/to/fragment_files}"
BED_FILE="${2:-metadata/ms_panel_regions_unique.bed}"
REF_2BIT="${3:-/path/to/hg19.2bit}"
OUT_DIR="${4:-results/intermediate/regional_mds/interval_endmotifs}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MAPQ=30
KMER=4
MIN_LEN=50

mkdir -p "$OUT_DIR"

echo "Fragment directory: $FRAG_DIR"
echo "BED file: $BED_FILE"
echo "Reference genome: $REF_2BIT"
echo "Output directory: $OUT_DIR"
echo "Threads: $THREADS"
echo

if ! command -v finaletoolkit >/dev/null 2>&1; then
    echo "Error: finaletoolkit not found in PATH."
    exit 1
fi

if [ ! -d "$FRAG_DIR" ]; then
    echo "Error: fragment directory not found: $FRAG_DIR"
    exit 1
fi

if [ ! -f "$BED_FILE" ]; then
    echo "Error: BED file not found: $BED_FILE"
    exit 1
fi

if [ ! -f "$REF_2BIT" ]; then
    echo "Error: reference .2bit file not found: $REF_2BIT"
    exit 1
fi

n_files=$(find "$FRAG_DIR" -maxdepth 1 -name "*.frag.gz" | wc -l | tr -d ' ')

if [ "$n_files" -eq 0 ]; then
    echo "Error: no .frag.gz files found in $FRAG_DIR"
    exit 1
fi

echo "Number of fragment files found: $n_files"
echo

for frag in "$FRAG_DIR"/*.frag.gz; do
    [ -e "$frag" ] || continue

    sample=$(basename "$frag" .frag.gz)

    echo "Processing $sample ..."

    finaletoolkit interval-end-motifs \
        "$frag" \
        "$REF_2BIT" \
        "$BED_FILE" \
        -k "$KMER" \
        -min "$MIN_LEN" \
        -q "$MAPQ" \
        -w "$THREADS" \
        -o "$OUT_DIR/${sample}.interval_endmotifs.tsv" \
        -v

    echo "Finished $sample"
    echo
done

echo "Regional interval end-motif extraction completed."
echo "Outputs saved in:"
echo "$OUT_DIR"
