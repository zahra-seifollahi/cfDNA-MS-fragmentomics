#!/usr/bin/env bash
#SBATCH --job-name=global_endmotifs
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=global_endmotifs_%j.out
#SBATCH --error=global_endmotifs_%j.err

# ============================================================
# Purpose:
#   Extract global 4-mer cfDNA end-motif frequencies from
#   fragment files using FinaleToolkit.
#
# Input:
#   Sorted and indexed fragment files (*.frag.gz)
#   Reference genome in 2bit format
#
# Output:
#   One end-motif frequency table per sample:
#   sample.endmotif.tsv
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

FRAG_DIR="${1:-/path/to/fragment_files}"
OUT_DIR="${2:-/path/to/output_endmotif_tables}"
REF_2BIT="${3:-/path/to/reference_genome.2bit}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
MAPQ=30
KMER=4

# -----------------------------
# Checks
# -----------------------------

if ! command -v finaletoolkit >/dev/null 2>&1; then
    echo "Error: finaletoolkit is not available in PATH." >&2
    exit 1
fi

if [ ! -d "$FRAG_DIR" ]; then
    echo "Error: fragment directory does not exist: $FRAG_DIR" >&2
    exit 1
fi

if [ ! -f "$REF_2BIT" ]; then
    echo "Error: reference 2bit file does not exist: $REF_2BIT" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "Fragment directory: $FRAG_DIR"
echo "Output directory: $OUT_DIR"
echo "Reference genome: $REF_2BIT"
echo "Threads: $THREADS"
echo "MAPQ: $MAPQ"
echo "K-mer size: $KMER"
echo

# -----------------------------
# Run FinaleToolkit
# -----------------------------

for FRAG in "$FRAG_DIR"/*.frag.gz; do
    [ -e "$FRAG" ] || {
        echo "No .frag.gz files found in: $FRAG_DIR"
        exit 0
    }

    SAMPLE=$(basename "$FRAG" .frag.gz)
    OUT_FILE="$OUT_DIR/${SAMPLE}.endmotif.tsv"

    echo "Processing: $SAMPLE"

    finaletoolkit end-motifs "$FRAG" "$REF_2BIT" \
        -q "$MAPQ" \
        -k "$KMER" \
        -o "$OUT_FILE" \
        -w "$THREADS" \
        -v

    echo "Done: $OUT_FILE"
    echo
done

echo "Global end-motif extraction completed."
