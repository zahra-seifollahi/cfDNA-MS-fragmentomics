#!/usr/bin/env bash
#SBATCH --job-name=make_fragments
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=make_fragments_%j.out
#SBATCH --error=make_fragments_%j.err

# ============================================================
# Purpose:
#   Convert paired-end deduplicated BAM files into compressed
#   cfDNA fragment files for downstream fragmentomic analyses.
#
# Input:
#   Deduplicated BAM files (*.bam)
#
# Output:
#   Sorted and indexed fragment files:
#   sample.frag.gz
#   sample.frag.gz.tbi
#
# Notes:
#   This workflow was used for bisulfite-converted cfDNA data.
#   Only properly paired reads are retained.
# ============================================================

set -euo pipefail

# -----------------------------
# User settings
# -----------------------------

IN_DIR="${1:-/path/to/deduplicated_bam_files}"
OUT_DIR="${2:-/path/to/output_fragment_files}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# -----------------------------
# Checks
# -----------------------------

if ! command -v samtools >/dev/null 2>&1; then
    echo "Error: samtools is not available in PATH." >&2
    exit 1
fi

if ! command -v bedtools >/dev/null 2>&1; then
    echo "Error: bedtools is not available in PATH." >&2
    exit 1
fi

if ! command -v bgzip >/dev/null 2>&1; then
    echo "Error: bgzip is not available in PATH." >&2
    exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
    echo "Error: tabix is not available in PATH." >&2
    exit 1
fi

if [ ! -d "$IN_DIR" ]; then
    echo "Error: input directory does not exist: $IN_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "Input BAM directory: $IN_DIR"
echo "Output fragment directory: $OUT_DIR"
echo "Threads: $THREADS"
echo

# -----------------------------
# Convert BAM files
# -----------------------------

for BAM in "$IN_DIR"/*.bam; do
    [ -e "$BAM" ] || {
        echo "No BAM files found in: $IN_DIR"
        exit 0
    }

    BASE=$(basename "$BAM" .bam)
    OUT_FRAG="$OUT_DIR/${BASE}.frag.gz"

    echo "Processing: $BASE"

    samtools collate -O -u -@ "$THREADS" "$BAM" \
        | samtools view -@ "$THREADS" -bf 0x2 -F 0x904 - \
        | bedtools bamtobed -bedpe -mate1 -i stdin \
        | awk 'BEGIN{OFS="\t"}
               $1==$4 && $2>=0 && $5>=0 {
                 start = ($2 < $5 ? $2 : $5)
                 stop  = ($3 > $6 ? $3 : $6)
                 mapq  = $8
                 strand = $9
                 print $1, start, stop, mapq, strand
               }' \
        | sort -k1,1 -k2,2n \
        | bgzip -@ "$THREADS" > "$OUT_FRAG"

    tabix -p bed "$OUT_FRAG"

    echo "Done: $OUT_FRAG"
    echo
done

echo "All BAM files processed."
