#!/usr/bin/env bash

# ============================================================
# Run FinaleToolkit fragment-length bins
#
# Purpose:
#   Generate *.frag_length.tsv files from tabix-indexed
#   .frag.gz fragment files using FinaleToolkit.
#
# Input:
#   data/example/fragments/*.frag.gz
#
# Output:
#   data/example/finale_length_results/*.frag_length.tsv
#
# Usage:
#   bash scripts/03_fragment_length/00_run_finaletoolkit_frag_length_bins.sh
#
# Or:
#   bash scripts/03_fragment_length/00_run_finaletoolkit_frag_length_bins.sh \
#   data/example/fragments \
#   data/example/finale_length_results
# ============================================================

set -euo pipefail

FRAG_DIR="${1:-data/example/fragments}"
OUT_DIR="${2:-data/example/finale_length_results}"

MIN_LENGTH=1
MAX_LENGTH=500
BIN_SIZE=1

if ! command -v finaletoolkit >/dev/null 2>&1; then
    echo "Error: finaletoolkit is not available in PATH." >&2
    echo "Activate your environment first, for example:" >&2
    echo "source .venv_finale_test/bin/activate" >&2
    exit 1
fi

if [ ! -d "$FRAG_DIR" ]; then
    echo "Error: fragment directory does not exist: $FRAG_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "Fragment directory: $FRAG_DIR"
echo "Output directory: $OUT_DIR"
echo "Min length: $MIN_LENGTH"
echo "Max length: $MAX_LENGTH"
echo "Bin size: $BIN_SIZE"
echo

for FRAG in "$FRAG_DIR"/*.frag.gz; do
    [ -e "$FRAG" ] || {
        echo "No .frag.gz files found in: $FRAG_DIR"
        exit 0
    }

    SAMPLE=$(basename "$FRAG" .frag.gz)
    OUT_FILE="$OUT_DIR/${SAMPLE}.frag_length.tsv"

    echo "Processing: $SAMPLE"

    finaletoolkit frag-length-bins "$FRAG" \
        -min "$MIN_LENGTH" \
        -max "$MAX_LENGTH" \
        --bin-size "$BIN_SIZE" \
        -o "$OUT_FILE" \
        -v

    echo "Done: $OUT_FILE"
    echo
done

echo "FinaleToolkit fragment-length bin extraction completed."
