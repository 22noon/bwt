#!/usr/bin/env bash
set -Eeuo pipefail

# Build both forward and reverse BWA BWT indexes for FM-index search experiments.
#
# Usage:
#   ./build_bwa_forward_bwt.sh reference.fa [output_prefix]
#
# Example:
#   ./build_bwa_forward_bwt.sh reference.fa reference
#
# Optional:
#   BWA=/path/to/bwa ./build_bwa_forward_bwt.sh reference.fa reference

usage() {
    echo "Usage: $0 <reference.fa[.gz]> [output_prefix]" >&2
    exit 2
}

[[ $# -ge 1 && $# -le 2 ]] || usage

FASTA=$1
RAW_PREFIX=${2:-"${FASTA%.*%.fa}"}
# Strip trailing .forward or .reverse if present to get clean base prefix
BASE_PREFIX="${RAW_PREFIX%.forward}"
BASE_PREFIX="${BASE_PREFIX%.reverse}"

FORWARD_PREFIX="${BASE_PREFIX}.forward"
REVERSE_PREFIX="${BASE_PREFIX}.reverse"

BWA=${BWA:-bwa}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REVCOMP="${SCRIPT_DIR}/revcomp.py"

if [[ ! -r "$FASTA" ]]; then
    echo "ERROR: Cannot read FASTA: $FASTA" >&2
    exit 1
fi

if ! command -v "$BWA" >/dev/null 2>&1; then
    echo "ERROR: BWA executable not found: $BWA" >&2
    exit 1
fi

if [[ ! -x "$REVCOMP" ]]; then
    echo "ERROR: revcomp.py not found or not executable at $REVCOMP" >&2
    exit 1
fi

# Avoid accidentally overwriting existing indexes.
for pfx in "$FORWARD_PREFIX" "$REVERSE_PREFIX"; do
    for ext in pac ann amb bwt sa; do
        if [[ -e "${pfx}.${ext}" ]]; then
            echo "ERROR: Output already exists: ${pfx}.${ext}" >&2
            echo "Choose a different output prefix or move existing files." >&2
            exit 1
        fi
    done
done

echo "Reference     : $FASTA"
echo "Forward Prefix: $FORWARD_PREFIX"
echo "Reverse Prefix: $REVERSE_PREFIX"
echo "BWA           : $("$BWA" 2>&1 | head -n 2 | tr '\n' ' ')"
echo

build_index() {
    local fa=$1
    local pfx=$2
    local label=$3

    echo "=== Building ${label} index (${pfx}) =="

    echo "[1/4] Packing reference..."
    "$BWA" fa2pac -f "$fa" "$pfx"

    [[ -s "${pfx}.pac" ]] || {
        echo "ERROR: Expected PAC file was not created for ${label}." >&2
        exit 1
    }

    echo "[2/4] Constructing BWT..."
    "$BWA" pac2bwt "${pfx}.pac" "${pfx}.bwt"

    [[ -s "${pfx}.bwt" ]] || {
        echo "ERROR: BWT file was not created for ${label}." >&2
        exit 1
    }

    echo "[3/4] Adding occurrence checkpoints..."
    "$BWA" bwtupdate "${pfx}.bwt"

    echo "[4/4] Adding SA locations..."
    "$BWA" bwt2sa -i 32  "${pfx}.bwt"  "${pfx}.sa"
    echo
}

# 1. Build Forward Index
build_index "$FASTA" "$FORWARD_PREFIX" "forward"

# 2. Generate Reverse Complement FASTA
REV_FASTA="${BASE_PREFIX}.rev_tmp.fa"
echo "Generating reverse-complemented reference FASTA..."
python3 "$REVCOMP" -i "$FASTA" -o "$REV_FASTA"

[[ -s "$REV_FASTA" ]] || {
    echo "ERROR: Failed to generate reverse-complemented FASTA." >&2
    exit 1
}

# Ensure cleanup of temporary file on exit
trap 'rm -f "$REV_FASTA"' EXIT

# 3. Build Reverse Index
build_index "$REV_FASTA" "$REVERSE_PREFIX" "reverse"

echo "Forward and reverse BWT index construction complete."
echo "Generated Files:"
for pfx in "$FORWARD_PREFIX" "$REVERSE_PREFIX"; do
    for ext in pac ann amb bwt sa; do
        [[ -f "${pfx}.${ext}" ]] &&
            printf '  %s\n' "${pfx}.${ext}"
    done
done
