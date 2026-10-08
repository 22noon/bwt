#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TMP="${ROOT}/tests/tmp_exact_search"
PREFIX="${TMP}/test"

rm -rf "${TMP}"
mkdir -p "${TMP}"

cat > "${TMP}/reference.fa" <<'EOF'
>test
ACGTTGCAACG
EOF

if [[ ! -x "${ROOT}/bwa" ]]; then
    echo "ERROR: expected BWA binary at ${ROOT}/bwa" >&2
    exit 1
fi

echo "Building test BWA/FMD index..."

"${ROOT}/bwa" index \
    -p "${PREFIX}" \
    "${TMP}/reference.fa"

echo "Running exact-search test..."

"${ROOT}/tests/test_exact_search" \
    "${PREFIX}.bwt" \
    "${PREFIX}.sa"
