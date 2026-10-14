#!/usr/bin/env bash
set -euo pipefail

benchmark=${1:-./build/bench_exact_search}
prefix=${2:-tests/tmp_exact_search/test}
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

printf '>read1\nACG\n>read2\nTTG\n>unsupported\nACN\n' > "$tmpdir/reads.fa"
"$benchmark" "$prefix" "$tmpdir/reads.fa" 3

printf '@read1\nACG\n+\n!!!\n@read2\nGCA\n+\n###\n' > "$tmpdir/reads.fq"
"$benchmark" "$prefix" "$tmpdir/reads.fq" 3

: > "$tmpdir/empty.fa"
"$benchmark" "$prefix" "$tmpdir/empty.fa" 1

printf '@broken\nACG\n+\n' > "$tmpdir/broken.fq"
if "$benchmark" "$prefix" "$tmpdir/broken.fq" 1; then
    echo 'Malformed FASTQ unexpectedly succeeded' >&2
    exit 1
fi

if "$benchmark" "$prefix" "$tmpdir/reads.fa" 0; then
    echo 'Invalid repetition count unexpectedly succeeded' >&2
    exit 1
fi
