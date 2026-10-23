#!/usr/bin/env bash
set -euo pipefail

fmsearch=$1
make_sa=$2
bwa=$3
work=$4
test_search=$5
run_loop_fmsearch=${6:-}
rm -rf "$work"
mkdir -p "$work"
prefix="$work/ref.forward"

cat > "$work/ref.fa" <<'EOF'
>toy
ACGTACGTAAAACCCCGGGGTTTTACGT
EOF
cat > "$work/reads.fq" <<'EOF'
@repeat
ACGT
+
abcd
@boundary
ACGT
+
efgh
@unique
TTTT
+
ijkl
@orientation
GTACGTA
+
abcdefg
@unmapped
AAAAAC
+
mnopqr
EOF

"$bwa" fa2pac -f "$work/ref.fa" "$prefix" >/dev/null
"$bwa" pac2bwt "$prefix.pac" "$prefix.bwt" >/dev/null
"$bwa" bwtupdate "$prefix.bwt" >/dev/null
"$make_sa" "$prefix.bwt" "$prefix.sa"
"$test_search" "$prefix.bwt" "$work/ref.fa" "$prefix.sa"

# Also exercise the same reversed-reference orientation used by the rank
# backend benchmark. The search test compares BWA and BatMIS intervals at
# every step, including the primary-row singleton probes.
reverse_prefix="$work/reverse.forward"
{
  echo '>toy_reverse'
  grep -v '^>' "$work/ref.fa" | tr -d '\n' | rev
  echo
} > "$work/reverse.fa"
"$bwa" fa2pac -f "$work/reverse.fa" "$reverse_prefix" >/dev/null
"$bwa" pac2bwt "$reverse_prefix.pac" "$reverse_prefix.bwt" >/dev/null
"$bwa" bwtupdate "$reverse_prefix.bwt" >/dev/null
"$test_search" "$reverse_prefix.bwt" "$work/reverse.fa"

# Search-only must not require an SA and must not emit SAM records.
cp "$prefix.bwt" "$work/no_sa.bwt"
cp "$prefix.ann" "$work/no_sa.ann"
cp "$prefix.amb" "$work/no_sa.amb"
cp "$prefix.pac" "$work/no_sa.pac"
"$fmsearch" -i "$work/no_sa" -q "$work/reads.fq" -n 0 --no-progress >"$work/search.stdout" 2>"$work/search.stderr"
grep -F 'fmsearch: exact hits=8 located=0' "$work/search.stderr" >/dev/null
test ! -s "$work/search.stdout"

# The experimental rank backend and its seeded path must retain the same
# all-hit count without requiring a suffix array.
"$fmsearch" -i "$work/no_sa" -q "$work/reads.fq" -n 0 \
  --rank-backend batmis --no-progress --timing \
  >"$work/batmis.stdout" 2>"$work/batmis.stderr"
grep -F 'fmsearch: exact hits=8 located=0' "$work/batmis.stderr" >/dev/null
grep -F 'rank_backend_setup=' "$work/batmis.stderr" >/dev/null
"$fmsearch" -i "$work/no_sa" -q "$work/reads.fq" -n 0 \
  --rank-backend batmis --batmis-singleton-fastpath --no-progress \
  >"$work/batmis-fast.stdout" 2>"$work/batmis-fast.stderr"
grep -F 'fmsearch: exact hits=8 located=0' "$work/batmis-fast.stderr" >/dev/null
"$fmsearch" -i "$work/no_sa" -q "$work/reads.fq" -n 0 \
  --rank-backend batmis --seeded-6mer --no-progress \
  >"$work/batmis-seeded.stdout" 2>"$work/batmis-seeded.stderr"
grep -F 'fmsearch: exact hits=8 located=0' "$work/batmis-seeded.stderr" >/dev/null
"$fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 \
  --rank-backend batmis --locate --no-progress \
  >"$work/batmis-locate.stdout" 2>"$work/batmis-locate.stderr"
grep -F 'fmsearch: exact hits=8 located=8' "$work/batmis-locate.stderr" >/dev/null

# SA loading does not locate; locate counts same exact occurrences; SAM emits all rows.
"$fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 --load-sa --no-progress --timing >"$work/load.stdout" 2>"$work/load.stderr"
grep -F 'fmsearch: exact hits=8 located=0' "$work/load.stderr" >/dev/null
"$fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 --locate --no-progress >"$work/locate.stdout" 2>"$work/locate.stderr"
grep -F 'fmsearch: exact hits=8 located=8' "$work/locate.stderr" >/dev/null
"$fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 --sam-output --timing -o "$work/out.sam" >"$work/sam.stdout" 2>"$work/sam.stderr"
grep -F 'fmsearch: exact hits=8 located=8' "$work/sam.stderr" >/dev/null
grep -F 'input_parse=' "$work/sam.stderr" >/dev/null

awk -F '\t' '
    /^@/ { if ($1 == "@SQ" && $2 == "SN:toy" && $3 == "LN:28") header=1; next }
    {
        records++
        if ($1 == "repeat") {
            if ($3 != "toy" || $6 != "4M" || $10 != "ACGT" || $11 != "abcd" || $12 != "NM:i:0" || $13 != "MD:Z:4" || $14 != "NH:i:3") exit 10
            if ($2 == 0) repeat_primary++
            else if ($2 == 256) repeat_secondary++
            else exit 16
            repeat[$4]=1; repeat_count++
        } else if ($1 == "boundary") {
            if ($3 != "toy" || $6 != "4M" || $10 != "ACGT" || $11 != "efgh" || $12 != "NM:i:0" || $13 != "MD:Z:4" || $14 != "NH:i:3") exit 11
            boundary[$4]=1; boundary_count++
        } else if ($1 == "unique") {
            if ($2 != 0 || $3 != "toy" || $4 != 21 || $5 != 60 || $6 != "4M" || $10 != "TTTT" || $11 != "ijkl" || $12 != "NM:i:0" || $13 != "MD:Z:4" || $14 != "NH:i:1") exit 12
            unique_count++
        } else if ($1 == "orientation") {
            if ($2 != 0 || $3 != "toy" || $4 != 3 || $5 != 60 || $6 != "7M" || $10 != "GTACGTA" || $11 != "abcdefg" || $12 != "NM:i:0" || $13 != "MD:Z:7" || $14 != "NH:i:1") exit 17
            orientation_count++
        } else if ($1 == "unmapped") {
            if ($2 != 4 || $3 != "*" || $4 != 0 || $6 != "*" || $10 != "AAAAAC" || $11 != "mnopqr") exit 13
            unmapped_count++
        } else exit 14
    }
    END {
        if (!header || records != 9 || repeat_count != 3 || boundary_count != 3 || repeat_primary != 1 || repeat_secondary != 2 ||
            !repeat[1] || !repeat[5] || !repeat[25] || !boundary[1] || !boundary[5] ||
            !boundary[25] || unique_count != 1 || orientation_count != 1 || unmapped_count != 1) exit 15
    }
' "$work/out.sam"

if [[ -n "$run_loop_fmsearch" ]]; then
  "$run_loop_fmsearch" -i "$work/no_sa" -q "$work/reads.fq" -n 0 \
    --no-progress --timing >"$work/run-loop-search.stdout" 2>"$work/run-loop-search.stderr"
  grep -F 'fmsearch: exact hits=8 located=0' "$work/run-loop-search.stderr" >/dev/null
  "$run_loop_fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 \
    --sam-output --timing --no-progress -o "$work/run-loop.sam" \
    >"$work/run-loop-sam.stdout" 2>"$work/run-loop-sam.stderr"
  grep -F 'fmsearch: exact hits=8 located=8' "$work/run-loop-sam.stderr" >/dev/null
  cmp "$work/out.sam" "$work/run-loop.sam"
fi

if "$fmsearch" -i "$prefix" -q "$work/reads.fq" -n 0 --sam-output --no-progress \
    -o "$work/missing/out.sam" >"$work/bad.out" 2>"$work/bad.err"; then
    echo "Expected unwritable SAM path to fail" >&2
    exit 1
fi
grep -F 'Could not open SAM output' "$work/bad.err" >/dev/null
echo "Forward-only SAM integration checks passed (8 mapped hits, 1 unmapped read)"
