# `fmsearch` command-line runner

This adds a small CLI around the current `Run_Mismatch_Pipeline` API. It does
not build BWA indexes; use an existing BWA index prefix.

## Build

Copy `include/mismatch_pipeline.hpp` to the project's `include/` directory and
`src/mismatch_pipeline.cpp` plus `src/fmsearch.cpp` to `src/`. Apply
`CMakeLists_cli.patch` to the repository's current `CMakeLists.txt`, then build:

```bash
cmake -S . -B build
cmake --build build --target fmsearch -j
```

The repository's normal CMake configuration may require `bwa` on `PATH` to
create test fixtures. No index-builder target is added here.

## Run

```bash
# Default: up to one mismatch, bidirectional zero-mismatch mode, first hit only
./build/fmsearch /path/to/index_prefix reads.fa

# Exact-only search
./build/fmsearch /path/to/index_prefix reads.fa --max-mismatches 0

# Try the Fast zero-mismatch path
./build/fmsearch /path/to/index_prefix reads.fa --mode fast

# Report all exact and one-mismatch hits
./build/fmsearch /path/to/index_prefix reads.fa --all-hits
```

`index_prefix` is the common prefix for BWA files such as `index_prefix.bwt`,
`index_prefix.sa`, `index_prefix.ann`, `index_prefix.amb`, and `index_prefix.pac`.
The `.bwt` and `.sa` files are required. If reference metadata files are present,
output uses reference name, 1-based position, and strand; otherwise it prints raw
SA coordinates in the indexed text.

Input may be FASTA, standard four-line FASTQ, or plain text with one sequence per
line. Sequences must contain only A/C/G/T (case-insensitive). The current
pipeline supports maximum mismatch counts 0 and 1 only.

Output is tab-separated. `mismatch_positions_0based` lists zero-based query
positions for reported mismatches; the field is empty for exact hits.
