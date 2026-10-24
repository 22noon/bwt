#pragma once

#include "bwa_fmd_index.hpp"
#include <array>
#include <string>
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
#include <vector>
#endif

// Pure exact search returning one SA interval.
SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);

#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
// Test-only trace from the dedicated singleton-run benchmark variant.
SA_Range bench_search_exact_unidirectional_trace(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction,
    std::vector<SA_Range>& trace);

// Test-only differential hook for the native multi-row left-extension body.
SA_Range bench_multirow_extend_left(
    const BwaFMDIndex& index,
    const SA_Range& range,
    uint8_t c);
#endif

#if defined(FM_BENCH_SINGLETON_NOINLINE) || defined(FM_BENCH_SINGLETON_ALWAYS_INLINE)
// Exposed only by opt-in benchmark builds for step-by-step differential checks.
SA_Range bench_singleton_extend_left(
    const BwaFMDIndex& index,
    const SA_Range& range,
    uint8_t c);
#endif

// Optional 6-mer seed accelerator for forward-only exact search. Keys encode
// bases left-to-right using A=0, C=1, G=2, T=3, two bits per base.
class ExactSearch6merTable {
public:
    explicit ExactSearch6merTable(const BwaFMDIndex& index);
    SA_Range search(const BwaFMDIndex& index, const std::string& read) const;

private:
    const BwaFMDIndex* index_ = nullptr;
    std::array<Interval, 4096> intervals_{};
};

// Pure exact search returning a paired FMD interval.
SA_Range search_exact_bidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);

// Extend an existing single SA interval by the supplied sequence.
// Direction::Right consumes extension left-to-right.
// Direction::Left consumes extension right-to-left.
SA_Range extend_exact_unidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction);

// Extend an existing paired FMD interval by the supplied sequence.
// Uses singleton-cardinality and sentinel-aware single-path
// optimizations before falling back to general FMD extension.
SA_Range extend_exact_bidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction);

