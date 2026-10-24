#include "bwa_fmd_index.hpp"
#include "exact_search.hpp"
#include "bwt.h"

#include <cctype>
#include <cstdint>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
uint64_t traced_multirow_extensions = 0;
uint64_t traced_multirow_continuations = 0;
uint64_t traced_multirow_to_singleton = 0;
uint64_t traced_multirow_to_empty = 0;
uint64_t multirow_interval_extensions = 0;
uint64_t multirow_row_zero_cases = 0;
uint64_t multirow_primary_spanning_cases = 0;
#endif

std::string read_reference(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open FASTA: " + path);
    std::string line, sequence;
    bool saw_header = false, saw_sequence = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            if (saw_header) throw std::runtime_error("FASTA must contain exactly one contig");
            saw_header = true;
            continue;
        }
        if (line.empty()) continue;
        if (!saw_header) throw std::runtime_error("FASTA sequence precedes header");
        for (unsigned char raw : line) {
            if (std::isspace(raw)) continue;
            char c = static_cast<char>(std::toupper(raw));
            if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                throw std::runtime_error("FASTA must contain only A/C/G/T");
            sequence.push_back(c);
            saw_sequence = true;
        }
    }
    if (!saw_header || !saw_sequence) throw std::runtime_error("FASTA has no contig sequence");
    return sequence;
}

uint64_t count_occurrences(const std::string& text, const std::string& query)
{
    uint64_t count = 0;
    if (query.size() > text.size()) return 0;
    for (size_t i = 0; i + query.size() <= text.size(); ++i)
        if (text.compare(i, query.size(), query) == 0) ++count;
    return count;
}

std::vector<uint64_t> occurrence_positions(const std::string& text, const std::string& query)
{
    std::vector<uint64_t> positions;
    if (query.size() > text.size()) return positions;
    for (size_t i = 0; i + query.size() <= text.size(); ++i)
        if (text.compare(i, query.size(), query) == 0) positions.push_back(i);
    return positions;
}

bool same_interval(const SA_Range& a, const SA_Range& b)
{
    if (a.size() != b.size()) return false;
    if (a.empty()) return true; // Empty intervals have no canonical insertion point.
    const Interval ai = a.interval();
    const Interval bi = b.interval();
    return ai.l == bi.l && ai.r == bi.r;
}

// Use for two BWA paths, where empty insertion points are part of the
// extension result and must remain identical step by step.
bool same_single_interval_bounds(const SA_Range& a, const SA_Range& b)
{
    if (!a.is_single() || !b.is_single() ||
        a.orientation() != b.orientation())
        return false;
    const Interval ai = a.interval();
    const Interval bi = b.interval();
    return ai.l == bi.l && ai.r == bi.r;
}

bool check(const BwaFMDIndex& index, const BwaFMDIndex& batmis_index,
           const BwaFMDIndex& batmis_fast_index,
           const std::string& reference,
           const ExactSearch6merTable& seeded,
           const ExactSearch6merTable& seeded_batmis,
           const ExactSearch6merTable& seeded_batmis_fast,
           const std::string& query,
           uint64_t& tested, bool check_locations)
{
    const auto expected_positions = occurrence_positions(reference, query);
    const SA_Range range = search_exact_unidirectional(index, query, Direction::Left);
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
    std::vector<SA_Range> run_loop_trace;
    const SA_Range traced_range = bench_search_exact_unidirectional_trace(
        index, query, Direction::Left, run_loop_trace);
    if (!same_single_interval_bounds(range, traced_range)) {
        std::cerr << "Run-loop traced result differs from normal run-loop search for "
                  << query << '\n';
        return false;
    }
#endif
    const SA_Range seeded_range = seeded.search(index, query);
    const SA_Range batmis_range = search_exact_unidirectional(
        batmis_index, query, Direction::Left);
    const SA_Range batmis_fast_range = search_exact_unidirectional(
        batmis_fast_index, query, Direction::Left);
    const SA_Range batmis_seeded_range = seeded_batmis.search(
        batmis_index, query);
    const SA_Range batmis_fast_seeded_range = seeded_batmis_fast.search(
        batmis_fast_index, query);
    const uint64_t expected = expected_positions.size();
    const uint64_t observed = range.size();
    ++tested;
    if (!same_interval(range, batmis_range) ||
        !same_interval(range, batmis_seeded_range) ||
        !same_interval(range, batmis_fast_range) ||
        !same_interval(range, batmis_fast_seeded_range)) {
        const Interval expected = range.interval();
        const Interval observed_batmis = batmis_range.interval();
        const Interval observed_seeded = batmis_seeded_range.interval();
        std::cerr << "Backend interval mismatch for query " << query
                  << ": BWA [" << expected.l << ',' << expected.r
                  << "), BatMIS [" << observed_batmis.l << ','
                  << observed_batmis.r << "), seeded BatMIS ["
                  << observed_seeded.l << ',' << observed_seeded.r << ")\n";
        std::cerr << "BatMIS fast-path interval [" << batmis_fast_range.interval().l
                  << ',' << batmis_fast_range.interval().r << ")\n";
        return false;
    }
    // Compare each extension step, not only the final occurrence count.
    SA_Range step_bwa = index.initial_single_range(
        query.back() == 'A' ? 0 : query.back() == 'C' ? 1 : query.back() == 'G' ? 2 : 3,
        SA_Range::Orientation::Left);
    SA_Range step_batmis = batmis_index.initial_single_range(
        query.back() == 'A' ? 0 : query.back() == 'C' ? 1 : query.back() == 'G' ? 2 : 3,
        SA_Range::Orientation::Left);
    SA_Range step_batmis_fast = batmis_fast_index.initial_single_range(
        query.back() == 'A' ? 0 : query.back() == 'C' ? 1 : query.back() == 'G' ? 2 : 3,
        SA_Range::Orientation::Left);
    step_bwa.set_position(0);
    step_bwa.set_valid(!step_bwa.empty());
    step_bwa.set_unidirectional(true);
    std::vector<SA_Range> generic_trace{step_bwa};
    for (size_t pos = query.size() - 1; pos > 0 && !step_bwa.empty(); --pos) {
        const char base = query[pos - 1];
        const uint8_t c = base == 'A' ? 0 : base == 'C' ? 1 : base == 'G' ? 2 : 3;
        step_bwa = index.extend_left(step_bwa, c);
        step_bwa.set_position(static_cast<uint32_t>(query.size() - pos));
        step_bwa.set_valid(!step_bwa.empty());
        step_bwa.set_unidirectional(true);
        generic_trace.push_back(step_bwa);
        step_batmis = batmis_index.extend_left(step_batmis, c);
        step_batmis_fast = batmis_fast_index.extend_left(step_batmis_fast, c);
        if (!same_interval(step_bwa, step_batmis)) {
            std::cerr << "Step interval mismatch for query " << query
                      << " at position " << (pos - 1) << '\n';
            return false;
        }
        if (!same_interval(step_batmis, step_batmis_fast)) {
            std::cerr << "BatMIS singleton fast-path mismatch for query " << query
                      << " at position " << (pos - 1) << '\n';
            return false;
        }
    }
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
    if (run_loop_trace.size() != generic_trace.size()) {
        std::cerr << "Run-loop/generic step-count mismatch for " << query
                  << ": run-loop " << run_loop_trace.size()
                  << ", generic " << generic_trace.size() << '\n';
        return false;
    }
    for (size_t step = 1; step < generic_trace.size(); ++step) {
        const SA_Range& before = generic_trace[step - 1];
        const SA_Range& after = generic_trace[step];
        if (before.size() > 1) {
            ++traced_multirow_extensions;
            if (after.empty())
                ++traced_multirow_to_empty;
            else if (after.size() == 1)
                ++traced_multirow_to_singleton;
            else
                ++traced_multirow_continuations;
        }
    }
    for (size_t step = 0; step < generic_trace.size(); ++step) {
        const SA_Range& actual = run_loop_trace[step];
        const SA_Range& expected_step = generic_trace[step];
        if (!same_single_interval_bounds(actual, expected_step) ||
            actual.index() != expected_step.index() ||
            actual.is_single() != expected_step.is_single() ||
            actual.orientation() != expected_step.orientation() ||
            actual.position() != expected_step.position() ||
            actual.valid() != expected_step.valid() ||
            actual.is_unidirectional() != expected_step.is_unidirectional()) {
            std::cerr << "Run-loop/generic state mismatch for " << query
                      << " at extension step " << step << ": run-loop ["
                      << actual.interval().l << ',' << actual.interval().r
                      << "), generic [" << expected_step.interval().l << ','
                      << expected_step.interval().r << ")\n";
            return false;
        }
    }
#endif
    if (!same_single_interval_bounds(range, step_bwa)) {
        std::cerr << "Production exact-search interval differs from generic "
                  << "extend_left() baseline for query " << query
                  << ": optimized [" << range.interval().l << ','
                  << range.interval().r << "), generic ["
                  << step_bwa.interval().l << ',' << step_bwa.interval().r
                  << ")\n";
        return false;
    }
    if (range.interval().l != seeded_range.interval().l ||
        range.interval().r != seeded_range.interval().r) {
        std::cerr << "Seeded interval mismatch for query " << query
                  << ": baseline [" << range.interval().l << ',' << range.interval().r
                  << "), seeded [" << seeded_range.interval().l << ','
                  << seeded_range.interval().r << ")\n";
        return false;
    }
    if (expected != observed) {
        const Interval bounds = range.is_single() ? range.interval() : range.primary_interval();
        std::cerr << "Mismatch for query " << query << ": expected " << expected
                  << ", observed " << observed << ", interval [" << bounds.l
                  << ',' << bounds.r << ")\n";
        return false;
    }
    if (check_locations) {
        auto observed_positions = index.locate_forward(range);
        auto batmis_positions = batmis_index.locate_forward(batmis_range);
        auto fast_positions = batmis_fast_index.locate_forward(batmis_fast_range);
        std::sort(observed_positions.begin(), observed_positions.end());
        std::sort(batmis_positions.begin(), batmis_positions.end());
        std::sort(fast_positions.begin(), fast_positions.end());
        if (observed_positions != expected_positions ||
            batmis_positions != expected_positions ||
            fast_positions != expected_positions) {
            const Interval bounds = range.is_single() ? range.interval() : range.primary_interval();
            auto print = [](const std::vector<uint64_t>& positions) {
                std::cerr << '[';
                for (size_t i = 0; i < positions.size(); ++i)
                    std::cerr << (i ? ", " : "") << positions[i];
                std::cerr << ']';
            };
            std::cerr << "Location mismatch for query " << query << ": expected positions ";
            print(expected_positions);
            std::cerr << ", observed positions ";
            print(observed_positions);
            std::cerr << ", BatMIS positions ";
            print(batmis_positions);
            std::cerr << ", BatMIS fast-path positions ";
            print(fast_positions);
            std::cerr << ", interval [" << bounds.l << ',' << bounds.r << ")\n";
            return false;
        }
    }
    return true;
}

bool enumerate_dna(const BwaFMDIndex& index, const BwaFMDIndex& batmis_index,
                   const BwaFMDIndex& batmis_fast_index,
                   const std::string& reference,
                   const ExactSearch6merTable& seeded,
                   const ExactSearch6merTable& seeded_batmis,
                   const ExactSearch6merTable& seeded_batmis_fast,
                   std::string& query,
                   size_t target, uint64_t& tested,
                   bool check_locations)
{
    static constexpr char bases[] = "ACGT";
    if (query.size() == target) return check(index, batmis_index, batmis_fast_index, reference,
        seeded, seeded_batmis, seeded_batmis_fast, query, tested, check_locations);
    for (char base : bases) {
        if (!base) break;
        query.push_back(base);
        if (!enumerate_dna(index, batmis_index, batmis_fast_index, reference, seeded,
                           seeded_batmis, seeded_batmis_fast, query, target, tested,
                           check_locations)) return false;
        query.pop_back();
    }
    return true;
}
}

int main(int argc, char** argv)
{
    if (argc != 3 && argc != 4) {
        std::cerr << "Usage: test_forward_only_search <index.bwt> <reference.fa> [index.sa]\n";
        return 2;
    }
    bwt_t* bwt = nullptr;
    try {
        const std::string reference = read_reference(argv[2]);
        bwt = bwt_restore_bwt(argv[1]);
        if (!bwt) throw std::runtime_error(std::string("cannot restore BWT: ") + argv[1]);
        if (static_cast<uint64_t>(bwt->seq_len) != reference.size())
            throw std::runtime_error("BWT sequence length does not match reference length");
        const bool check_locations = argc == 4;
        if (check_locations) bwt_restore_sa(argv[3], bwt);
        BwaFMDIndex index(bwt);
        BwaFMDIndex batmis_index(bwt, RankBackend::BatMIS);
        BwaFMDIndex batmis_fast_index(bwt, RankBackend::BatMIS, true);
        const ExactSearch6merTable seeded(index);
        const ExactSearch6merTable seeded_batmis(batmis_index);
        const ExactSearch6merTable seeded_batmis_fast(batmis_fast_index);
        for (uint8_t c = 0; c < 4; ++c) {
            const SA_Range a = index.initial_single_range(
                c, SA_Range::Orientation::Left);
            const SA_Range b = batmis_index.initial_single_range(
                c, SA_Range::Orientation::Left);
            if (!same_interval(a, b))
                throw std::runtime_error("initial interval mismatch");
        }
        uint64_t single_steps = 0;
        for (uint64_t row = 0; row <= bwt->seq_len; ++row) {
            for (uint8_t c = 0; c < 4; ++c) {
                const SA_Range a = index.extend_left(
                    SA_Range::single(&index, {row, row + 1}), c);
                const SA_Range b = batmis_index.extend_left(
                    SA_Range::single(&batmis_index, {row, row + 1}), c);
                const SA_Range fast = batmis_fast_index.extend_left(
                    SA_Range::single(&batmis_fast_index, {row, row + 1}), c);
                if (!same_interval(a, b) || !same_interval(b, fast))
                    throw std::runtime_error("singleton/primary-row extension mismatch");
                ++single_steps;
            }
        }
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
        const uint64_t row_count = static_cast<uint64_t>(bwt->seq_len) + 1;
        for (uint64_t left = 0; left < row_count; ++left) {
            for (uint64_t right = left + 2; right <= row_count; ++right) {
                const SA_Range parent = SA_Range::single(
                    &index, {left, right}, SA_Range::Orientation::Left);
                const bool spans_primary =
                    left <= static_cast<uint64_t>(bwt->primary) &&
                    static_cast<uint64_t>(bwt->primary) < right;
                for (uint8_t c = 0; c < 4; ++c) {
                    const SA_Range expected = index.extend_left(parent, c);
                    const SA_Range actual =
                        bench_multirow_extend_left(index, parent, c);
                    if (!same_single_interval_bounds(actual, expected) ||
                        actual.index() != expected.index() ||
                        actual.is_single() != expected.is_single() ||
                        actual.orientation() != expected.orientation() ||
                        actual.position() != expected.position() ||
                        actual.valid() != expected.valid() ||
                        actual.is_unidirectional() != expected.is_unidirectional()) {
                        std::cerr << "Multi-row extension mismatch for interval ["
                                  << left << ',' << right << "), base "
                                  << static_cast<unsigned>(c) << ": expected ["
                                  << expected.interval().l << ','
                                  << expected.interval().r << "), actual ["
                                  << actual.interval().l << ','
                                  << actual.interval().r << ")\n";
                        bwt_destroy(bwt);
                        return 1;
                    }
                    ++multirow_interval_extensions;
                    multirow_row_zero_cases += left == 0;
                    multirow_primary_spanning_cases += spans_primary;
                }
            }
        }
#endif
        uint64_t tested = 0;
        std::string query;
        for (size_t length = 1; length <= 6; ++length)
            if (!enumerate_dna(index, batmis_index, batmis_fast_index, reference,
                               seeded, seeded_batmis, seeded_batmis_fast,
                               query, length, tested,
                               check_locations)) {
                bwt_destroy(bwt);
                return 1;
            }
        for (size_t length = 1; length <= 32 && length <= reference.size(); ++length)
            for (size_t start = 0; start + length <= reference.size(); ++start)
                if (!check(index, batmis_index, batmis_fast_index, reference, seeded,
                           seeded_batmis, seeded_batmis_fast,
                           reference.substr(start, length),
                           tested, check_locations)) {
                    bwt_destroy(bwt);
                    return 1;
                }

        // Include a query with a present final 6-mer seed but an absent prefix.
        bool tested_seed_present_miss = false;
        if (reference.size() >= 7) {
            for (size_t start = 0; start + 7 <= reference.size() && !tested_seed_present_miss; ++start) {
                std::string candidate = reference.substr(start, 7);
                static constexpr char bases[] = "ACGT";
                for (size_t b = 0; b < 4; ++b) {
                    candidate[0] = bases[b];
                    if (count_occurrences(reference, candidate) == 0 &&
                        count_occurrences(reference, candidate.substr(1)) > 0) {
                        if (!check(index, batmis_index, batmis_fast_index, reference, seeded,
                                   seeded_batmis, seeded_batmis_fast,
                                   candidate, tested,
                                   check_locations)) {
                            bwt_destroy(bwt);
                            return 1;
                        }
                        tested_seed_present_miss = true;
                        break;
                    }
                }
            }
        }
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
        if (!traced_multirow_extensions || !traced_multirow_continuations ||
            !traced_multirow_to_singleton || !traced_multirow_to_empty ||
            !multirow_row_zero_cases || !multirow_primary_spanning_cases) {
            throw std::runtime_error(
                "forward-only fixture did not cover all multi-row extension states");
        }
#endif
        std::cout << "Passed: " << tested << " queries tested; "
                  << single_steps << " singleton-symbol extensions checked\n";
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
        std::cout << "Multi-row differential: " << multirow_interval_extensions
                  << " interval/base requests; " << traced_multirow_extensions
                  << " traced steps (" << traced_multirow_continuations
                  << " remain multi-row, " << traced_multirow_to_singleton
                  << " become singleton, " << traced_multirow_to_empty
                  << " become empty)\n";
#endif
        bwt_destroy(bwt);
        return 0;
    } catch (const std::exception& e) {
        if (bwt) bwt_destroy(bwt);
        std::cerr << "test_forward_only_search: " << e.what() << '\n';
        return 2;
    }
}
