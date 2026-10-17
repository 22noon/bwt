#include "exact_search.hpp"
#include "mismatch_search.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

uint8_t code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
        std::cerr << "Invalid DNA base: " << c << '\n';
        std::abort();
    }
}

SA_Range exact_range(
    const BwaFMDIndex& index,
    const std::string& pattern,
    Direction direction)
{
    assert(!pattern.empty());

    SA_Range range;
    if (direction == Direction::Right) {
        range = index.initial_range(code(pattern.front()));
        for (size_t i = 1; i < pattern.size(); ++i)
            range = index.extend_right(range, code(pattern[i]));
    } else {
        range = index.initial_range(code(pattern.back()));
        for (size_t i = pattern.size() - 1; i-- > 0;)
            range = index.extend_left(range, code(pattern[i]));
    }
    return range;
}

void generate_patterns(
    size_t length,
    std::string& current,
    std::vector<std::string>& output)
{
    if (current.size() == length) {
        output.push_back(current);
        return;
    }

    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    for (char c : bases) {
        current.push_back(c);
        generate_patterns(length, current, output);
        current.pop_back();
    }
}

unsigned hamming_distance(
    const std::string& a,
    const std::string& b)
{
    assert(a.size() == b.size());
    unsigned distance = 0;
    for (size_t i = 0; i < a.size(); ++i)
        distance += (a[i] != b[i]);
    return distance;
}

std::string reverse_complement(const std::string& sequence)
{
    std::string result;
    result.reserve(sequence.size());
    for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
        switch (*it) {
        case 'A': result.push_back('T'); break;
        case 'C': result.push_back('G'); break;
        case 'G': result.push_back('C'); break;
        case 'T': result.push_back('A'); break;
        default: std::abort();
        }
    }
    return result;
}

std::string read_single_fasta(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open reference FASTA: " + path);
    std::string sequence, line;
    bool saw_header = false;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (saw_header)
                throw std::runtime_error("expected a single-record reference FASTA");
            saw_header = true;
            continue;
        }
        for (char c : line) {
            if (c == ' ' || c == '\t' || c == '\r') continue;
            if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
            code(c); // Reject unsupported reference symbols in this fixture.
            sequence.push_back(c);
        }
    }
    if (!saw_header || sequence.empty())
        throw std::runtime_error("reference FASTA is missing a header or sequence");
    return sequence;
}

std::vector<uint64_t> naive_positions(
    const std::string& text, const std::string& pattern)
{
    std::vector<uint64_t> positions;
    if (pattern.size() > text.size()) return positions;
    for (size_t i = 0; i + pattern.size() <= text.size(); ++i)
        if (text.compare(i, pattern.size(), pattern) == 0)
            positions.push_back(static_cast<uint64_t>(i));
    return positions;
}

bool same_interval(const Interval& a, const Interval& b)
{
    return a.l == b.l && a.r == b.r;
}

bool same_search_result(const SA_Range& a, const SA_Range& b)
{
    if (!a.is_bidirectional() || !b.is_bidirectional() ||
        !same_interval(a.primary_interval(), b.primary_interval()) ||
        !same_interval(a.companion_interval(), b.companion_interval()) ||
        a.valid() != b.valid() || a.position() != b.position() ||
        a.mismatch_count() != b.mismatch_count() ||
        a.forbidden_base() != b.forbidden_base() ||
        a.is_unidirectional() != b.is_unidirectional() ||
        a.mismatches().size() != b.mismatches().size())
        return false;
    for (size_t i = 0; i < a.mismatches().size(); ++i)
        if (a.mismatches()[i].Pos != b.mismatches()[i].Pos ||
            a.mismatches()[i].Char != b.mismatches()[i].Char)
            return false;
    return true;
}

bool same_ordered_results(
    const std::vector<SA_Range>& a, const std::vector<SA_Range>& b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!same_search_result(a[i], b[i])) return false;
    return true;
}

std::string scanned_prefix_to_pattern(
    const std::string& scanned_prefix, Direction direction)
{
    if (direction == Direction::Right) return scanned_prefix;
    return std::string(scanned_prefix.rbegin(), scanned_prefix.rend());
}

SA_Range range_for_scanned_prefix(
    const BwaFMDIndex& index,
    const std::string& query,
    const std::string& scanned_prefix,
    const std::vector<char>& observed_bases,
    Direction direction,
    bool capture_reference_bases,
    char forbidden_base = '\0')
{
    SA_Range range = exact_range(
        index, scanned_prefix_to_pattern(scanned_prefix, direction), direction);
    range.set_position(static_cast<uint32_t>(scanned_prefix.size() - 1));
    range.set_valid(!range.empty());
    range.set_forbidden_base(forbidden_base);
    for (size_t i = 1; i < observed_bases.size(); ++i) {
        if (observed_bases[i] == query[static_cast<size_t>(
                direction == Direction::Right ? i : query.size() - i - 1)])
            continue;
        if (capture_reference_bases) {
            const uint32_t original_one_based = direction == Direction::Right
                ? static_cast<uint32_t>(i + 1)
                : static_cast<uint32_t>(query.size() - i);
            range.add_mismatch(original_one_based, observed_bases[i]);
        } else {
            range.add_mismatch(
                static_cast<uint32_t>(i),
                query[static_cast<size_t>(
                    direction == Direction::Right ? i : query.size() - i - 1)]);
        }
    }
    return range;
}

struct SearchOracleCoverage {
    uint64_t threshold_states = 0;
    uint64_t unchanged_cardinality = 0;
    uint64_t changed_cardinality = 0;
    uint64_t empty_expected_extension = 0;
    uint64_t sentinel_adjusted_full_extension = 0;
    uint64_t partial_without_sentinel = 0;
    uint64_t forbidden_threshold_rejects = 0;
    uint64_t max_hits_one_stops = 0;
    uint64_t max_hits_two_stops = 0;
    std::string first_sentinel_full_query;
    std::string first_sentinel_full_prefix;
    Direction first_sentinel_full_direction = Direction::Right;
    Interval first_sentinel_parent_primary;
    Interval first_sentinel_parent_companion;
    Interval first_sentinel_child_primary;
    Interval first_sentinel_child_companion;
    uint8_t first_sentinel_nonempty_branches = 0;
    std::string first_non_sentinel_partial_query;
    std::string first_non_sentinel_partial_prefix;
    Direction first_non_sentinel_partial_direction = Direction::Right;
    Interval first_partial_parent_primary;
    Interval first_partial_parent_companion;
    Interval first_partial_child_primary;
    Interval first_partial_child_companion;
    uint8_t first_partial_nonempty_branches = 0;
};

void brute_force_search_oracle(
    const BwaFMDIndex& index,
    const std::string& query,
    const std::string& indexed_text,
    unsigned mismatch_budget,
    int max_hits,
    Direction direction,
    bool capture_reference_bases,
    std::vector<SA_Range>& expected_hits,
    std::vector<SA_Range>& expected_rejects,
    SearchOracleCoverage& coverage)
{
    struct State {
        std::string scanned;
        std::vector<char> observed;
    };

    const char first = direction == Direction::Right
        ? query.front() : query.back();
    std::vector<State> stack{{std::string(1, first), std::vector<char>{first}}};
    uint64_t hits_found = 0;
    const auto hit_limit_reached = [&]() {
        if (max_hits <= 0 || hits_found < static_cast<uint64_t>(max_hits))
            return false;
        if (max_hits == 1) ++coverage.max_hits_one_stops;
        if (max_hits == 2) ++coverage.max_hits_two_stops;
        return true;
    };

    const auto push_hit = [&](const State& state) {
        const std::string pattern =
            scanned_prefix_to_pattern(state.scanned, direction);
        SA_Range hit = range_for_scanned_prefix(
            index, query, state.scanned, state.observed, direction,
            capture_reference_bases);
        const std::vector<uint64_t> expected_positions =
            naive_positions(indexed_text, pattern);
        std::vector<uint64_t> actual_positions = index.locate(hit);
        std::sort(actual_positions.begin(), actual_positions.end());
        if (actual_positions != expected_positions ||
            hit.size() != expected_positions.size() ||
            hamming_distance(query, pattern) != mismatch_budget)
            throw std::runtime_error(
                "brute-force hit oracle disagrees for pattern " + pattern);
        expected_hits.push_back(std::move(hit));
        ++hits_found;
    };

    while (!stack.empty()) {
        State current = std::move(stack.back());
        stack.pop_back();
        const size_t position = current.scanned.size() - 1;
        if (position + 1 >= query.size()) continue;

        const char expected_base = query[static_cast<size_t>(
            direction == Direction::Right ? position + 1
                                          : query.size() - position - 2)];
        unsigned mismatch_count = 0;
        for (size_t i = 1; i < current.observed.size(); ++i) {
            const char query_base = query[static_cast<size_t>(
                direction == Direction::Right ? i : query.size() - i - 1)];
            mismatch_count += current.observed[i] != query_base;
        }
        if (mismatch_count == mismatch_budget) {
            ++coverage.threshold_states;
            const SA_Range parent = range_for_scanned_prefix(
                index, query, current.scanned, current.observed, direction,
                capture_reference_bases);
            State next = current;
            next.scanned.push_back(expected_base);
            next.observed.push_back(expected_base);
            const SA_Range exact = range_for_scanned_prefix(
                index, query, next.scanned, next.observed, direction,
                capture_reference_bases);

            const Interval relevant = direction == Direction::Left
                ? parent.primary_interval() : parent.companion_interval();
            const bwtint_t sentinel_row = index.bwt()->primary;
            const bool contains_sentinel =
                relevant.l <= sentinel_row && sentinel_row < relevant.r;
            const uint64_t sentinel_offset = contains_sentinel ? 1 : 0;
            if (parent.size() == exact.size())
                ++coverage.unchanged_cardinality;
            else
                ++coverage.changed_cardinality;
            if (exact.empty()) ++coverage.empty_expected_extension;
            if (contains_sentinel && exact.size() + sentinel_offset == parent.size())
                ++coverage.sentinel_adjusted_full_extension;
            if (contains_sentinel && exact.size() + sentinel_offset == parent.size() &&
                parent.size() != exact.size() &&
                coverage.first_sentinel_full_query.empty()) {
                coverage.first_sentinel_full_query = query;
                coverage.first_sentinel_full_prefix = current.scanned;
                coverage.first_sentinel_full_direction = direction;
                coverage.first_sentinel_parent_primary = parent.primary_interval();
                coverage.first_sentinel_parent_companion = parent.companion_interval();
                coverage.first_sentinel_child_primary = exact.primary_interval();
                coverage.first_sentinel_child_companion = exact.companion_interval();
                const BranchSet branches = parent.branch(direction);
                for (uint8_t c = 0; c < 4; ++c)
                    coverage.first_sentinel_nonempty_branches +=
                        branches.has_branch(c) ? 1 : 0;
            }
            if (!contains_sentinel && exact.size() < parent.size())
                ++coverage.partial_without_sentinel;
            if (!contains_sentinel && exact.size() < parent.size() &&
                !exact.empty() && coverage.first_non_sentinel_partial_query.empty()) {
                coverage.first_non_sentinel_partial_query = query;
                coverage.first_non_sentinel_partial_prefix = current.scanned;
                coverage.first_non_sentinel_partial_direction = direction;
                coverage.first_partial_parent_primary = parent.primary_interval();
                coverage.first_partial_parent_companion = parent.companion_interval();
                coverage.first_partial_child_primary = exact.primary_interval();
                coverage.first_partial_child_companion = exact.companion_interval();
                const BranchSet branches = parent.branch(direction);
                for (uint8_t c = 0; c < 4; ++c)
                    coverage.first_partial_nonempty_branches +=
                        branches.has_branch(c) ? 1 : 0;
            }

            // This intentionally models the current raw-cardinality test.
            // Sentinel-aware semantics are measured separately above.
            if (parent.size() != exact.size()) {
                SA_Range reject = parent;
                reject.set_forbidden_base(expected_base);
                expected_rejects.push_back(std::move(reject));
                ++coverage.forbidden_threshold_rejects;
            }
            if (exact.empty()) continue;
            if (position + 1 == query.size() - 1) {
                push_hit(next);
                if (hit_limit_reached()) return;
            } else {
                stack.push_back(std::move(next));
            }
            continue;
        }

        for (char branch_base : {'A', 'C', 'G', 'T'}) {
            State child = current;
            child.scanned.push_back(branch_base);
            child.observed.push_back(branch_base);
            SA_Range child_range = range_for_scanned_prefix(
                index, query, child.scanned, child.observed, direction,
                capture_reference_bases);
            if (child_range.empty()) continue;

            unsigned child_mismatches = mismatch_count +
                (branch_base != expected_base);
            if (child_mismatches > mismatch_budget) {
                expected_rejects.push_back(std::move(child_range));
                continue;
            }
            if (position + 1 == query.size() - 1) {
                if (child_mismatches == mismatch_budget) {
                    push_hit(child);
                    if (hit_limit_reached()) return;
                }
            } else {
                stack.push_back(std::move(child));
            }
        }
    }
}

void check_complete_search_results(
    const BwaFMDIndex& index,
    const std::string& query,
    const std::string& indexed_text,
    unsigned mismatch_budget,
    int max_hits,
    Direction direction,
    bool capture_reference_bases,
    SearchOracleCoverage& coverage)
{
    const char first = direction == Direction::Right
        ? query.front() : query.back();
    SA_Range seed = index.initial_range(code(first));
    seed.set_position(0);
    seed.set_valid(true);

    std::vector<SA_Range> expected_hits, expected_rejects;
    brute_force_search_oracle(
        index, query, indexed_text, mismatch_budget, max_hits, direction,
        capture_reference_bases, expected_hits, expected_rejects, coverage);

    std::vector<SA_Range> actual_hits, actual_rejects;
    search_unidirectional(
        query, seed, static_cast<int>(mismatch_budget), max_hits,
        actual_hits, actual_rejects, direction, capture_reference_bases);

    if (!same_ordered_results(actual_hits, expected_hits) ||
        !same_ordered_results(actual_rejects, expected_rejects)) {
        std::cerr << "Complete mismatch-search result failure: query=" << query
                  << " direction="
                  << (direction == Direction::Right ? "right" : "left")
                  << " budget=" << mismatch_budget << " max_hits=" << max_hits
                  << " capture=" << capture_reference_bases
                  << " actual hits/rejects=" << actual_hits.size() << '/'
                  << actual_rejects.size() << " expected="
                  << expected_hits.size() << '/' << expected_rejects.size()
                  << '\n';
        std::abort();
    }
}

void test_complete_ordered_results(
    const BwaFMDIndex& index,
    const std::vector<std::string>& queries,
    const std::string& indexed_text)
{
    SearchOracleCoverage coverage;
    for (const std::string& query : queries) {
        for (Direction direction : {Direction::Right, Direction::Left}) {
            for (unsigned budget : {0u, 1u}) {
                for (bool capture : {false, true}) {
                    check_complete_search_results(
                        index, query, indexed_text, budget, 0, direction,
                        capture, coverage);
                }
            }
        }
    }

    // Exercise early termination after the first and second ordered hits.
    for (const std::string& query : {"ACGT", "TGCA", "AAAA", "TTTT"}) {
        for (Direction direction : {Direction::Right, Direction::Left}) {
            for (bool capture : {false, true}) {
                check_complete_search_results(
                    index, query, indexed_text, 1, 1, direction,
                    capture, coverage);
                check_complete_search_results(
                    index, query, indexed_text, 1, 2, direction,
                    capture, coverage);
            }
        }
    }
    // Search the exhaustive set until a real two-hit stop is exercised.
    for (const std::string& query : queries) {
        if (coverage.max_hits_two_stops != 0) break;
        for (Direction direction : {Direction::Right, Direction::Left}) {
            for (bool capture : {false, true}) {
                if (coverage.max_hits_two_stops != 0) break;
                check_complete_search_results(
                    index, query, indexed_text, 1, 2, direction,
                    capture, coverage);
            }
        }
    }

    if (coverage.threshold_states == 0 ||
        coverage.unchanged_cardinality == 0 ||
        coverage.changed_cardinality == 0 ||
        coverage.empty_expected_extension == 0 ||
        coverage.sentinel_adjusted_full_extension == 0 ||
        coverage.partial_without_sentinel == 0 ||
        coverage.forbidden_threshold_rejects == 0 ||
        coverage.max_hits_one_stops == 0 ||
        coverage.max_hits_two_stops == 0) {
        std::cerr << "Threshold coverage incomplete: states="
                  << coverage.threshold_states
                  << " same=" << coverage.unchanged_cardinality
                  << " changed=" << coverage.changed_cardinality
                  << " empty=" << coverage.empty_expected_extension
                  << " sentinel-adjusted-full="
                  << coverage.sentinel_adjusted_full_extension
                  << " partial-without-sentinel="
                  << coverage.partial_without_sentinel
                  << " forbidden-rejects="
                  << coverage.forbidden_threshold_rejects
                  << " max-hit-stops=" << coverage.max_hits_one_stops
                  << '/' << coverage.max_hits_two_stops << '\n';
        std::abort();
    }
    std::cout << "Complete ordered mismatch results passed; threshold states="
              << coverage.threshold_states
              << " unchanged=" << coverage.unchanged_cardinality
              << " changed=" << coverage.changed_cardinality
              << " empty=" << coverage.empty_expected_extension
              << " sentinel-adjusted-full="
              << coverage.sentinel_adjusted_full_extension
              << " partial-without-sentinel="
              << coverage.partial_without_sentinel
              << " forbidden-rejects="
              << coverage.forbidden_threshold_rejects
              << " max-hit-stops=" << coverage.max_hits_one_stops
              << '/' << coverage.max_hits_two_stops << '\n';
    if (!coverage.first_sentinel_full_query.empty()) {
        const auto print_interval = [](const char* label, const Interval& interval) {
            std::cout << ' ' << label << "=[" << interval.l << ','
                      << interval.r << ')';
        };
        std::cout << "Threshold sentinel-only cardinality example: query="
                  << coverage.first_sentinel_full_query
                  << " scanned-prefix=" << coverage.first_sentinel_full_prefix
                  << " direction="
                  << (coverage.first_sentinel_full_direction == Direction::Right
                          ? "right" : "left");
        print_interval("parent-primary", coverage.first_sentinel_parent_primary);
        print_interval("parent-companion", coverage.first_sentinel_parent_companion);
        print_interval("child-primary", coverage.first_sentinel_child_primary);
        print_interval("child-companion", coverage.first_sentinel_child_companion);
        std::cout << " nonempty-branches="
                  << static_cast<unsigned>(coverage.first_sentinel_nonempty_branches);
        std::cout << '\n';
    }
    if (!coverage.first_non_sentinel_partial_query.empty()) {
        const auto print_interval = [](const char* label, const Interval& interval) {
            std::cout << ' ' << label << "=[" << interval.l << ','
                      << interval.r << ')';
        };
        std::cout << "Threshold ordinary partial example: query="
                  << coverage.first_non_sentinel_partial_query
                  << " scanned-prefix=" << coverage.first_non_sentinel_partial_prefix
                  << " direction="
                  << (coverage.first_non_sentinel_partial_direction == Direction::Right
                          ? "right" : "left");
        print_interval("parent-primary", coverage.first_partial_parent_primary);
        print_interval("parent-companion", coverage.first_partial_parent_companion);
        print_interval("child-primary", coverage.first_partial_child_primary);
        print_interval("child-companion", coverage.first_partial_child_companion);
        std::cout << " nonempty-branches="
                  << static_cast<unsigned>(coverage.first_partial_nonempty_branches)
                  << '\n';
    }
}

void test_over_budget_seed_rejections(const BwaFMDIndex& index)
{
    const std::string query = "AC";
    for (Direction direction : {Direction::Right, Direction::Left}) {
        const char expected_seed = direction == Direction::Right
            ? query.front() : query.back();
        const char observed_seed = direction == Direction::Right ? 'C' : 'A';
        const char expected_next = direction == Direction::Right
            ? query.back() : query.front();
        for (bool capture : {false, true}) {
            SA_Range seed = index.initial_range(code(observed_seed));
            seed.set_position(0);
            seed.set_valid(true);
            seed.add_mismatch(
                capture
                    ? static_cast<uint32_t>(direction == Direction::Right
                          ? 1 : query.size())
                    : 0,
                capture ? observed_seed : expected_seed);

            std::vector<SA_Range> actual_hits, actual_rejects;
            search_unidirectional(
                query, seed, 0, 0, actual_hits, actual_rejects,
                direction, capture);
            if (!actual_hits.empty())
                throw std::runtime_error(
                    "over-budget seed unexpectedly produced a hit");

            std::vector<SA_Range> expected_rejects;
            for (char branch_base : {'A', 'C', 'G', 'T'}) {
                const std::string scanned{observed_seed, branch_base};
                const std::string pattern =
                    scanned_prefix_to_pattern(scanned, direction);
                SA_Range rejected = exact_range(index, pattern, direction);
                if (rejected.empty()) continue;
                rejected.set_position(1);
                rejected.set_valid(true);
                rejected.add_mismatch(
                    capture
                        ? static_cast<uint32_t>(direction == Direction::Right
                              ? 1 : query.size())
                        : 0,
                    capture ? observed_seed : expected_seed);
                if (branch_base != expected_next) {
                    if (capture) {
                        rejected.add_mismatch(
                            direction == Direction::Right ? 2 : 1,
                            branch_base);
                    } else {
                        rejected.add_mismatch(1, expected_next);
                    }
                }
                expected_rejects.push_back(std::move(rejected));
            }
            if (!same_ordered_results(actual_rejects, expected_rejects))
                throw std::runtime_error(
                    "over-budget branch rejection differs from expected ordered states");
        }
    }
}

void test_threshold_cardinality_examples(const BwaFMDIndex& index)
{
    struct Example {
        std::string query;
        std::string current_pattern;
        Direction direction;
        uint32_t parent_position;
        char expected_base;
        Interval parent_primary;
        Interval parent_companion;
        Interval child_primary;
        Interval child_companion;
        uint8_t expected_nonempty_branches;
        bool sentinel_adjusts_full;
    };
    const std::vector<Example> examples{
        {"AAC", "AC", Direction::Left, 1, 'A',
         {3, 6}, {15, 18}, {1, 3}, {16, 18}, 1, true},
        {"AA", "A", Direction::Right, 0, 'A',
         {1, 6}, {18, 23}, {1, 3}, {21, 23}, 2, false}
    };

    for (const Example& example : examples) {
        SA_Range parent = exact_range(
            index, example.current_pattern, example.direction);
        parent.set_position(example.parent_position);
        parent.set_valid(true);
        const SA_Range exact = example.direction == Direction::Left
            ? parent.extend_left(code(example.expected_base))
            : parent.extend_right(code(example.expected_base));
        if (!same_interval(parent.primary_interval(), example.parent_primary) ||
            !same_interval(parent.companion_interval(), example.parent_companion) ||
            !same_interval(exact.primary_interval(), example.child_primary) ||
            !same_interval(exact.companion_interval(), example.child_companion))
            throw std::runtime_error("threshold interval fixture changed");

        const Interval relevant = example.direction == Direction::Left
            ? parent.primary_interval() : parent.companion_interval();
        const bool contains_sentinel =
            relevant.l <= index.bwt()->primary &&
            index.bwt()->primary < relevant.r;
        const uint64_t sentinel_offset = contains_sentinel ? 1 : 0;
        if (contains_sentinel != example.sentinel_adjusts_full ||
            (exact.size() + sentinel_offset == parent.size()) !=
                example.sentinel_adjusts_full)
            throw std::runtime_error("threshold sentinel classification changed");

        const BranchSet branches = parent.branch(example.direction);
        uint8_t nonempty = 0;
        for (uint8_t c = 0; c < 4; ++c)
            nonempty += branches.has_branch(c) ? 1 : 0;
        if (nonempty != example.expected_nonempty_branches ||
            !branches.has_branch(code(example.expected_base)) ||
            !same_interval(branches[code(example.expected_base)].primary_interval(),
                           exact.primary_interval()) ||
            !same_interval(branches[code(example.expected_base)].companion_interval(),
                           exact.companion_interval()))
            throw std::runtime_error("threshold branch oracle differs");

        const char first = example.direction == Direction::Right
            ? example.query.front() : example.query.back();
        SA_Range seed = index.initial_range(code(first));
        seed.set_position(0);
        seed.set_valid(true);
        std::vector<SA_Range> hits, rejects;
        search_unidirectional(example.query, seed, 0, 0, hits, rejects,
                              example.direction);

        SA_Range expected_reject = parent;
        expected_reject.set_forbidden_base(example.expected_base);
        bool found_reject = false;
        for (const SA_Range& reject : rejects)
            if (same_search_result(reject, expected_reject))
                found_reject = true;
        if (!found_reject)
            throw std::runtime_error("threshold reject state was not emitted");

        // The raw condition emits a reject in both cases. In the sentinel
        // example this documents accepted redundant-reject behaviour: the
        // sentinel-adjusted cardinality and BWA branch oracle show no
        // alternative nucleotide child to preserve. This is not a false hit
        // and does not establish a search correctness bug. The ordinary
        // partial-extension example has another nonempty child, so its
        // threshold parent reject preserves a real alternative.
        const bool raw_condition = parent.size() != exact.size();
        const bool adjusted_condition =
            exact.size() + sentinel_offset != parent.size();
        if (!raw_condition ||
            adjusted_condition != !example.sentinel_adjusts_full)
            throw std::runtime_error("threshold cardinality comparison changed");
    }
}


uint64_t expected_count(
    const BwaFMDIndex& index,
    const std::string& query,
    unsigned mismatches,
    Direction direction,
    const std::vector<std::string>& patterns,
    size_t seed_length = 1)
{
    uint64_t count = 0;

    assert(seed_length > 0 && seed_length <= query.size());
    const size_t seed_begin = direction == Direction::Right
        ? 0
        : query.size() - seed_length;

    for (const auto& candidate : patterns) {
        if (candidate.compare(seed_begin, seed_length,
                              query, seed_begin, seed_length) != 0)
            continue;

        if (hamming_distance(query, candidate) == mismatches)
            count += exact_range(index, candidate, direction).size();
    }

    return count;
}


uint64_t observed_count(
    const BwaFMDIndex& index,
    const std::string& query,
    unsigned mismatches,
    Direction direction)
{
    // The seed represents the first character consumed in the selected
    // search direction. position() is the zero-based consumed position.
    const char first = direction == Direction::Right
        ? query.front()
        : query.back();

    SA_Range seed = index.initial_range(code(first));
    seed.set_position(0);
    seed.set_valid(true);

    std::vector<SA_Range> hits;
    std::vector<SA_Range> rejects;
    search_unidirectional(
        query, seed, static_cast<int>(mismatches), 0,
        hits, rejects, direction);

    uint64_t count = 0;
    for (const auto& hit : hits) {
        assert(!hit.empty());
        assert(hit.mismatch_count() == mismatches);
        count += hit.size();
    }
    return count;
}

uint64_t observed_count_from_exact_seed(
    const BwaFMDIndex& index,
    const std::string& query,
    unsigned mismatches,
    Direction direction,
    size_t seed_length)
{
    assert(seed_length > 0 && seed_length < query.size());

    const size_t seed_begin = direction == Direction::Right
        ? 0
        : query.size() - seed_length;
    const std::string seed_read = query.substr(seed_begin, seed_length);
    SA_Range seed = search_exact_bidirectional(
        index, seed_read, direction);

    if (!seed.is_bidirectional() ||
        seed.position() != seed_length - 1) {
        std::cerr << "Exact bidirectional seed has invalid search state\n";
        std::abort();
    }

    std::vector<SA_Range> hits;
    std::vector<SA_Range> rejects;
    search_unidirectional(
        query, seed, static_cast<int>(mismatches), 0,
        hits, rejects, direction);

    uint64_t count = 0;
    for (const auto& hit : hits) {
        assert(!hit.empty());
        assert(hit.mismatch_count() == mismatches);
        count += hit.size();
    }
    return count;
}

void test_query(
    const BwaFMDIndex& index,
    const std::string& query,
    const std::vector<std::string>& patterns)
{
    for (Direction direction : {Direction::Right, Direction::Left}) {
        for (unsigned mismatches = 0; mismatches <= 1; ++mismatches) {
            const uint64_t expected = expected_count(
                index, query, mismatches, direction, patterns);
            const uint64_t observed = observed_count(
                index, query, mismatches, direction);

            if (expected != observed) {
                std::cerr
                    << "Mismatch-search count failure\n"
                    << "  query: " << query << '\n'
                    << "  direction: "
                    << (direction == Direction::Right ? "right" : "left")
                    << '\n'
                    << "  mismatches: " << mismatches << '\n'
                    << "  expected total interval size: " << expected << '\n'
                    << "  observed total interval size: " << observed << '\n';
                std::abort();
            }

            // A multi-base seed must remain a paired FMD range when passed
            // to mismatch search.  This covers exact continuation and the
            // one-mismatch branching path from an exact bidirectional seed.
            if (query.size() > 2) {
                constexpr size_t seed_length = 2;
                const uint64_t seeded_expected = expected_count(
                    index, query, mismatches, direction, patterns, seed_length);
                const uint64_t seeded_observed = observed_count_from_exact_seed(
                    index, query, mismatches, direction, seed_length);

                if (seeded_expected != seeded_observed) {
                    std::cerr
                        << "Exact-seed mismatch-search count failure\n"
                        << "  query: " << query << '\n'
                        << "  direction: "
                        << (direction == Direction::Right ? "right" : "left")
                        << '\n'
                        << "  mismatches: " << mismatches << '\n'
                        << "  expected total interval size: " << seeded_expected << '\n'
                        << "  observed total interval size: " << seeded_observed << '\n';
                    std::abort();
                }
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <index.bwt> <index.sa> <reference.fa>\n";
        return 2;
    }

    const std::string reference = read_single_fasta(argv[3]);
    const std::string indexed_text = reference + reverse_complement(reference);

    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) {
        std::cerr << "Failed to restore BWT: " << argv[1] << '\n';
        return 1;
    }

    bwt_restore_sa(argv[2], bwt);
    if (!bwt->sa) {
        std::cerr << "Failed to restore SA: " << argv[2] << '\n';
        bwt_destroy(bwt);
        return 1;
    }

    BwaFMDIndex index(bwt);
    if (bwt->seq_len != indexed_text.size()) {
        std::cerr << "Reference FASTA length does not match the indexed text\n";
        bwt_destroy(bwt);
        return 1;
    }
    size_t queries_tested = 0;
    std::vector<std::string> complete_result_queries;

    // Length 1 is deliberately omitted: the API seed is already one
    // consumed character, so there is no position at which to introduce
    // a mismatch for a one-character query.
    for (size_t length = 2; length <= 4; ++length) {
        std::vector<std::string> patterns;
        std::string current;
        generate_patterns(length, current, patterns);

        for (const auto& query : patterns) {
            test_query(index, query, patterns);
            complete_result_queries.push_back(query);
            ++queries_tested;
        }

        std::cout << "Mismatch-search length " << length << ": passed\n";
    }

    test_complete_ordered_results(
        index, complete_result_queries, indexed_text);
    test_over_budget_seed_rejections(index);
    test_threshold_cardinality_examples(index);

    bwt_destroy(bwt);
    std::cout << "Mismatch-search tests passed; queries tested: "
              << queries_tested << '\n';
    return 0;
}
