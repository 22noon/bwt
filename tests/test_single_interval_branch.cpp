
#include "bwa_fmd_index.hpp"
#include "exact_search.hpp"
#include "bwt.h"

#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

uint8_t code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: std::abort();
    }
}

char base(uint8_t c)
{
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    return bases[c];
}

[[noreturn]] void fail(
    const std::string& message,
    const std::string& pattern,
    Direction direction,
    uint8_t c)
{
    std::cerr
        << "FAIL: " << message
        << "\n  pattern: " << pattern
        << "\n  direction: "
        << (direction == Direction::Left ? "left" : "right")
        << "\n  branch: " << base(c) << '\n';
    std::abort();
}

bool same_nonempty_interval(
    const Interval& a,
    const Interval& b)
{
    if (a.size() != b.size())
        return false;

    // Empty intervals may have different insertion boundaries.
    return a.empty() || (a.l == b.l && a.r == b.r);
}

void check_single_branch(
    const BwaFMDIndex& index,
    const std::string& pattern,
    Direction direction)
{
    const SA_Range single = search_exact_unidirectional(
        index, pattern, direction);

    const SA_Range paired = search_exact_bidirectional(
        index, pattern, direction);

    if (!single.is_single())
        throw std::runtime_error("Exact search did not return a single interval");

    if (!paired.is_bidirectional())
        throw std::runtime_error("Exact search did not return a bidirectional interval");

    const auto expected_orientation =
        direction == Direction::Left
            ? SA_Range::Orientation::Left
            : SA_Range::Orientation::Right;

    if (single.orientation() != expected_orientation) {
        std::cerr << "Wrong orientation for pattern " << pattern << '\n';
        std::abort();
    }

    const BranchSet actual = index.branch(single, direction);
    const BranchSet paired_branches = index.branch(paired, direction);

    const Interval paired_active =
        direction == Direction::Left
            ? paired.primary_interval()
            : paired.companion_interval();

    if (single.size() != paired.size()) {
        std::cerr << "Parent interval sizes differ for " << pattern << '\n';
        std::abort();
    }

    for (uint8_t c = 0; c < 4; ++c) {
        const SA_Range& child = actual[c];

        // The branch must retain the single-interval representation.
        if (!child.is_single())
            fail("Branch changed interval mode", pattern, direction, c);

        if (child.orientation() != expected_orientation)
            fail("Branch changed orientation", pattern, direction, c);

        // Compare with the individual extension implementation.
        const SA_Range individual =
            direction == Direction::Left
                ? index.extend_left(single, c)
                : index.extend_right(single, c);

        if (!individual.is_single() ||
            individual.size() != child.size() ||
            !same_nonempty_interval(
                individual.interval(), child.interval())) {
            fail("Four-way branch differs from individual extension",
                 pattern, direction, c);
        }

        // Compare the single interval with the active side of the
        // equivalent bidirectional branch.
        const SA_Range& paired_child = paired_branches[c];

        const Interval paired_child_active =
            direction == Direction::Left
                ? paired_child.primary_interval()
                : paired_child.companion_interval();

        if (child.size() != paired_child.size() ||
            !same_nonempty_interval(
                child.interval(), paired_child_active)) {
            fail("Single branch differs from bidirectional branch",
                 pattern, direction, c);
        }

        const bool present = !child.empty();

        if (actual.has_branch(c) != present)
            fail("Incorrect branch mask", pattern, direction, c);
    }
}

void enumerate(
    const BwaFMDIndex& index,
    size_t target_length,
    std::string& pattern,
    uint64_t& tested)
{
    if (pattern.size() == target_length) {
        check_single_branch(index, pattern, Direction::Left);
        check_single_branch(index, pattern, Direction::Right);
        ++tested;
        return;
    }

    for (uint8_t c = 0; c < 4; ++c) {
        pattern.push_back(base(c));
        enumerate(index, target_length, pattern, tested);
        pattern.pop_back();
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr
            << "Usage: " << argv[0]
            << " <index.bwt> <index.sa>\n";
        return 2;
    }

    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) {
        std::cerr << "Could not restore BWT: " << argv[1] << '\n';
        return 1;
    }

    bwt_restore_sa(argv[2], bwt);
    if (!bwt->sa) {
        std::cerr << "Could not restore suffix array: " << argv[2] << '\n';
        bwt_destroy(bwt);
        return 1;
    }

    try {
        BwaFMDIndex index(bwt);
        uint64_t tested = 0;

        for (size_t length = 1; length <= 6; ++length) {
            std::string pattern;
            pattern.reserve(length);

            enumerate(index, length, pattern, tested);

            std::cout
                << "Single-interval branching length "
                << length << ": passed\n";
        }

        std::cout
            << "Single-interval branching passed for "
            << tested << " patterns in both directions.\n";
    } catch (...) {
        bwt_destroy(bwt);
        throw;
    }

    bwt_destroy(bwt);
    return 0;
}

