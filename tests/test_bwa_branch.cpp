#include "bwa_fmd_index.hpp"
#include "bitmask_iterator.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
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
        assert(false);
        return 0;
    }
}

const char* base_name(uint8_t c)
{
    static const char* names[] = {"A", "C", "G", "T"};
    return names[c];
}


/*
 * Compare two SA_Range objects.
 *
 * For this test we deliberately compare the actual primary and
 * companion intervals rather than only the interval size.
 */
bool same_range(
    const SA_Range& a,
    const SA_Range& b)
{
    if (a.is_bidirectional() != b.is_bidirectional())
        return false;

    if (a.is_bidirectional()) {
        return
            a.primary_interval().l ==
                b.primary_interval().l &&
            a.primary_interval().r ==
                b.primary_interval().r &&
            a.companion_interval().l ==
                b.companion_interval().l &&
            a.companion_interval().r ==
                b.companion_interval().r;
    }

    return
        a.interval().l == b.interval().l &&
        a.interval().r == b.interval().r;
}


/*
 * Construct the FMD state for P by extending from left to right.
 */
SA_Range build_right(
    const BwaFMDIndex& index,
    const std::string& pattern)
{
    assert(!pattern.empty());

    SA_Range range =
        index.initial_range(code(pattern[0]));

    for (size_t i = 1; i < pattern.size(); ++i)
        range = index.extend_right(
            range,
            code(pattern[i]));

    return range;
}


/*
 * Construct the same logical pattern by extending from right to left.
 *
 * Starting with the last character, prepend the preceding characters.
 */
SA_Range build_left(
    const BwaFMDIndex& index,
    const std::string& pattern)
{
    assert(!pattern.empty());

    const size_t last = pattern.size() - 1;

    SA_Range range =
        index.initial_range(code(pattern[last]));

    for (size_t i = last; i-- > 0;)
        range = index.extend_left(
            range,
            code(pattern[i]));

    return range;
}


/*
 * Test the generic bit-mask iterator independently of the FM index.
 */
void test_bitmask_iterator()
{
    struct Case {
        uint8_t mask;
        std::vector<unsigned> expected;
    };

    const std::vector<Case> cases = {
        {0x00, {}},
        {0x01, {0}},
        {0x02, {1}},
        {0x04, {2}},
        {0x08, {3}},
        {0x03, {0, 1}},
        {0x05, {0, 2}},
        {0x0A, {1, 3}},
        {0x0B, {0, 1, 3}},
        {0x0F, {0, 1, 2, 3}},
        {0x55, {0, 2, 4, 6}},
        {0xAA, {1, 3, 5, 7}},
    };

    for (const auto& test : cases) {
        std::vector<unsigned> actual;

        for (unsigned bit : bitmask_range(test.mask))
            actual.push_back(bit);

        assert(actual == test.expected);
    }

    /*
     * Also check that incrementing the iterator removes exactly
     * the lowest set bit each time.
     */
    {
        auto it = BitMaskIterator<uint8_t>(0b10110100);
        auto end = BitMaskIterator<uint8_t>(0);

        const unsigned expected[] = {2, 4, 5, 7};
        size_t i = 0;

        while (it != end) {
            assert(i < 4);
            assert(*it == expected[i]);
            ++it;
            ++i;
        }

        assert(i == 4);
    }

    std::cout << "BitMaskIterator tests passed.\n";
}


/*
 * Verify one BranchSet against the corresponding individual
 * extension primitive.
 */
void check_branch_set(
    const BwaFMDIndex& index,
    const SA_Range& range,
    Direction direction,
    const std::string& pattern)
{
    const BranchSet branches =
        index.branch(range, direction);

    uint8_t expected_mask = 0;

    for (uint8_t c = 0; c < 4; ++c) {
        const SA_Range expected =
            direction == Direction::Left
                ? index.extend_left(range, c)
                : index.extend_right(range, c);

        const bool expected_present =
            !expected.empty();

        const bool actual_present =
            branches.has_branch(c);

        if (expected_present != actual_present) {
            std::cerr
                << "Branch mask mismatch\n"
                << "  pattern:   " << pattern << '\n'
                << "  direction: "
                << (direction == Direction::Left
                    ? "left"
                    : "right")
                << '\n'
                << "  base:      " << base_name(c) << '\n'
                << "  expected:  "
                << expected_present << '\n'
                << "  actual:    "
                << actual_present << '\n';

            std::abort();
        }

        if (expected_present)
            expected_mask |=
                static_cast<uint8_t>(1u << c);

        /*
         * The BranchSet must contain exactly the same SA range as
         * the ordinary extension primitive.
         *
         * We check this even for empty branches because branch()
         * deliberately constructs all four child ranges.
         */
        const SA_Range& actual =
            branches[c];

        if (!same_range(actual, expected)) {
            std::cerr
                << "Branch range mismatch\n"
                << "  pattern:   " << pattern << '\n'
                << "  direction: "
                << (direction == Direction::Left
                    ? "left"
                    : "right")
                << '\n'
                << "  base:      " << base_name(c) << '\n';

            std::abort();
        }
    }

    /*
     * Check the complete mask.
     */
    assert(branches.branch_mask() == expected_mask);

    /*
     * Check the iterator.
     *
     * It must return precisely the set bits in increasing order.
     */
    uint8_t iterator_mask = 0;

    unsigned previous = 0;
    bool first = true;

    for (unsigned c : branches.possible_branches()) {
        assert(c < 4);

        if (!first)
            assert(c > previous);

        first = false;
        previous = c;

        iterator_mask |=
            static_cast<uint8_t>(1u << c);

        /*
         * Every value produced by the iterator must correspond
         * to an actual branch.
         */
        assert(branches.has_branch(
            static_cast<uint8_t>(c)));
    }

    assert(iterator_mask == expected_mask);
}


/*
 * Generate all DNA strings of a given length.
 */
void generate_patterns(
    size_t length,
    std::string& pattern,
    std::vector<std::string>& patterns)
{
    if (pattern.size() == length) {
        patterns.push_back(pattern);
        return;
    }

    static const char bases[] = {'A', 'C', 'G', 'T'};

    for (char c : bases) {
        pattern.push_back(c);

        generate_patterns(
            length,
            pattern,
            patterns);

        pattern.pop_back();
    }
}


void test_branching(
    const BwaFMDIndex& index)
{
    size_t states_tested = 0;

    /*
     * Test every pattern from length 1 through 6.
     *
     * This gives:
     *
     *   4 + 16 + 64 + 256 + 1024 + 4096
     *   = 5460 patterns
     */
    for (size_t length = 1; length <= 6; ++length) {
        std::vector<std::string> patterns;
        std::string pattern;

        generate_patterns(
            length,
            pattern,
            patterns);

        for (const std::string& p : patterns) {
            /*
             * Skip empty FM states. There are no meaningful
             * branches from an empty search interval.
             */
            const SA_Range right_state =
                build_right(index, p);

            if (!right_state.empty()) {
                check_branch_set(
                    index,
                    right_state,
                    Direction::Right,
                    p);

                ++states_tested;
            }

            const SA_Range left_state =
                build_left(index, p);

            if (!left_state.empty()) {
                check_branch_set(
                    index,
                    left_state,
                    Direction::Left,
                    p);

                ++states_tested;
            }
        }

        std::cout
            << "Branching length "
            << length
            << ": passed\n";
    }

    std::cout
        << "Branching tests passed.\n"
        << "Non-empty states tested: "
        << states_tested
        << "\n";
}

} // namespace


int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr
            << "Usage: "
            << argv[0]
            << " <index.bwt> <index.sa>\n";

        return 2;
    }

    bwt_t* bwt =
        bwt_restore_bwt(argv[1]);

    if (!bwt) {
        std::cerr
            << "Failed to restore BWT: "
            << argv[1] << '\n';

        return 1;
    }

    bwt_restore_sa(argv[2], bwt);

    if (!bwt->sa) {
        std::cerr
            << "Failed to restore SA: "
            << argv[2] << '\n';

        bwt_destroy(bwt);
        return 1;
    }

    BwaFMDIndex index(bwt);

    test_bitmask_iterator();
    test_branching(index);

    bwt_destroy(bwt);

    std::cout
        << "\nAll branching tests passed.\n";

    return 0;
}
