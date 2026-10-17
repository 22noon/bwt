#include "mismatch_search.hpp"

#include <cstdint>
#include <stdexcept>
#include <utility>

namespace {

uint8_t base_code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
        throw std::invalid_argument(
            "search_unidirectional: read must contain only A/C/G/T");
    }
}

char base_char(uint8_t c)
{
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    if (c >= 4) throw std::out_of_range("base_char: code must be 0..3");
    return bases[c];
}

int read_index(uint32_t position, int length, Direction direction)
{
    return direction == Direction::Left
        ? length - static_cast<int>(position) - 1
        : static_cast<int>(position);
}

SA_Range child_with_state(
    const SA_Range& child_interval,
    uint32_t next_position,
    char expected,
    char branch,
    Direction direction,
    bool capture_reference_bases,
    int read_length)
{
    SA_Range child = child_interval;
    child.set_position(next_position);
    child.set_valid(!child.empty());

    if (branch != expected) {
        if (capture_reference_bases) {
            const uint32_t read_pos_one_based =
                direction == Direction::Left
                    ? static_cast<uint32_t>(read_length - static_cast<int>(next_position))
                    : static_cast<uint32_t>(next_position + 1);
            child.add_mismatch(read_pos_one_based, branch);
        } else {
            child.add_mismatch(next_position, expected);
        }
    }
    return child;
}


SA_Range extend_expected_base(
    const SA_Range& range,
    uint8_t base,
    Direction direction)
{
    return direction == Direction::Left
        ? range.extend_left(base)
        : range.extend_right(base);
}

BranchSet branch_range(
    const SA_Range& range,
    Direction direction)
{
    return range.branch(direction);
}


} // namespace

static void extend_impl(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases)
{
    if (Count < 0)
        throw std::invalid_argument("Count must be non-negative");
    if (read.empty()) return;
    if (!seed.index())
    throw std::invalid_argument(
        "extend requires a valid index");
    if (seed.is_single()) {
        const auto expected_orientation =
            direction == Direction::Left
                ? SA_Range::Orientation::Left
                : SA_Range::Orientation::Right;

        if (seed.orientation() != expected_orientation)
            throw std::invalid_argument(
                "extend: single seed orientation does not match "
                "extension direction");
    } else if (!seed.is_bidirectional()) {
        throw std::invalid_argument(
            "extend: unsupported SA_Range mode");
    }
    if (seed.position() >= read.size())
        throw std::invalid_argument("seed position is outside the read");
    if (seed.position() == read.size() - 1) {
        if (seed.mismatch_count() == static_cast<uint32_t>(Count))
            Hit_List.push_back(seed);
        return;
    }
    if (!seed.valid() && seed.position() != 0) return;

    std::vector<SA_Range> stack;
    stack.reserve(4 * read.size());
    stack.push_back(seed);
    if (seed.position() == 0)
        stack.back().set_valid(true);
    int hits_found = 0;
    const int read_length = static_cast<int>(read.size());

    while (!stack.empty()) {
        SA_Range current = std::move(stack.back());
        stack.pop_back();
        if (!current.valid() || current.empty()) continue;

        const uint32_t next_position = current.position() + 1;
        if (next_position >= read.size()) continue;
        const char expected = read[static_cast<size_t>(
            read_index(next_position, read_length, direction))];
        const uint8_t expected_symbol = base_code(expected);
        const bool at_threshold =
            current.mismatch_count() == static_cast<uint32_t>(Count);

        if (at_threshold) {
            // Frozen design: exact extension first; branch only if its
            // cardinality differs from the parent.
            SA_Range exact = extend_expected_base( current, expected_symbol, direction);

            if (current.size() != exact.size()) {
                // This is a conservative threshold deferral, not an assertion
                // that another nucleotide child exists. Keep the parent
                // paired range and metadata, and forbid only the expected
                // base so a later pass can consider alternatives. The raw
                // cardinality comparison is intentional existing behaviour:
                // sentinel adjustment can make this state redundant, and
                // such entries remain part of the Rejects contract.
                SA_Range reject = current;
                reject.set_forbidden_base(expected);
                Rejects.push_back(std::move(reject));
            }

            if (exact.empty()) continue;
            exact.set_position(next_position);
            exact.set_valid(true);
            exact.set_forbidden_base('\0');
            if (next_position == read.size() - 1) {
                if (exact.mismatch_count() == static_cast<uint32_t>(Count)) {
                    Hit_List.push_back(std::move(exact));
                    if (MAXHITS > 0 && ++hits_found >= MAXHITS) return;
                }
            } else {
                stack.push_back(std::move(exact));
            }
            continue;
        }

        // Below threshold, explore all available nucleotide branches.
        const BranchSet branches = branch_range(current, direction);
        for (uint8_t c = 0; c < 4; ++c) {
            if (!branches.has_branch(c)) continue;
            if (current.has_forbidden_base() &&
                base_char(c) == current.forbidden_base())
                continue;

            SA_Range child = child_with_state(
                branches[c], next_position, expected,
                base_char(c), direction, capture_reference_bases,
                read_length);
            child.set_forbidden_base('\0');
            if (child.empty()) continue;

            if (child.mismatch_count() > static_cast<uint32_t>(Count)) {
                // Unlike a threshold deferral above, this is the actual
                // branch child that exceeded the mismatch budget. Preserve
                // its child intervals and mismatch metadata as constructed;
                // it is not a parent state with one base forbidden.
                Rejects.push_back(std::move(child));
                continue;
            }
            if (next_position == read.size() - 1) {
                if (child.mismatch_count() == static_cast<uint32_t>(Count)) {
                    Hit_List.push_back(std::move(child));
                    if (MAXHITS > 0 && ++hits_found >= MAXHITS) return;
                }
            } else {
                stack.push_back(std::move(child));
            }
        }
    }
}


void extend(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases)
{
    extend_impl(
        read,
        seed,
        Count,
        MAXHITS,
        Hit_List,
        Rejects,
        direction,
        capture_reference_bases);
}


void search_unidirectional(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases)
{
    extend(
        read,
        seed,
        Count,
        MAXHITS,
        Hit_List,
        Rejects,
        direction,
        capture_reference_bases);
}
