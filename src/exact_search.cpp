#include "exact_search.hpp"

#include <cstdint>
#include <stdexcept>

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
            "search_exact_unidirectional: read must contain only A/C/G/T");
    }
}

SA_Range empty_bidirectional_range(
    const BwaFMDIndex& index,
    uint32_t position)
{
    SA_Range result = SA_Range::bidirectional(&index, {}, {});
    result.set_position(position);
    result.set_valid(false);
    result.set_unidirectional(true);
    return result;
}

} // namespace

SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction)
{
    if (read.empty())
        throw std::invalid_argument(
            "search_exact_unidirectional: read must not be empty");

    const bool scan_right = direction == Direction::Right;
    const auto orientation = scan_right
        ? SA_Range::Orientation::Right
        : SA_Range::Orientation::Left;
    const size_t first_index = scan_right ? 0 : read.size() - 1;

    SA_Range range = index.initial_single_range(
        base_code(read[first_index]), orientation);
    range.set_position(0);
    range.set_valid(!range.empty());
    range.set_unidirectional(true);

    if (range.empty())
        return range;

    for (size_t consumed = 1; consumed < read.size(); ++consumed) {
        const size_t read_pos = scan_right
            ? consumed
            : read.size() - consumed - 1;
        const uint8_t c = base_code(read[read_pos]);

        range = scan_right
            ? index.extend_right(range, c)
            : index.extend_left(range, c);
        range.set_position(static_cast<uint32_t>(consumed));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);

        if (range.empty())
            return range;
    }

    return range;
}

SA_Range search_exact_bidirectional_seed(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction)
{
    if (read.empty())
        throw std::invalid_argument(
            "search_exact_bidirectional_seed: read must not be empty");

    const bool scan_right = direction == Direction::Right;
    const size_t first_index = scan_right ? 0 : read.size() - 1;

    SA_Range range = index.initial_range(base_code(read[first_index]));
    range.set_position(0);
    range.set_valid(!range.empty());
    range.set_unidirectional(true);

    if (range.empty())
        return range;

    for (size_t consumed = 1; consumed < read.size(); ++consumed) {
        const size_t read_pos = scan_right
            ? consumed
            : read.size() - consumed - 1;
        const uint8_t c = base_code(read[read_pos]);

        if (range.empty())
            return empty_bidirectional_range(
                index, static_cast<uint32_t>(consumed));

        if (range.size() == 1) {
            SA_Range next;
            const bool extended = scan_right
                ? index.extend_right_singleton(range, c, next)
                : index.extend_left_singleton(range, c, next);

            if (!extended)
                return empty_bidirectional_range(
                    index, static_cast<uint32_t>(consumed));

            range = next;
        } else {
            range = scan_right
                ? index.extend_right(range, c)
                : index.extend_left(range, c);
        }

        range.set_position(static_cast<uint32_t>(consumed));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
        if (range.empty()) return range;
    }

    return range;
}
