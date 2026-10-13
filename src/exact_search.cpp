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

SA_Range empty_range(const BwaFMDIndex& index, uint32_t position)
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

    const size_t first_index = direction == Direction::Right
        ? 0
        : read.size() - 1;

    SA_Range range = index.initial_range(base_code(read[first_index]));
    range.set_position(0);
    range.set_valid(!range.empty());
    range.set_unidirectional(true);

    if (range.empty())
        return range;

    for (size_t consumed = 1; consumed < read.size(); ++consumed) {
        const size_t read_pos = direction == Direction::Right
            ? consumed
            : read.size() - consumed - 1;
        const uint8_t c = base_code(read[read_pos]);

        if (range.empty())
            return empty_range(index, static_cast<uint32_t>(consumed));

        if (range.size() == 1) {
            SA_Range next;
            const bool extended = direction == Direction::Right
                ? index.extend_right_singleton(range, c, next)
                : index.extend_left_singleton(range, c, next);

            if (!extended)
                return empty_range(index, static_cast<uint32_t>(consumed));

            range = next;
            range.set_position(static_cast<uint32_t>(consumed));
            range.set_valid(true);
            range.set_unidirectional(true);
        } else {
            range = direction == Direction::Right
                ? index.extend_right(range, c)
                : index.extend_left(range, c);
            range.set_position(static_cast<uint32_t>(consumed));
            range.set_valid(!range.empty());
            range.set_unidirectional(true);

            if (range.empty())
                return range;
        }
    }

    return range;
}
