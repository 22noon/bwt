
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
            "Exact search: sequence must contain only A/C/G/T");
    }
}


SA_Range empty_bidirectional_range(
    const BwaFMDIndex& index,
    uint32_t position)
{
    SA_Range result =
        SA_Range::bidirectional(&index, {}, {});
    result.set_position(position);
    result.set_valid(false);
    result.set_unidirectional(true);
    return result;
}

void validate_seed(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    bool require_single,
    Direction direction)
{
    if (seed.index() != &index)
        throw std::invalid_argument(
            "Exact extension: seed belongs to a different index");

    if (seed.empty())
        return;

    if (require_single && !seed.is_single())
        throw std::invalid_argument(
            "Unidirectional extension requires a single SA interval");

    if (!require_single && !seed.is_bidirectional())
        throw std::invalid_argument(
            "Bidirectional extension requires a paired FMD interval");

    if (seed.is_single()) {
        const auto expected =
            direction == Direction::Right
                ? SA_Range::Orientation::Right
                : SA_Range::Orientation::Left;

        if (seed.orientation() != expected)
            throw std::invalid_argument(
                "Seed orientation does not match extension direction");
    }
}

/*
 * Attempt a paired FMD extension without calculating all four branches.
 *
 * This succeeds only if the requested nucleotide child contains every
 * parent row except a possible sentinel row:
 *
 *     child_size + sentinel_present == parent_size
 *
 * For left extension, the primary interval supplies the BWT ranks.
 * For right extension, the companion interval supplies the ranks and
 * the BWA alphabet symbol is complement(c).
 *
 * On success, derive the other interval using the FMD interval-boundary
 * formula. On failure, the caller must use the general extension.
 */
bool try_extend_single_path(
    const BwaFMDIndex& index,
    const SA_Range& parent,
    uint8_t c,
    Direction direction,
    SA_Range& out)
{
    if (!parent.is_bidirectional())
        throw std::logic_error(
            "Single-path extension requires a paired FMD range");

    if (parent.empty() || parent.size() <= 1)
        return false;

    if (c > 3)
        throw std::invalid_argument(
            "Single-path extension: base must be 0..3");

    const bwt_t* bwt = index.bwt();
    const Interval primary = parent.primary_interval();
    const Interval companion = parent.companion_interval();

    const bool extend_left = direction == Direction::Left;

    // The interval on which the requested BWT ranks are calculated.
    const Interval source = extend_left ? primary : companion;

    // Logical right extension is represented by the complement symbol
    // on BWA's companion/reverse-complement side.
    const uint8_t bwt_c = extend_left
        ? c
        : static_cast<uint8_t>(3 - c);

    const bwtint_t l = static_cast<bwtint_t>(source.l);
    const bwtint_t r = static_cast<bwtint_t>(source.r);

    bwtint_t before = 0;
    bwtint_t through = 0;

    bwt_2occ(
        bwt,
        l - 1,
        r - 1,
        bwt_c,
        &before,
        &through);

    const uint64_t child_size =
        static_cast<uint64_t>(through - before);

    // The sentinel row belongs to the source interval when its row
    // lies in the half-open interval [source.l, source.r).
    const bool sentinel_present =
        source.l <= static_cast<uint64_t>(bwt->primary) &&
        static_cast<uint64_t>(bwt->primary) < source.r;

    const uint64_t parent_size = parent.size();
    const uint64_t sentinel_adjustment =
        sentinel_present ? 1u : 0u;

    // Require a non-empty child and proof that no other nucleotide
    // branch can contain a row.
    if (child_size == 0 ||
        child_size + sentinel_adjustment != parent_size)
        return false;

    const bwtint_t child_start =
        bwt->L2[bwt_c] + 1 + before;

    Interval new_primary;
    Interval new_companion;

    if (extend_left) {
        new_primary = {
            static_cast<uint64_t>(child_start),
            static_cast<uint64_t>(child_start) + child_size
        };

        // BWA's left-extension boundary shifts by the sentinel row.
        const uint64_t other_start =
            companion.l + sentinel_adjustment;

        new_companion = {
            other_start,
            other_start + child_size
        };
    } else {
        new_companion = {
            static_cast<uint64_t>(child_start),
            static_cast<uint64_t>(child_start) + child_size
        };

        // BWA's right-extension boundary shifts by the sentinel row.
        const uint64_t other_start =
            primary.l + sentinel_adjustment;

        new_primary = {
            other_start,
            other_start + child_size
        };
    }

    out = SA_Range::bidirectional(
        &index, new_primary, new_companion);

    return true;
}

SA_Range extend_paired_one(
    const BwaFMDIndex& index,
    const SA_Range& parent,
    uint8_t c,
    Direction direction,
    bool use_single_path)
{
    if (parent.empty())
        return parent;

    // First use the existing cardinality-one fast path.
    if (parent.size() == 1) {
        SA_Range next;

        const bool ok = direction == Direction::Right
            ? index.extend_right_bidirectional_singleton(parent, c, next)
            : index.extend_left_bidirectional_singleton(parent, c, next);

        if (!ok)
            return empty_bidirectional_range(
                index, parent.position() + 1);

        return next;
    }

    // The additional single-path optimization is enabled only for
    // seeded bidirectional extension.
    if (use_single_path) {
        SA_Range next;

        if (try_extend_single_path(
                index, parent, c, direction, next))
            return next;
    }

    // General paired extension: calculates the FMD extension using
    // BWA's standard bwt_extend() implementation.
    return direction == Direction::Right
        ? index.extend_right(parent, c)
        : index.extend_left(parent, c);
}

} // namespace

namespace {

SA_Range search_exact_impl(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction,
    bool bidirectional,
    const char* empty_message)
{
    if (read.empty())
        throw std::invalid_argument(empty_message);

    const bool scan_right = direction == Direction::Right;

    const auto orientation = scan_right
        ? SA_Range::Orientation::Right
        : SA_Range::Orientation::Left;

    const size_t first = scan_right ? 0 : read.size() - 1;
    const uint8_t first_base = base_code(read[first]);

    SA_Range range = bidirectional
        ? index.initial_range(first_base)
        : index.initial_single_range(first_base, orientation);

    range.set_position(0);
    range.set_valid(!range.empty());
    range.set_unidirectional(true);

    for (size_t consumed = 1;
         consumed < read.size() && !range.empty();
         ++consumed) {

        const size_t pos = scan_right
            ? consumed
            : read.size() - consumed - 1;

        const uint8_t c = base_code(read[pos]);

        if (bidirectional) {
            // Preserve the paired singleton fast path.
            // Do not enable seeded single-path detection here.
            range = extend_paired_one(
                index, range, c, direction, false);
        } else {
            range = direction == Direction::Right
                ? index.extend_right(range, c)
                : index.extend_left(range, c);
        }

        range.set_position(
            static_cast<uint32_t>(consumed));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
    }

    return range;
}

} // namespace


SA_Range SA_Range::search_exact(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction,
    Mode mode)
{
    switch (mode) {
    case Mode::Single:
        return search_exact_impl(
            index,
            read,
            direction,
            false,
            "SA_Range::search_exact: read must not be empty");

    case Mode::Bidirectional:
        return search_exact_impl(
            index,
            read,
            direction,
            true,
            "SA_Range::search_exact: read must not be empty");
    }

    throw std::invalid_argument(
        "SA_Range::search_exact: invalid interval mode");
}

SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction)
{
    return SA_Range::search_exact(
        index, read, direction, SA_Range::Mode::Single);
}

SA_Range search_exact_bidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction)
{
    return SA_Range::search_exact(
        index, read, direction, SA_Range::Mode::Bidirectional);
}

SA_Range extend_exact_unidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction)
{
    validate_seed(index, seed, true, direction);

    SA_Range range = seed;

    if (extension.empty())
        return range;

    if (range.empty())
        return range;

    const bool scan_right = direction == Direction::Right;

    for (size_t consumed = 0;
         consumed < extension.size() && !range.empty();
         ++consumed) {
        const size_t pos = scan_right
            ? consumed
            : extension.size() - consumed - 1;

        range = direction == Direction::Right
            ? index.extend_right(range, base_code(extension[pos]))
            : index.extend_left(range, base_code(extension[pos]));

        range.set_position(seed.position() +
                           static_cast<uint32_t>(consumed + 1));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
    }

    return range;
}


SA_Range extend_exact_bidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction)
{
    validate_seed(index, seed, false, direction);

    SA_Range range = seed;

    if (extension.empty() || range.empty())
        return range;

    const bool scan_right = direction == Direction::Right;

    for (size_t consumed = 0;
         consumed < extension.size() && !range.empty();
         ++consumed) {
        const size_t pos = scan_right
            ? consumed
            : extension.size() - consumed - 1;

        range = extend_paired_one(
            index,
            range,
            base_code(extension[pos]),
            direction,
            true); // Seeded search enables single-path detection.

        range.set_position(seed.position() +
                           static_cast<uint32_t>(consumed + 1));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
    }

    return range;
}

