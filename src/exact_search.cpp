
#include "exact_search.hpp"

#include <cstdint>
#include <stdexcept>
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
#include <vector>
#endif

namespace {

#if defined(FM_BENCH_SINGLETON_NOINLINE)
#define FM_SINGLETON_HELPER_ATTR __attribute__((noinline))
#else
#if defined(FM_BENCH_SINGLETON_ALWAYS_INLINE)
#define FM_SINGLETON_HELPER_ATTR inline __attribute__((always_inline))
#else
#define FM_SINGLETON_HELPER_ATTR inline
#endif
#endif

static FM_SINGLETON_HELPER_ATTR
SA_Range singleton_extension_body(
    const BwaFMDIndex& index,
    const SA_Range& range,
    uint8_t bwt_c)
{
    const bwt_t* bwt = index.bwt();
    const Interval interval = range.interval();
    const bwtint_t row = static_cast<bwtint_t>(interval.l);
    const bwtint_t rank_end = row == 0
        ? static_cast<bwtint_t>(-1)
        : row - 1;
    const bwtint_t rank = bwt_occ(bwt, rank_end, bwt_c);
    const bwtint_t start = bwt->L2[bwt_c] + 1 + rank;

    bool matches = false;
    if (row != bwt->primary) {
        const bwtint_t packed_row = row - (row > bwt->primary);
        matches = static_cast<uint8_t>(bwt_B0(bwt, packed_row)) == bwt_c;
    }

    const bwtint_t end = matches ? start + 1 : start;
    return SA_Range::single(
        &index,
        Interval{static_cast<uint64_t>(start), static_cast<uint64_t>(end)},
        range.orientation());
}

static inline Interval multirow_extension_body(
    const BwaFMDIndex& index,
    const Interval& interval,
    uint8_t bwt_c)
{
    const bwt_t* bwt = index.bwt();

    // bwt_2occ() takes inclusive BWT row endpoints. Represent the prefix
    // before row zero explicitly as BWA's unsigned -1 sentinel.
    const bwtint_t left_end = interval.l == 0
        ? static_cast<bwtint_t>(-1)
        : static_cast<bwtint_t>(interval.l - 1);
    const bwtint_t right_end = static_cast<bwtint_t>(interval.r - 1);

    bwtint_t before = 0;
    bwtint_t through = 0;
    bwt_2occ(bwt, left_end, right_end, bwt_c, &before, &through);

    const bwtint_t offset = bwt->L2[bwt_c] + 1;
    return Interval{
        static_cast<uint64_t>(offset + before),
        static_cast<uint64_t>(offset + through)
    };
}
#undef FM_SINGLETON_HELPER_ATTR

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

#if defined(FM_BENCH_SINGLETON_NOINLINE) || defined(FM_BENCH_SINGLETON_ALWAYS_INLINE)
SA_Range bench_singleton_extend_left(
    const BwaFMDIndex& index,
    const SA_Range& range,
    uint8_t c)
{
    if (index.rank_backend() == RankBackend::Bwa &&
        !range.empty() && range.is_single() && range.size() == 1 &&
        range.orientation() == SA_Range::Orientation::Left)
        return singleton_extension_body(index, range, c);
    return index.extend_left(range, c);
}
#endif

namespace {

SA_Range search_exact_impl(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction,
    bool bidirectional,
    const char* empty_message
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
    , std::vector<SA_Range>* trace = nullptr
#endif
    )
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
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
    if (trace)
        trace->push_back(range);
#endif

#if !defined(FM_BENCH_SINGLETON_STEP_BASELINE)
    // Production uses a dedicated run loop once the eligible BWA left-singleton
    // state is reached. The separately compiled baseline keeps the prior
    // per-step dispatch for controlled comparisons.
    size_t consumed = 1;
    while (consumed < read.size() && !range.empty()) {
        const size_t pos = scan_right
            ? consumed
            : read.size() - consumed - 1;
        uint8_t c = base_code(read[pos]);

        if (!bidirectional && direction == Direction::Left &&
            index.rank_backend() == RankBackend::Bwa && range.is_single() &&
            range.size() == 1) {
            // A singleton extension can only produce an empty interval or a
            // singleton interval. Once eligible, keep extending directly
            // without repeating generic mode/direction/backend/cardinality
            // dispatch; empty is the only loop condition needed.
            while (true) {
                range = singleton_extension_body(index, range, c);
                range.set_position(static_cast<uint32_t>(consumed));
                range.set_valid(!range.empty());
                range.set_unidirectional(true);
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
                if (trace)
                    trace->push_back(range);
#endif
                ++consumed;
                if (range.empty() || consumed >= read.size())
                    break;
                const size_t next_pos = scan_right
                    ? consumed
                    : read.size() - consumed - 1;
                c = base_code(read[next_pos]);
            }
            continue;
        }

#if !defined(FM_BENCH_DISABLE_MULTIROW_FAST_PATH)
        if (!bidirectional && direction == Direction::Left &&
            index.rank_backend() == RankBackend::Bwa && range.is_single() &&
            range.size() > 1) {
            // Stay in the native paired-rank path while the child remains
            // multi-row. A singleton is handed back to the singleton loop on
            // the next outer iteration; an empty child terminates the search.
            Interval interval = range.interval();
            const SA_Range::Orientation single_orientation = range.orientation();
            while (true) {
                interval = multirow_extension_body(index, interval, c);
                range = SA_Range::single(
                    &index, interval, single_orientation);
                range.set_position(static_cast<uint32_t>(consumed));
                range.set_valid(!range.empty());
                range.set_unidirectional(true);
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
                if (trace)
                    trace->push_back(range);
#endif
                ++consumed;
                if (interval.size() <= 1 ||
                    consumed >= read.size())
                    break;
                c = base_code(read[read.size() - consumed - 1]);
            }
            continue;
        }
#endif

        if (bidirectional) {
            range = extend_paired_one(
                index, range, c, direction, false);
        } else {
            range = direction == Direction::Right
                ? index.extend_right(range, c)
                : index.extend_left(range, c);
        }
        range.set_position(static_cast<uint32_t>(consumed));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
        if (trace)
            trace->push_back(range);
#endif
        ++consumed;
    }
#else
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
            if (direction == Direction::Left &&
                index.rank_backend() == RankBackend::Bwa && range.is_single() &&
                range.size() == 1)
                range = singleton_extension_body(index, range, c);
            else
                range = direction == Direction::Right
                    ? index.extend_right(range, c)
                    : index.extend_left(range, c);
        }

        range.set_position(
            static_cast<uint32_t>(consumed));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
        if (trace)
            trace->push_back(range);
#endif
    }
#endif

    return range;
}

} // namespace

#if defined(FM_BENCH_SINGLETON_RUN_LOOP_TRACE)
SA_Range bench_search_exact_unidirectional_trace(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction,
    std::vector<SA_Range>& trace)
{
    trace.clear();
    return search_exact_impl(
        index,
        read,
        direction,
        false,
        "bench_search_exact_unidirectional_trace: read must not be empty",
        &trace);
}

SA_Range bench_multirow_extend_left(
    const BwaFMDIndex& index,
    const SA_Range& range,
    uint8_t c)
{
    return SA_Range::single(
        &index,
        multirow_extension_body(index, range.interval(), c),
        range.orientation());
}
#endif


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

namespace {
uint16_t sixmer_key(const std::string& read, size_t start)
{
    uint16_t key = 0;
    for (size_t i = 0; i < 6; ++i)
        key = static_cast<uint16_t>((key << 2) | base_code(read[start + i]));
    return key;
}
}

ExactSearch6merTable::ExactSearch6merTable(const BwaFMDIndex& index)
    : index_(&index)
{
    static constexpr char bases[] = "ACGT";
    for (uint16_t key = 0; key < intervals_.size(); ++key) {
        std::string seed(6, 'A');
        uint16_t value = key;
        for (int i = 5; i >= 0; --i) {
            seed[static_cast<size_t>(i)] = bases[value & 3];
            value >>= 2;
        }
        intervals_[key] = search_exact_unidirectional(index, seed, Direction::Left).interval();
    }
}

SA_Range ExactSearch6merTable::search(const BwaFMDIndex& index, const std::string& read) const
{
    if (index_ != &index)
        throw std::invalid_argument("6-mer table belongs to a different BWT index");
    if (read.size() < 6)
        return search_exact_unidirectional(index, read, Direction::Left);

    const uint16_t key = sixmer_key(read, read.size() - 6);
    const Interval seed_interval = intervals_[key];
    // Preserve the baseline's precise empty interval bounds when the seed is
    // absent; the baseline can stop before reaching six consumed bases.
    if (seed_interval.empty())
        return search_exact_unidirectional(index, read, Direction::Left);

    SA_Range range = SA_Range::single(&index, seed_interval, SA_Range::Orientation::Left);
    range.set_position(5);
    range.set_valid(true);
    range.set_unidirectional(true);
    for (size_t pos = read.size() - 6; pos > 0 && !range.empty();) {
        --pos;
        range = index.extend_left(range, base_code(read[pos]));
        range.set_position(static_cast<uint32_t>(read.size() - pos - 1));
        range.set_valid(!range.empty());
        range.set_unidirectional(true);
    }
    return range;
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

