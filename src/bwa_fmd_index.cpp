#include "bwa_fmd_index.hpp"

#include <stdexcept>

namespace {

Interval from_bwa_primary(const bwtintv_t& v)
{
    return Interval{
        static_cast<uint64_t>(v.x[0]),
        static_cast<uint64_t>(v.x[0] + v.x[2])
    };
}

Interval from_bwa_companion(const bwtintv_t& v)
{
    return Interval{
        static_cast<uint64_t>(v.x[1]),
        static_cast<uint64_t>(v.x[1] + v.x[2])
    };
}

bwtintv_t to_bwa(const SA_Range& r)
{
    bwtintv_t v{};
    const auto p = r.primary_interval();
    const auto q = r.companion_interval();

    v.x[0] = static_cast<bwtint_t>(p.l);
    v.x[1] = static_cast<bwtint_t>(q.l);
    v.x[2] = static_cast<bwtint_t>(p.size());
    v.info = 0;
    return v;
}

/*
 * For a singleton FMD interval, bwt_extend() reduces to one BWT
 * character test plus the corresponding C-array/rank calculation.
 *
 * BWA's bwt_2occ4() uses inclusive Occ endpoints:
 *
 *     tk = Occ(*, l - 1)
 *     tl = Occ(*, l)
 *
 * Therefore tl[c] - tk[c] is exactly the number of occurrences of
 * character c in the singleton row.
 */
bool singleton_extension_sizes(
    const bwt_t* bwt,
    bwtint_t l,
    bwtint_t size[4])
{
    bwtint_t tk[4], tl[4];

    /*
     * l == 0 gives l - 1 == UINT64_MAX, which is BWA's
     * representation of -1 and is explicitly handled by bwt_occ4().
     */
    bwt_2occ4(
        bwt,
        l - 1,
        l,
        tk,
        tl);

    for (int c = 0; c < 4; ++c)
        size[c] = tl[c] - tk[c];

    return true;
}

/*
 * Reproduce the x[is_back] boundary construction in bwt_extend():
 *
 *   ok[3].x[is_back] = old_x[is_back] + contains_dollar;
 *   ok[2] = ok[3] + size[3];
 *   ok[1] = ok[2] + size[2];
 *   ok[0] = ok[1] + size[1];
 *
 * Thus the start of the child interval for base c is:
 *
 *   old_start + contains_dollar + sum(size[j], j > c)
 */
bwtint_t child_start(
    bwtint_t old_start,
    bwtint_t contains_dollar,
    const bwtint_t size[4],
    uint8_t c)
{
    bwtint_t start = old_start + contains_dollar;

    for (int j = 3; j > static_cast<int>(c); --j)
        start += size[j];

    return start;
}

} // namespace

SA_Range BwaFMDIndex::initial_range(uint8_t c) const
{
    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::initial_range: base must be 0..3");

    bwtintv_t v{};
    bwt_set_intv(bwt_, c, v);

    return SA_Range::bidirectional(
        from_bwa_primary(v),
        from_bwa_companion(v));
}

SA_Range BwaFMDIndex::extend_all_one(
    const SA_Range& range,
    uint8_t c,
    int is_back) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex extension requires a bidirectional/FMD range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex extension: base must be 0..3");

    bwtintv_t ik = to_bwa(range);
    bwtintv_t ok[4]{};

    /*
     * BWA convention:
     *
     *   is_back = 1:
     *       logical left extension.
     *
     *   is_back = 0:
     *       logical right extension.
     *
     * bwt_extend() handles the paired FMD interval internally.
     */
    bwt_extend(bwt_, &ik, ok, is_back);

    return SA_Range::bidirectional(
        from_bwa_primary(ok[c]),
        from_bwa_companion(ok[c]));
}

SA_Range BwaFMDIndex::extend_left(
    const SA_Range& range,
    uint8_t c) const
{
    return extend_all_one(range, c, 1);
}

SA_Range BwaFMDIndex::extend_right(
    const SA_Range& range,
    uint8_t c) const
{
    return extend_all_one(range, c, 0);
}

void BwaFMDIndex::extend_left_all(
    const SA_Range& range,
    SA_Range out[4]) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_left_all requires a bidirectional range");

    bwtintv_t ik = to_bwa(range);
    bwtintv_t ok[4]{};

    bwt_extend(bwt_, &ik, ok, 1);

    for (int c = 0; c < 4; ++c) {
        out[c] = SA_Range::bidirectional(
            from_bwa_primary(ok[c]),
            from_bwa_companion(ok[c]));
    }
}

void BwaFMDIndex::extend_right_all(
    const SA_Range& range,
    SA_Range out[4]) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_right_all requires a bidirectional range");

    bwtintv_t ik = to_bwa(range);
    bwtintv_t ok[4]{};

    bwt_extend(bwt_, &ik, ok, 0);

    for (int c = 0; c < 4; ++c) {
        out[c] = SA_Range::bidirectional(
            from_bwa_primary(ok[c]),
            from_bwa_companion(ok[c]));
    }
}

bool BwaFMDIndex::extend_left_singleton(
    const SA_Range& range,
    uint8_t c,
    SA_Range& out) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_left_singleton requires a bidirectional range");

    if (range.size() != 1)
        throw std::logic_error(
            "BwaFMDIndex::extend_left_singleton requires a singleton range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::extend_left_singleton: base must be 0..3");

    const auto p = range.primary_interval();
    const auto q = range.companion_interval();

    const bwtint_t primary = static_cast<bwtint_t>(p.l);
    const bwtint_t companion = static_cast<bwtint_t>(q.l);

    bwtint_t size[4];
    singleton_extension_sizes(bwt_, primary, size);

    /*
     * For left extension, bwt_extend(..., is_back=1) calculates
     * the new primary interval from Occ on the primary side.
     */
    if (size[c] != 1)
        return false;

    const bwtint_t new_primary =
        bwt_->L2[c] + 1 + bwt_occ(
            bwt_,
            primary == 0 ? (bwtint_t)-1 : primary - 1,
            c);

    /*
     * The companion start is the corresponding cumulative boundary.
     * The '$' row contributes one to the boundary if the primary
     * interval contains bwt->primary.
     */
    const bwtint_t contains_dollar =
        (primary <= bwt_->primary && bwt_->primary <= primary) ? 1 : 0;

    const bwtint_t new_companion =
        child_start(companion, contains_dollar, size, c);

    out = SA_Range::bidirectional(
        Interval{
            static_cast<uint64_t>(new_primary),
            static_cast<uint64_t>(new_primary + 1)
        },
        Interval{
            static_cast<uint64_t>(new_companion),
            static_cast<uint64_t>(new_companion + 1)
        });

    return true;
}

bool BwaFMDIndex::extend_right_singleton(
    const SA_Range& range,
    uint8_t c,
    SA_Range& out) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_right_singleton requires a bidirectional range");

    if (range.size() != 1)
        throw std::logic_error(
            "BwaFMDIndex::extend_right_singleton requires a singleton range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::extend_right_singleton: base must be 0..3");

    const auto p = range.primary_interval();
    const auto q = range.companion_interval();

    const bwtint_t primary = static_cast<bwtint_t>(p.l);
    const bwtint_t companion = static_cast<bwtint_t>(q.l);

    bwtint_t size[4];
    singleton_extension_sizes(bwt_, companion, size);

    /*
     * For right extension, bwt_extend(..., is_back=0) calculates
     * the new companion interval from Occ on the companion side.
     */
    if (size[c] != 1)
        return false;

    const bwtint_t new_companion =
        bwt_->L2[c] + 1 + bwt_occ(
            bwt_,
            companion == 0 ? (bwtint_t)-1 : companion - 1,
            c);

    /*
     * The primary start is the corresponding cumulative boundary.
     */
    const bwtint_t contains_dollar =
        (companion <= bwt_->primary && bwt_->primary <= companion) ? 1 : 0;

    const bwtint_t new_primary =
        child_start(primary, contains_dollar, size, c);

    out = SA_Range::bidirectional(
        Interval{
            static_cast<uint64_t>(new_primary),
            static_cast<uint64_t>(new_primary + 1)
        },
        Interval{
            static_cast<uint64_t>(new_companion),
            static_cast<uint64_t>(new_companion + 1)
        });

    return true;
}
