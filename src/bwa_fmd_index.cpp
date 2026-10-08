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
} // namespace

SA_Range BwaFMDIndex::initial_range(uint8_t c) const
{
    if (c > 3)
        throw std::invalid_argument("BwaFMDIndex::initial_range: base must be 0..3");

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
        throw std::invalid_argument("BwaFMDIndex extension: base must be 0..3");

    bwtintv_t ik = to_bwa(range);
    bwtintv_t ok[4]{};

    /*
     * BWA convention:
     *
     *   is_back = 1:
     *       extend the logical pattern to the LEFT.
     *
     *   is_back = 0:
     *       extend the logical pattern to the RIGHT.
     *
     * bwt_extend() also performs the required complement handling
     * for the paired FMD interval.
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

    const bwtint_t l = static_cast<bwtint_t>(p.l);

    // The BWT character at the unique primary row must be c.
    if (bwt_B0(bwt_, l) != c)
        return false;

    const bwtint_t new_l =
        bwt_->L2[c] + 1 + bwt_occ(bwt_, l, c);

    out = SA_Range::bidirectional(
        Interval{
            static_cast<uint64_t>(new_l),
            static_cast<uint64_t>(new_l + 1)
        },
        q);

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

    const bwtint_t l = static_cast<bwtint_t>(q.l);
    const uint8_t rc = static_cast<uint8_t>(3 - c);

    // RC(c) must occur immediately before RC(P) in the BWT sense.
    if (bwt_B0(bwt_, l) != rc)
        return false;

    const bwtint_t new_l =
        bwt_->L2[rc] + 1 + bwt_occ(bwt_, l, rc);

    out = SA_Range::bidirectional(
        p,
        Interval{
            static_cast<uint64_t>(new_l),
            static_cast<uint64_t>(new_l + 1)
        });

    return true;
}
