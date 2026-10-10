#include "bwa_fmd_index.hpp"

#include <stdexcept>
#include <vector>

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

BranchSet BwaFMDIndex::branch(
    const SA_Range& range,
    Direction direction) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::branch requires a bidirectional range");

    bwtintv_t ik = to_bwa(range);
    bwtintv_t ok[4]{};

    const bool left = direction == Direction::Left;

    bwt_extend(
        bwt_,
        &ik,
        ok,
        left ? 1 : 0);

    BranchSet result;

    for (uint8_t c = 0; c < 4; ++c) {
        /*
         * BWA's bwt_extend() uses the logical base directly for
         * left extension, but the complementary base for right
         * extension.
         */
        const uint8_t fmd_c =
            left ? c : static_cast<uint8_t>(3 - c);

        const auto& v = ok[fmd_c];

        result.branches_[c] =
            SA_Range::bidirectional(
                from_bwa_primary(v),
                from_bwa_companion(v));

        if (v.x[2] != 0)
            result.mask_ |=
                static_cast<uint8_t>(1u << c);
    }

    return result;
}

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
    return extend_all_one(range, static_cast<uint8_t>(3 - c), 0);
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

    /*
     * BWA's right-extension result is indexed by the base on the
     * companion/reverse-complement side.  The public API is expressed
     * in terms of the logical base appended to the right of P, so:
     *
     *     logical c -> BWA ok[complement(c)] = ok[3-c]
     */
    for (uint8_t c = 0; c < 4; ++c) {
        const uint8_t fmd_c = static_cast<uint8_t>(3 - c);

        out[c] = SA_Range::bidirectional(
            from_bwa_primary(ok[fmd_c]),
            from_bwa_companion(ok[fmd_c]));
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

    /*
     * For a singleton, the primary-side BWT interval contains exactly
     * one row.  The requested extension therefore either succeeds with
     * one row or fails.  Check that row directly.
     */
    uint8_t bwt_c = 0;
    if (primary == bwt_->primary)
        return false; // The BWT row at '$' has no nucleotide character.

    const bwtint_t packed_primary =
        primary - (primary > bwt_->primary);
    bwt_c = static_cast<uint8_t>(bwt_B0(bwt_, packed_primary));

    if (bwt_c != c)
        return false;

    /*
     * Left extension uses LF on the primary interval.  bwt_occ() is
     * inclusive, so the occurrence count must end at primary - 1.
     */
    const bwtint_t new_primary =
        bwt_->L2[c] + 1 + bwt_occ(
            bwt_,
            primary == 0 ? (bwtint_t)-1 : primary - 1,
            c);

    /*
     * The general bwt_extend() boundary is
     *
     *     old_start + contains_dollar + sum(size[j], j > c).
     *
     * Here the extension succeeded from a singleton, so the selected
     * child has size 1 and every other child has size 0.  The '$' row
     * was rejected above, so contains_dollar is also zero.
     */
    const bwtint_t new_companion = companion;

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

    /*
     * For right extension, the companion interval is the BWT side on
     * which the character test and LF calculation are performed.
     */
    if (companion == bwt_->primary)
        return false; // The BWT row at '$' has no nucleotide character.

    const bwtint_t packed_companion =
        companion - (companion > bwt_->primary);
    const uint8_t bwt_c =
        static_cast<uint8_t>(bwt_B0(bwt_, packed_companion));

    if (bwt_c != static_cast<uint8_t>(3 - c))
        return false;

    /*
     * bwt_occ() is inclusive, so the occurrence count must end at
     * companion - 1.
     */
    const uint8_t fmd_c = static_cast<uint8_t>(3 - c);

    const bwtint_t new_companion =
        bwt_->L2[fmd_c] + 1 + bwt_occ(
            bwt_,
            companion == 0 ? (bwtint_t)-1 : companion - 1,
            fmd_c);
    /*
     * The same singleton simplification applies to the paired primary
     * boundary.  The '$' row was rejected above, so its correction is
     * zero, and all cumulative child-size terms for bases greater than
     * c are zero.
     */
    const bwtint_t new_primary = primary;

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

std::vector<uint64_t> BwaFMDIndex::locate(
    const SA_Range& range) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::locate requires a bidirectional/FMD range");

    /*
     * bwt_restore_sa() must have been called before locate().
     */
    if (!bwt_->sa || bwt_->sa_intv <= 0)
        throw std::logic_error(
            "BwaFMDIndex::locate requires the BWA suffix array "
            "to be restored");

    const auto p = range.primary_interval();

    const bwtint_t l =
        static_cast<bwtint_t>(p.l);
    const bwtint_t r =
        static_cast<bwtint_t>(p.r);

    std::vector<uint64_t> positions;
    positions.reserve(static_cast<size_t>(r - l));

    for (bwtint_t k = l; k < r; ++k)
        positions.push_back(
            static_cast<uint64_t>(bwt_sa(bwt_, k)));

    return positions;
}
