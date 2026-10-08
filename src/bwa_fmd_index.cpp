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

// Return the BWT character at a logical BWA BWT row.
// BWA's packed BWT omits the '$' row, which is represented by
// bwt->primary, so logical rows after primary need their packed
// index decremented by one.
bool bwt_char_at(const bwt_t* bwt, bwtint_t k, uint8_t& c)
{
    if (k == bwt->primary)
        return false; // '$' row has no nucleotide BWT character

    const bwtint_t packed_k = k - (k > bwt->primary);
    c = static_cast<uint8_t>(bwt_B0(bwt, packed_k));
    return true;
}

// LF mapping for a logical BWT row containing character c.
// bwt_occ() is inclusive, so Occ(c, k-1) gives the number of c's
// strictly before row k.
bwtint_t lf(const bwt_t* bwt, bwtint_t k, uint8_t c)
{
    const bwtint_t occ =
        bwt_occ(bwt, k == 0 ? static_cast<bwtint_t>(-1) : k - 1, c);

    return bwt->L2[c] + 1 + occ;
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
     * bwt_extend() updates both the primary and paired FMD intervals.
     * The paired interval is maintained by BWA's FMD representation.
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

    // For BWA's bwt_extend(..., is_back=1), a singleton can only
    // produce a non-empty child for the character stored in the BWT
    // at the unique primary row.
    uint8_t bwt_c = 0;
    if (!bwt_char_at(bwt_, l, bwt_c) || bwt_c != c)
        return false;

    // Primary interval: LF maps the unique row containing c.
    const bwtint_t new_primary = lf(bwt_, l, c);

    // Companion interval: reproduce bwt_extend()'s x[is_back] boundary
    // calculation. For a singleton only one child has size 1, namely c.
    // Therefore no other alphabet bucket contributes to the start of c.
    const bwtint_t crosses_dollar =
        (p.l <= bwt_->primary && p.l + 1 - 1 >= bwt_->primary) ? 1 : 0;
    const bwtint_t new_companion =
        static_cast<bwtint_t>(q.l) + crosses_dollar;

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
    const bwtint_t l = static_cast<bwtint_t>(q.l);

    // For BWA's bwt_extend(..., is_back=0), a singleton can only
    // produce a non-empty child for the character stored in the BWT
    // at the unique companion row.
    uint8_t bwt_c = 0;
    if (!bwt_char_at(bwt_, l, bwt_c) || bwt_c != c)
        return false;

    // Companion interval: LF maps the unique row containing c.
    const bwtint_t new_companion = lf(bwt_, l, c);

    // Primary interval: reproduce bwt_extend()'s x[is_back] boundary
    // calculation. For a singleton only one child has size 1, namely c.
    const bwtint_t crosses_dollar =
        (q.l <= bwt_->primary && q.l + 1 - 1 >= bwt_->primary) ? 1 : 0;
    const bwtint_t new_primary =
        static_cast<bwtint_t>(p.l) + crosses_dollar;

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
