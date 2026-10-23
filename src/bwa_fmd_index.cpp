#include "bwa_fmd_index.hpp"

#include <stdexcept>
#include <vector>
#include <limits>
#include <utility>

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
    if (!range.is_single() && !range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::branch: unsupported range mode");

    const bool left = direction == Direction::Left;
    const int is_back = left ? 1 : 0;

    bwtintv_t ik{};

    if (range.is_bidirectional()) {
        ik = to_bwa(range);
    } else {
        const auto expected_orientation =
            left
                ? SA_Range::Orientation::Left
                : SA_Range::Orientation::Right;

        if (range.orientation() != expected_orientation)
            throw std::invalid_argument(
                "BwaFMDIndex::branch: single interval orientation "
                "does not match extension direction");

        const Interval interval = range.interval();

        // Populate only the active interval. The inactive start is
        // irrelevant to the active child intervals returned by bwt_extend.
        ik.x[!is_back] = static_cast<bwtint_t>(interval.l);
        ik.x[2] = static_cast<bwtint_t>(interval.size());
    }

    bwtintv_t ok[4]{};
    bwt_extend(bwt_, &ik, ok, is_back);

    BranchSet result;

    for (uint8_t c = 0; c < 4; ++c) {
        // BWA indexes right-extension branches by the complemented base.
        const uint8_t fmd_c =
            left ? c : static_cast<uint8_t>(3 - c);

        const auto& v = ok[fmd_c];

        if (range.is_bidirectional()) {
            result.branches_[c] =
                SA_Range::bidirectional(
                    this,
                    from_bwa_primary(v),
                    from_bwa_companion(v));
        } else {
            const bwtint_t start = v.x[!is_back];
            const bwtint_t size = v.x[2];

            result.branches_[c] =
                SA_Range::single(
                    this,
                    Interval{
                        static_cast<uint64_t>(start),
                        static_cast<uint64_t>(start + size)
                    },
                    range.orientation());
        }

        if (v.x[2] != 0)
            result.mask_ |= static_cast<uint8_t>(1u << c);
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
        this,
        from_bwa_primary(v),
        from_bwa_companion(v));
}

SA_Range BwaFMDIndex::initial_single_range(
    uint8_t c,
    SA_Range::Orientation orientation) const
{
    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::initial_single_range: base must be 0..3");

    /*
     * bwt_set_intv(bwt_, c, v) would initialise these two intervals:
     *
     *   primary:   [L2[c] + 1,     L2[c + 1] + 1)
     *   companion: [L2[3 - c] + 1, L2[4 - c] + 1)
     *
     * A right-oriented single range stores only the latter, which is the
     * interval for RC(c).  Constructing it directly avoids materialising the
     * companion/primary interval that a one-direction search will not use.
     */
    const uint8_t bwt_c = orientation == SA_Range::Orientation::Left
        ? c
        : static_cast<uint8_t>(3 - c);

    return SA_Range::single(
        this,
        Interval{
            static_cast<uint64_t>(bwt_->L2[bwt_c] + 1),
            static_cast<uint64_t>(bwt_->L2[bwt_c + 1] + 1)
        },
        orientation);
}


SA_Range BwaFMDIndex::extend_single_interval_one(
    const SA_Range& range,
    uint8_t bwt_c) const
{
    if (!range.is_single())
        throw std::logic_error(
            "BwaFMDIndex::extend_single_interval_one requires a single range");

    if (bwt_c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::extend_single_interval_one: base must be 0..3");

    const Interval interval = range.interval();

    // The experimental backend uses BatMIS-style exclusive-prefix ranks
    // over an in-memory BatMIS layout constructed from this BWA BWT.
    // BWA's native path below remains unchanged when it is selected.
    if (bat_rank_) {
        const bwtint_t start = bwt_->L2[bwt_c] + 1 +
            bat_rank_->occ(interval.l, bwt_c);
        if (batmis_singleton_fast_path_ && interval.size() == 1) {
            uint8_t row_symbol = 0;
            const bool has_symbol = bat_rank_->symbol_at(interval.l, row_symbol);
            const bool matches = has_symbol && row_symbol == bwt_c;
            const bwtint_t end = matches ? start + 1 : start;
            return SA_Range::single(
                this,
                Interval{static_cast<uint64_t>(start),
                         static_cast<uint64_t>(end)},
                range.orientation());
        }
        const bwtint_t end = bwt_->L2[bwt_c] + 1 +
            bat_rank_->occ(interval.r, bwt_c);
        if (end < start)
            return SA_Range::single(
                this, Interval{0, 0}, range.orientation());
        return SA_Range::single(
            this,
            Interval{static_cast<uint64_t>(start),
                     static_cast<uint64_t>(end)},
            range.orientation());
    }

    /*
     * Fast path for a singleton SA interval.
     *
     * The interval contains one BWT row. Its preceding BWT symbol
     * either matches bwt_c, producing one row after LF extension,
     * or does not match, producing an empty interval.
     *
     * The BWT representation omits the sentinel symbol from its
     * packed BWT string, so the primary row must be checked before
     * accessing bwt_B0().
     */
    if (interval.size() == 1) {
        const bwtint_t row =
            static_cast<bwtint_t>(interval.l);

        const bwtint_t rank_end =
            row == 0
                ? static_cast<bwtint_t>(-1)
                : row - 1;

        const bwtint_t rank =
            bwt_occ(bwt_, rank_end, bwt_c);

        const bwtint_t start =
            bwt_->L2[bwt_c] + 1 + rank;

        bool matches = false;

        if (row != bwt_->primary) {
            const bwtint_t packed_row =
                row - (row > bwt_->primary);

            matches =
                static_cast<uint8_t>(
                    bwt_B0(bwt_, packed_row)) == bwt_c;
        }

        const bwtint_t end =
            matches ? start + 1 : start;

        return SA_Range::single(
            this,
            Interval{
                static_cast<uint64_t>(start),
                static_cast<uint64_t>(end)
            },
            range.orientation());
    }

    /*
     * General single-interval extension.
     *
     * Keep the existing two-rank implementation for intervals
     * containing more than one row.
     */
    bwtint_t before = 0;
    bwtint_t through = 0;

    bwt_2occ(
        bwt_,
        static_cast<bwtint_t>(interval.l - 1),
        static_cast<bwtint_t>(interval.r - 1),
        bwt_c,
        &before,
        &through);

    const bwtint_t start =
        bwt_->L2[bwt_c] + 1 + before;

    return SA_Range::single(
        this,
        Interval{
            static_cast<uint64_t>(start),
            static_cast<uint64_t>(
                bwt_->L2[bwt_c] + 1 + through)
        },
        range.orientation());
}


SA_Range BwaFMDIndex::extend_bidirectional_one(
    const SA_Range& range,
    uint8_t c,
    bool left) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex extension requires a bidirectional range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex extension: base must be 0..3");

    const bwtintv_t ik = to_bwa(range);

    // BWA uses the primary interval for left extension and the
    // companion interval for right extension.
    const int is_back = left ? 1 : 0;
    const uint8_t fmd_c =
        left ? c : static_cast<uint8_t>(3 - c);

    bwtint_t tk[4], tl[4];

    bwt_2occ4(
        bwt_,
        ik.x[!is_back] - 1,
        ik.x[!is_back] - 1 + ik.x[2],
        tk,
        tl);

    // Calculate the requested child's interval on the side
    // being backward-extended.
    const bwtint_t child_size = tl[fmd_c] - tk[fmd_c];
    const bwtint_t child_start =
        bwt_->L2[fmd_c] + 1 + tk[fmd_c];

    // BWA's sentinel correction applies to the other interval.
    const bool contains_sentinel =
        ik.x[!is_back] <= bwt_->primary &&
        ik.x[!is_back] + ik.x[2] - 1 >= bwt_->primary;

    bwtint_t other_start =
        ik.x[is_back] + (contains_sentinel ? 1 : 0);

    // Position the selected child in the companion interval.
    // BWA orders the child intervals T, G, C, A, so the
    // sizes of all branches with larger alphabet codes precede it.
    for (int j = 3; j > fmd_c; --j)
        other_start += tl[j] - tk[j];

    const Interval active{
        static_cast<uint64_t>(child_start),
        static_cast<uint64_t>(child_start + child_size)
    };

    const Interval other{
        static_cast<uint64_t>(other_start),
        static_cast<uint64_t>(other_start + child_size)
    };

    return left
        ? SA_Range::bidirectional(this, active, other)
        : SA_Range::bidirectional(this, other, active);
}


SA_Range BwaFMDIndex::extend_left(
    const SA_Range& range,
    uint8_t c) const
{
    if (c > 3)
        throw std::invalid_argument("BwaFMDIndex::extend_left: base must be 0..3");

    if (range.empty())
        return range.is_single()
            ? SA_Range::single(this, range.interval(), range.orientation())
            : SA_Range::bidirectional(this, range.primary_interval(), range.companion_interval());

    if (range.size() == 1)
        return extend_cardinality_one(range, c, true);

    if (range.is_single()) {
        if (range.orientation() != SA_Range::Orientation::Left)
            throw std::logic_error(
                "BwaFMDIndex::extend_left requires a left-oriented single range");
        return extend_single_interval_one(range, c);
    }
    return extend_bidirectional_one(range, c, true);
}

SA_Range BwaFMDIndex::extend_right(
    const SA_Range& range,
    uint8_t c) const
{
    if (c > 3)
        throw std::invalid_argument("BwaFMDIndex::extend_right: base must be 0..3");

    if (range.empty())
        return range.is_single()
            ? SA_Range::single(this, range.interval(), range.orientation())
            : SA_Range::bidirectional(this, range.primary_interval(), range.companion_interval());

    if (range.size() == 1)
        return extend_cardinality_one(range, c, false);

    if (range.is_single()) {
        if (range.orientation() != SA_Range::Orientation::Right)
            throw std::logic_error(
                "BwaFMDIndex::extend_right requires a right-oriented single range");
        return extend_single_interval_one(range, static_cast<uint8_t>(3 - c));
    }
    return extend_bidirectional_one(range, c, false);
}

SA_Range BwaFMDIndex::extend_cardinality_one(
    const SA_Range& range,
    uint8_t c,
    bool left) const
{
    if (range.size() != 1 || range.empty())
        throw std::logic_error("extend_cardinality_one requires a non-empty cardinality-one range");

    if (range.is_single()) {
        if (range.orientation() != (left ? SA_Range::Orientation::Left : SA_Range::Orientation::Right))
            throw std::logic_error("single interval orientation does not match extension direction");
        return extend_single_interval_one(
            range, left ? c : static_cast<uint8_t>(3 - c));
    }

    SA_Range out;
    const bool ok = left
        ? extend_left_bidirectional_singleton(range, c, out)
        : extend_right_bidirectional_singleton(range, c, out);
    if (ok)
        return out;
    if (ok)
        return out;

    // No matching row: insertion-point coordinates are irrelevant.
    return SA_Range::bidirectional(
        this,
        Interval{0, 0},
        Interval{0, 0});
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
            this,
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
            this,
            from_bwa_primary(ok[fmd_c]),
            from_bwa_companion(ok[fmd_c]));
    }
}

bool BwaFMDIndex::extend_left_bidirectional_singleton(
    const SA_Range& range,
    uint8_t c,
    SA_Range& out) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_left_bidirectional_singleton requires a bidirectional range");

    if (range.size() != 1)
        throw std::logic_error(
            "BwaFMDIndex::extend_left_bidirectional_singleton requires a singleton range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::extend_left_bidirectional_singleton: base must be 0..3");

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
        this,
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

bool BwaFMDIndex::extend_right_bidirectional_singleton(
    const SA_Range& range,
    uint8_t c,
    SA_Range& out) const
{
    if (!range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::extend_right_bidirectional_singleton requires a bidirectional range");

    if (range.size() != 1)
        throw std::logic_error(
            "BwaFMDIndex::extend_right_bidirectional_singleton requires a singleton range");

    if (c > 3)
        throw std::invalid_argument(
            "BwaFMDIndex::extend_right_bidirectional_singleton: base must be 0..3");

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
        this,
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
    if (!range.is_single() && !range.is_bidirectional())
        throw std::logic_error(
            "BwaFMDIndex::locate: unsupported range mode");

    if (range.index() != this)
        throw std::invalid_argument(
            "BwaFMDIndex::locate: range belongs to a different index");

    if (!bwt_->sa || bwt_->sa_intv <= 0)
        throw std::logic_error(
            "BwaFMDIndex::locate requires the BWA suffix array "
            "to be restored");

    const Interval interval = range.is_single()
        ? range.interval()
        : range.primary_interval();

    const bwtint_t l = static_cast<bwtint_t>(interval.l);
    const bwtint_t r = static_cast<bwtint_t>(interval.r);

    std::vector<uint64_t> positions;
    positions.reserve(static_cast<size_t>(r - l));

    for (bwtint_t k = l; k < r; ++k)
        positions.push_back(
            static_cast<uint64_t>(bwt_sa(bwt_, k)));

    return positions;
}

std::vector<uint64_t> BwaFMDIndex::locate_forward(
    const SA_Range& range) const
{
    if (!range.is_single() && !range.is_bidirectional())
        throw std::logic_error("locate_forward: unsupported range mode");
    if (range.index() != this)
        throw std::invalid_argument("locate_forward: range belongs to a different index");
    if (!bwt_->sa || bwt_->sa_intv <= 0)
        throw std::logic_error("locate_forward requires the BWA sampled suffix array");

    const Interval interval = range.is_single()
        ? range.interval()
        : range.primary_interval();
    std::vector<uint64_t> positions;
    positions.reserve(static_cast<size_t>(interval.size()));
    for (bwtint_t row = static_cast<bwtint_t>(interval.l);
         row < static_cast<bwtint_t>(interval.r); ++row) {
        const bwtint_t coordinate = bwt_sa(bwt_, row);
        if (coordinate < bwt_->seq_len)
            positions.push_back(static_cast<uint64_t>(coordinate));
    }
    return positions;
}


#include <limits>
#include <utility>

std::vector<LocatedHit> BwaFMDIndex::locate(
    const SA_Range& range,
    uint64_t pattern_length) const
{
    if (!range.is_single() && !range.is_bidirectional())
        throw std::logic_error(
            "locate: unsupported range mode");

    if (range.index() != this)
        throw std::invalid_argument(
            "locate: range belongs to a different index");

    if (!bns_)
        throw std::logic_error(
            "locate requires BWA reference metadata");

    if (pattern_length == 0)
        throw std::invalid_argument(
            "locate: pattern length must be greater than zero");

    if (pattern_length >
        static_cast<uint64_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument(
            "locate: pattern length exceeds BWA metadata API limits");

    if (!bwt_->sa || bwt_->sa_intv <= 0)
        throw std::logic_error(
            "locate requires the BWA suffix array to be restored");

    const Interval interval = range.is_single()
        ? range.interval()
        : range.primary_interval();

    // A right-oriented single range stores RC(pattern), not pattern.
    // Coordinates are located normally, but its reported strand must
    // be flipped to describe the original query.
    const bool flip_strand =
        range.is_single() &&
        range.orientation() == SA_Range::Orientation::Right;

    const bwtint_t l = static_cast<bwtint_t>(interval.l);
    const bwtint_t r = static_cast<bwtint_t>(interval.r);

    std::vector<LocatedHit> hits;
    hits.reserve(static_cast<size_t>(r - l));

    for (bwtint_t k = l; k < r; ++k) {
        const int64_t pac =
            static_cast<int64_t>(bwt_sa(bwt_, k));

        int is_reverse = 0;
        int64_t pos = bns_depos(bns_, pac, &is_reverse);

        // For reverse-half hits, bns_depos() returns the mirrored
        // endpoint. Convert it to the leftmost reference coordinate.
        if (is_reverse) {
            if (pattern_length - 1 >
                static_cast<uint64_t>(pos))
                continue;

            pos -= static_cast<int64_t>(pattern_length - 1);
        }

        if (pos < 0 ||
            static_cast<uint64_t>(pos) >=
                static_cast<uint64_t>(bns_->l_pac))
            continue;

        if (pattern_length >
            static_cast<uint64_t>(bns_->l_pac - pos))
            continue;

        const int64_t end =
            pos + static_cast<int64_t>(pattern_length);

        const int ref_id = bns_intv2rid(bns_, pos, end);

        if (ref_id < 0 || ref_id >= bns_->n_seqs)
            continue;

        const bntann1_t& ref = bns_->anns[ref_id];
        const int64_t ref_end = ref.offset + ref.len;

        if (pos < ref.offset || end > ref_end)
            continue;

        int ambi_ref_id = -1;
        const int ambiguous = bns_cnt_ambi(
            bns_,
            pos,
            static_cast<int>(pattern_length),
            &ambi_ref_id);

        if (ambiguous != 0)
            continue;

        LocatedHit hit;
        hit.ref_id = static_cast<uint32_t>(ref_id);
        hit.ref_name = ref.name ? ref.name : "";
        hit.position =
            static_cast<uint64_t>(pos - ref.offset);

        const bool reported_reverse =
            static_cast<bool>(is_reverse) != flip_strand;

        hit.strand = reported_reverse
            ? Strand::Reverse
            : Strand::Forward;

        hits.push_back(std::move(hit));
    }

    return hits;
}

