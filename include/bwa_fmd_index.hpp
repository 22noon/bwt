#pragma once

extern "C" {
#include "bwt.h"
}

#include "sa_range.hpp"
#include <vector>

class BwaFMDIndex {
public:
    explicit BwaFMDIndex(const bwt_t* bwt)
        : bwt_(bwt)
    {
        if (!bwt_)
            throw std::invalid_argument("BwaFMDIndex: null bwt_t");
    }

    /*
     * Create the FMD state for a single base.
     *
     * BWA's bwt_set_intv() initializes both coupled intervals.
     */
    SA_Range initial_range(uint8_t c) const;

    /*
     * Logical extension operations.
     *
     * These always mean:
     *
     *     extend_left(c)  -> cP
     *     extend_right(c) -> Pc
     *
     * BWA's bwt_extend() handles the FMD pairing and complement
     * transformation internally.
     */
    SA_Range extend_left(
        const SA_Range& range,
        uint8_t c) const;

    SA_Range extend_right(
        const SA_Range& range,
        uint8_t c) const;

    /*
     * Four-way branching.
     *
     * Returned order is BWA's alphabet order:
     *     0=A, 1=C, 2=G, 3=T
     */
    void extend_left_all(
        const SA_Range& range,
        SA_Range out[4]) const;

    void extend_right_all(
        const SA_Range& range,
        SA_Range out[4]) const;

    /* Find branches at this point */
    BranchSet branch( const SA_Range& range, Direction direction) const;
    /*
     * Direct access for integration with the existing BWA API.
     */
    const bwt_t* bwt() const noexcept { return bwt_; }

    bool extend_left_singleton( const SA_Range& range, uint8_t c, SA_Range& out) const;
    bool extend_right_singleton( const SA_Range& range, uint8_t c, SA_Range& out) const;

    /*
     * Locate every suffix-array row in the primary interval.
     *
     * The returned positions are SA coordinates in BWA's indexed text
     * (reference + reverse complement), in the same order as the rows
     * in the primary interval.
     *
     * Requires the BWA sampled suffix array to have been restored.
     */
    std::vector<uint64_t> locate(
        const SA_Range& range) const;


private:
    const bwt_t* bwt_;

    SA_Range extend_all_one(
        const SA_Range& range,
        uint8_t c,
        int is_back) const;
};
