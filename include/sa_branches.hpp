#pragma once

#include "bitmask_iterator.hpp"
#include "sa_range.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

class BwaFMDIndex;

class BranchSet {
public:
    static constexpr uint8_t ALPHABET_SIZE = 4;

    using BranchIterator = BitMaskIterator<uint8_t>;

    BranchSet() = default;

    uint8_t branch_mask() const noexcept
    {
        return mask_;
    }

    bool has_branch(uint8_t c) const
    {
        if (c >= ALPHABET_SIZE)
            throw std::out_of_range(
                "BranchSet: base must be 0..3");

        return (mask_ & (uint8_t{1} << c)) != 0;
    }

    const SA_Range& operator[](uint8_t c) const
    {
        if (c >= ALPHABET_SIZE)
            throw std::out_of_range(
                "BranchSet: base must be 0..3");

        return branches_[c];
    }

    SA_Range& operator[](uint8_t c)
    {
        if (c >= ALPHABET_SIZE)
            throw std::out_of_range(
                "BranchSet: base must be 0..3");

        return branches_[c];
    }

    BitMaskRange<uint8_t> possible_branches() const noexcept
    {
        return bitmask_range(mask_);
    }

private:
    friend class BwaFMDIndex;

    std::array<SA_Range, ALPHABET_SIZE> branches_{};
    uint8_t mask_ = 0;
};
