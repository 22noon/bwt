#pragma once

#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>

template <typename Mask>
class BitMaskIterator {
    static_assert(
        std::is_unsigned<Mask>::value,
        "BitMaskIterator requires an unsigned integer mask type");

public:
    using value_type = unsigned;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    constexpr BitMaskIterator() noexcept = default;

    explicit constexpr BitMaskIterator(Mask mask) noexcept
        : remaining_(mask)
    {}

    constexpr value_type operator*() const noexcept
    {
        /*
         * The iterator is never dereferenced when remaining_ == 0.
         */
        return static_cast<value_type>(
            __builtin_ctzll(
                static_cast<unsigned long long>(remaining_)));
    }

    constexpr BitMaskIterator& operator++() noexcept
    {
        remaining_ &= static_cast<Mask>(remaining_ - 1);
        return *this;
    }

    constexpr BitMaskIterator operator++(int) noexcept
    {
        BitMaskIterator previous = *this;
        ++(*this);
        return previous;
    }

    constexpr bool operator==(
        const BitMaskIterator& other) const noexcept
    {
        return remaining_ == other.remaining_;
    }

    constexpr bool operator!=(
        const BitMaskIterator& other) const noexcept
    {
        return remaining_ != other.remaining_;
    }

private:
    Mask remaining_ = 0;
};


template <typename Mask>
class BitMaskRange {
    static_assert(
        std::is_unsigned<Mask>::value,
        "BitMaskRange requires an unsigned integer mask type");

public:
    explicit constexpr BitMaskRange(Mask mask) noexcept
        : mask_(mask)
    {}

    constexpr BitMaskIterator<Mask> begin() const noexcept
    {
        return BitMaskIterator<Mask>(mask_);
    }

    constexpr BitMaskIterator<Mask> end() const noexcept
    {
        return BitMaskIterator<Mask>(0);
    }

private:
    Mask mask_;
};


template <typename Mask>
constexpr BitMaskRange<Mask> bitmask_range(Mask mask) noexcept
{
    return BitMaskRange<Mask>(mask);
}
