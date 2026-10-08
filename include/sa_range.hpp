#pragma once

#include <cstdint>
#include <stdexcept>

struct Interval {
    uint64_t l = 0;
    uint64_t r = 0;  // half-open [l,r)

    uint64_t size() const noexcept { return r - l; }
    bool empty() const noexcept { return l == r; }
};

class SA_Range {
public:
    enum class Mode { Single, Bidirectional };
    enum class Orientation { Left, Right };

private:
    Mode mode_ = Mode::Single;
    Orientation orientation_ = Orientation::Left;

    // Single-index/FMD representation. In this first BWA backend,
    // a single range is one BWT interval and does not retain a pattern.
    Interval interval_;

    // BWA/FMD representation:
    // primary_   = x[0]
    // companion_ = x[1]
    // size       = x[2]
    //
    // The pair represents the current logical pattern P and its
    // reverse-complement relationship in BWA's FMD representation.
    Interval primary_;
    Interval companion_;

public:
    SA_Range() = default;

    static SA_Range single(
        Interval interval,
        Orientation orientation = Orientation::Left);

    static SA_Range bidirectional(
        Interval primary,
        Interval companion);

    bool is_single() const noexcept;
    bool is_bidirectional() const noexcept;
    Orientation orientation() const noexcept;

    Interval interval() const;
    Interval primary_interval() const;
    Interval companion_interval() const;

    uint64_t size() const;
    bool empty() const;

    SA_Range to_single(Orientation orientation) const;
};
