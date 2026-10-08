#include "sa_range.hpp"

SA_Range SA_Range::single(
    Interval interval,
    Orientation orientation)
{
    SA_Range r;
    r.mode_ = Mode::Single;
    r.orientation_ = orientation;
    r.interval_ = interval;
    return r;
}

SA_Range SA_Range::bidirectional(
    Interval primary,
    Interval companion)
{
    SA_Range r;
    r.mode_ = Mode::Bidirectional;
    r.primary_ = primary;
    r.companion_ = companion;
    return r;
}

bool SA_Range::is_single() const noexcept
{
    return mode_ == Mode::Single;
}

bool SA_Range::is_bidirectional() const noexcept
{
    return mode_ == Mode::Bidirectional;
}

SA_Range::Orientation SA_Range::orientation() const noexcept
{
    return orientation_;
}

Interval SA_Range::interval() const
{
    if (!is_single())
        throw std::logic_error(
            "SA_Range::interval(): range is bidirectional");
    return interval_;
}

Interval SA_Range::primary_interval() const
{
    if (!is_bidirectional())
        throw std::logic_error(
            "SA_Range::primary_interval(): range is single-index");
    return primary_;
}

Interval SA_Range::companion_interval() const
{
    if (!is_bidirectional())
        throw std::logic_error(
            "SA_Range::companion_interval(): range is single-index");
    return companion_;
}

uint64_t SA_Range::size() const
{
    return is_single() ? interval_.size() : primary_.size();
}

bool SA_Range::empty() const
{
    return size() == 0;
}

SA_Range SA_Range::to_single(Orientation orientation) const
{
    if (is_single())
        return *this;

    if (orientation == Orientation::Left)
        return SA_Range::single(primary_, Orientation::Left);

    return SA_Range::single(companion_, Orientation::Right);
}
