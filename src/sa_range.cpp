#include "sa_range.hpp"
#include "bwa_fmd_index.hpp"
#include "sa_branches.hpp"

#include <stdexcept>

void SA_Range::copy_search_state_to(SA_Range& destination) const
{
    destination.index_ = index_;
    destination.Pos = Pos;
    destination.Valid = Valid;
    destination.Mismatches = Mismatches;
    destination.MM = MM;
    destination.forbidden_char = forbidden_char;
    destination.unidirectional_search = unidirectional_search;
#ifdef DEBUG
    destination.Scanned_Sequence = Scanned_Sequence;
#endif
}

SA_Range SA_Range::single(
    const BwaFMDIndex* index,
    Interval interval,
    Orientation orientation)
{
    SA_Range r;
    r.index_ = index;
    r.mode_ = Mode::Single;
    r.orientation_ = orientation;
    r.interval_ = interval;
    r.Valid = !interval.empty();
    return r;
}

SA_Range SA_Range::bidirectional(
    const BwaFMDIndex* index,
    Interval primary,
    Interval companion)
{
    if (primary.size() != companion.size())
        throw std::invalid_argument(
            "SA_Range::bidirectional: primary and companion intervals must have equal size");

    SA_Range r;
    r.index_ = index;
    r.mode_ = Mode::Bidirectional;
    r.primary_ = primary;
    r.companion_ = companion;
    r.Valid = !primary.empty();
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

    SA_Range res;
    res.index_ = index_;
    res.mode_ = Mode::Single;
    res.orientation_ = orientation;
    if (orientation == Orientation::Left) {
        res.interval_ = primary_;
    } else {
        res.interval_ = companion_;
    }
    copy_search_state_to(res);
    return res;
}

SA_Range SA_Range::extend_left(uint8_t c) const
{
    if (!index_)
        throw std::logic_error("SA_Range::extend_left: null index pointer");
    SA_Range next = index_->extend_left(*this, c);
    copy_search_state_to(next);
    return next;
}

SA_Range SA_Range::extend_right(uint8_t c) const
{
    if (!index_)
        throw std::logic_error("SA_Range::extend_right: null index pointer");
    SA_Range next = index_->extend_right(*this, c);
    copy_search_state_to(next);
    return next;
}

BranchSet SA_Range::branch(Direction direction) const
{
    if (!index_)
        throw std::logic_error("SA_Range::branch: null index pointer");
    BranchSet bs = index_->branch(*this, direction);
    for (uint8_t c = 0; c < BranchSet::ALPHABET_SIZE; ++c) {
        copy_search_state_to(bs[c]);
    }
    return bs;
}

std::vector<LocatedHit> SA_Range::locate() const
{
    if (!index_)
        throw std::logic_error("SA_Range::locate: null index pointer");
    if (Pos == 0)
        throw std::logic_error("SA_Range::locate: Pos cannot be zero");
    return index_->locate(*this, static_cast<uint64_t>(Pos));
}
