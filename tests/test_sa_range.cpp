#include "sa_range.hpp"
#include "bwa_fmd_index.hpp"
#include "sa_branches.hpp"

#include <cassert>
#include <iostream>

int main()
{
    // Test construction & basic accessors
    Interval a{10, 15};
    auto s = SA_Range::single(nullptr, a, SA_Range::Orientation::Left);

    assert(s.is_single());
    assert(s.size() == 5);
    assert(s.index() == nullptr);

    Interval p{20, 25};
    Interval q{30, 35};
    auto b = SA_Range::bidirectional(nullptr, p, q);

    assert(b.is_bidirectional());
    assert(b.size() == 5);
    assert(b.index() == nullptr);

    // Test primary/companion size validation
    Interval q_bad{30, 36};
    bool threw_size = false;
    try {
        (void)SA_Range::bidirectional(nullptr, p, q_bad);
    } catch (const std::invalid_argument&) {
        threw_size = true;
    }
    assert(threw_size);

    // Test search state methods
    SA_Range empty_state = SA_Range::single(nullptr, {});
    assert(!empty_state.valid());

    SA_Range state = SA_Range::single(nullptr, a);
    assert(state.valid());
    state.set_valid(false);
    assert(!state.valid());
    state.set_valid(true);
    assert(state.valid());

    assert(state.position() == 0);
    state.set_position(5);
    assert(state.position() == 5);
    state.advance_position();
    assert(state.position() == 6);

    assert(state.mismatch_count() == 0);
    assert(state.mismatches().empty());
    state.add_mismatch(2, 'A');
    assert(state.mismatch_count() == 1);
    assert(state.mismatches().front().Pos == 2);
    assert(state.mismatches().front().Char == 'A');

    assert(!state.has_forbidden_base());
    state.set_forbidden_base('C');
    assert(state.has_forbidden_base());
    assert(state.forbidden_base() == 'C');

    assert(!state.is_unidirectional());
    state.set_unidirectional(true);
    assert(state.is_unidirectional());

    // Test to_single state propagation
    b.set_position(12);
    b.set_valid(true);
    b.add_mismatch(1, 'G');
    b.set_forbidden_base('T');
    b.set_unidirectional(true);

    auto l_prop = b.to_single(SA_Range::Orientation::Left);
    assert(l_prop.is_single());
    assert(l_prop.interval().l == 20);
    assert(l_prop.interval().r == 25);
    assert(l_prop.position() == 12);
    assert(l_prop.valid() == true);
    assert(l_prop.mismatch_count() == 1);
    assert(l_prop.mismatches()[0].Pos == 1);
    assert(l_prop.mismatches()[0].Char == 'G');
    assert(l_prop.forbidden_base() == 'T');
    assert(l_prop.is_unidirectional() == true);

    auto r_prop = b.to_single(SA_Range::Orientation::Right);
    assert(r_prop.is_single());
    assert(r_prop.interval().l == 30);
    assert(r_prop.interval().r == 35);
    assert(r_prop.position() == 12);
    assert(r_prop.valid() == true);
    assert(r_prop.mismatch_count() == 1);
    assert(r_prop.forbidden_base() == 'T');
    assert(r_prop.is_unidirectional() == true);

    std::cout << "SA_Range refactor & state propagation tests passed successfully\n";
}
