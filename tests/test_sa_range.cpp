#include "sa_range.hpp"

#include <cassert>
#include <iostream>

int main()
{
    Interval a{10, 15};
    auto s = SA_Range::single(a, SA_Range::Orientation::Left);

    assert(s.is_single());
    assert(s.size() == 5);

    Interval p{20, 25};
    Interval q{30, 35};
    auto b = SA_Range::bidirectional(p, q);

    assert(b.is_bidirectional());
    assert(b.size() == 5);

    auto l = b.to_single(SA_Range::Orientation::Left);
    assert(l.is_single());
    assert(l.interval().l == 20);
    assert(l.interval().r == 25);

    auto r = b.to_single(SA_Range::Orientation::Right);
    assert(r.is_single());
    assert(r.interval().l == 30);
    assert(r.interval().r == 35);

    std::cout << "SA_Range tests passed\n";
}
