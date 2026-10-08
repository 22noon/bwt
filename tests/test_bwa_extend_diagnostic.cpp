#include "bwt.h"
#include "bwa_fmd_index.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>

static void print_bwa(const char* s, const bwtintv_t& v)
{
    std::cout << s
              << ": x[0]=" << v.x[0]
              << " x[1]=" << v.x[1]
              << " size=" << v.x[2]
              << '\n';
}

static void print_range(const char* s, const SA_Range& r)
{
    const auto p = r.primary_interval();
    const auto q = r.companion_interval();

    std::cout << s
              << ": primary=[" << p.l << ',' << p.r << ')'
              << " companion=[" << q.l << ',' << q.r << ')'
              << " size=" << r.size()
              << '\n';
}

static void assert_same(
    const SA_Range& r,
    const bwtintv_t& v)
{
    const auto p = r.primary_interval();
    const auto q = r.companion_interval();

    assert(p.l == static_cast<uint64_t>(v.x[0]));
    assert(p.r == static_cast<uint64_t>(v.x[0] + v.x[2]));

    assert(q.l == static_cast<uint64_t>(v.x[1]));
    assert(q.r == static_cast<uint64_t>(v.x[1] + v.x[2]));

    assert(r.size() == static_cast<uint64_t>(v.x[2]));
}

static bwtintv_t direct_extend(
    const bwt_t* bwt,
    const bwtintv_t& in,
    uint8_t base,
    bool left)
{
    bwtintv_t out[4]{};

    // BWA:
    //   is_back = 1 -> logical left extension
    //   is_back = 0 -> logical right extension
    bwt_extend(bwt, &in, out, left ? 1 : 0);

    return out[base];
}

static void check_initial(
    const BwaFMDIndex& index,
    const bwt_t* bwt,
    uint8_t base,
    const char* name)
{
    bwtintv_t expected{};
    bwt_set_intv(bwt, base, expected);

    const SA_Range actual = index.initial_range(base);

    std::cout << "\nInitial " << name << '\n';
    print_bwa("BWA", expected);
    print_range("SA_Range", actual);

    assert_same(actual, expected);
}

static SA_Range check_extend(
    const BwaFMDIndex& index,
    const bwt_t* bwt,
    const SA_Range& before,
    uint8_t base,
    bool left,
    const char* operation,
    const char* pattern)
{
    // Convert the wrapper state back to exactly the bwtintv_t consumed
    // by bwt_extend(). This tests the SA_Range <-> BWA representation.
    bwtintv_t in{};

    const auto p = before.primary_interval();
    const auto q = before.companion_interval();

    in.x[0] = static_cast<bwtint_t>(p.l);
    in.x[1] = static_cast<bwtint_t>(q.l);
    in.x[2] = static_cast<bwtint_t>(p.size());
    in.info = 0;

    const bwtintv_t expected =
        direct_extend(bwt, in, base, left);

    const SA_Range actual =
        left
            ? index.extend_left(before, base)
            : index.extend_right(before, base);

    std::cout << '\n'
              << operation
              << "  pattern=" << pattern
              << '\n';

    print_bwa("BWA", expected);
    print_range("SA_Range", actual);

    assert_same(actual, expected);

    return actual;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " <bwa-index.bwt>\n";
        return 2;
    }

    // BWA's restore routine loads the same .bwt file produced by
    // `bwa index`.
    bwt_t* bwt = bwt_restore_bwt(argv[1]);

    if (bwt == nullptr) {
        std::cerr
            << "Could not load BWA BWT: "
            << argv[1] << '\n';
        return 1;
    }

    {
        BwaFMDIndex index(bwt);

        check_initial(index, bwt, 0, "A");
        check_initial(index, bwt, 1, "C");
        check_initial(index, bwt, 2, "G");
        check_initial(index, bwt, 3, "T");

        // Logical sequence:
        //
        //     A -> AC -> ACG -> TACG -> TACGT
        //
        SA_Range r = index.initial_range(0);

        r = check_extend(
            index, bwt, r, 1, false,
            "extend_right(C)", "AC");

        r = check_extend(
            index, bwt, r, 2, false,
            "extend_right(G)", "ACG");

        r = check_extend(
            index, bwt, r, 3, true,
            "extend_left(T)", "TACG");

        r = check_extend(
            index, bwt, r, 3, false,
            "extend_right(T)", "TACGT");
    }

    bwt_destroy(bwt);

    std::cout
        << "\nAll BWA/SA_Range diagnostic checks passed.\n";

    return 0;
}
