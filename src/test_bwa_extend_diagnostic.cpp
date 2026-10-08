#include "bwt.h"
#include "bwa_fmd_index.hpp"

#include <cassert>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

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

static bwtintv_t to_bwa(const SA_Range& r)
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
    const bwtintv_t in = to_bwa(before);
    const bwtintv_t expected = direct_extend(bwt, in, base, left);

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

static std::string key(const SA_Range& r)
{
    const auto p = r.primary_interval();
    const auto q = r.companion_interval();

    return std::to_string(p.l) + ":" + std::to_string(p.r) + ":" +
           std::to_string(q.l) + ":" + std::to_string(q.r);
}

static std::vector<SA_Range> collect_singletons(
    const bwt_t* bwt,
    unsigned max_depth)
{
    std::vector<SA_Range> current;
    std::vector<SA_Range> singletons;
    std::unordered_set<std::string> seen;

    for (uint8_t c = 0; c < 4; ++c) {
        bwtintv_t v{};
        bwt_set_intv(bwt, c, v);
        SA_Range r = SA_Range::bidirectional(
            Interval{v.x[0], v.x[0] + v.x[2]},
            Interval{v.x[1], v.x[1] + v.x[2]});

        current.push_back(r);
    }

    for (unsigned depth = 1; depth <= max_depth && !current.empty(); ++depth) {
        std::vector<SA_Range> next;

        for (const SA_Range& r : current) {
            if (r.size() == 1) {
                const std::string k = key(r);
                if (seen.insert(k).second)
                    singletons.push_back(r);
                continue;
            }

            const bwtintv_t in = to_bwa(r);
            bwtintv_t out[4]{};
            bwt_extend(bwt, &in, out, 1);

            for (uint8_t c = 0; c < 4; ++c) {
                if (out[c].x[2] == 0)
                    continue;

                next.push_back(SA_Range::bidirectional(
                    Interval{out[c].x[0], out[c].x[0] + out[c].x[2]},
                    Interval{out[c].x[1], out[c].x[1] + out[c].x[2]}));
            }
        }

        current.swap(next);
    }

    for (const SA_Range& r : current) {
        if (r.size() == 1) {
            const std::string k = key(r);
            if (seen.insert(k).second)
                singletons.push_back(r);
        }
    }

    return singletons;
}

static void check_singleton(
    const BwaFMDIndex& index,
    const bwt_t* bwt,
    const SA_Range& range)
{
    assert(range.size() == 1);

    const bwtintv_t in = to_bwa(range);

    bwtintv_t left_expected[4]{};
    bwtintv_t right_expected[4]{};
    bwt_extend(bwt, &in, left_expected, 1);
    bwt_extend(bwt, &in, right_expected, 0);

    for (uint8_t c = 0; c < 4; ++c) {
        const uint64_t sentinel = UINT64_C(0xdeadbeefcafebabe);
        SA_Range out = SA_Range::single(Interval{sentinel, sentinel + 1});

        const bool ok = index.extend_left_singleton(range, c, out);
        const bool expected_ok = left_expected[c].x[2] == 1;

        assert(ok == expected_ok);
        if (ok) {
            assert_same(out, left_expected[c]);
        } else {
            const auto p = out.interval();
            assert(p.l == sentinel);
            assert(p.r == sentinel + 1);
        }
    }

    for (uint8_t c = 0; c < 4; ++c) {
        const uint64_t sentinel = UINT64_C(0x123456789abcdef0);
        SA_Range out = SA_Range::single(Interval{sentinel, sentinel + 1});

        const bool ok = index.extend_right_singleton(range, c, out);
        const bool expected_ok = right_expected[c].x[2] == 1;

        assert(ok == expected_ok);
        if (ok) {
            assert_same(out, right_expected[c]);
        } else {
            const auto p = out.interval();
            assert(p.l == sentinel);
            assert(p.r == sentinel + 1);
        }
    }
}

static void check_singleton_fast_paths(
    const BwaFMDIndex& index,
    const bwt_t* bwt)
{
    const auto singletons = collect_singletons(bwt, 12);

    std::cout << "\nChecking singleton fast paths for "
              << singletons.size() << " singleton states...\n";

    for (const SA_Range& r : singletons)
        check_singleton(index, bwt, r);

    std::cout << "Singleton fast-path checks passed.\n";
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " <bwa-index.bwt>\n";
        return 2;
    }

    // This must be the .bwt file produced by `bwa index`, not the FASTA.
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

        check_singleton_fast_paths(index, bwt);
    }

    bwt_destroy(bwt);

    std::cout
        << "\nAll BWA/SA_Range diagnostic checks passed.\n";

    return 0;
}
