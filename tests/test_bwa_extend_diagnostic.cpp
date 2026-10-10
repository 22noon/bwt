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

static void assert_same(const SA_Range& r, const bwtintv_t& v)
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
    bwtintv_t in{};

    const auto p = before.primary_interval();
    const auto q = before.companion_interval();

    in.x[0] = static_cast<bwtint_t>(p.l);
    in.x[1] = static_cast<bwtint_t>(q.l);
    in.x[2] = static_cast<bwtint_t>(p.size());
    in.info = 0;

    /*
     * The public API uses logical bases.  BWA's right-extension
     * result is indexed by the complementary base because the
     * companion interval represents the reverse-complement side.
     */
    const uint8_t bwa_base =
        left ? base : static_cast<uint8_t>(3 - base);

    const bwtintv_t expected =
        direct_extend(bwt, in, bwa_base, left);

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

static std::string interval_key(const SA_Range& r)
{
    const auto p = r.primary_interval();
    const auto q = r.companion_interval();

    return std::to_string(p.l) + ":" +
           std::to_string(q.l) + ":" +
           std::to_string(r.size());
}

/*
 * Generate singleton states using the real bwt_extend() routine.
 * This avoids assuming that a particular hand-written sequence
 * happens to produce singleton intervals.
 */
static std::vector<SA_Range> collect_singletons(
    const BwaFMDIndex& index,
    const bwt_t* bwt,
    unsigned max_depth)
{
    std::vector<SA_Range> result;
    std::unordered_set<std::string> seen;

    std::vector<SA_Range> frontier;

    for (uint8_t c = 0; c < 4; ++c)
        frontier.push_back(index.initial_range(c));

    for (unsigned depth = 1; depth <= max_depth; ++depth) {
        std::vector<SA_Range> next;

        for (const SA_Range& r : frontier) {
            if (r.size() == 1) {
                const std::string k = interval_key(r);
                if (seen.insert(k).second)
                    result.push_back(r);
            }

            bwtintv_t in{};
            const auto p = r.primary_interval();
            const auto q = r.companion_interval();

            in.x[0] = static_cast<bwtint_t>(p.l);
            in.x[1] = static_cast<bwtint_t>(q.l);
            in.x[2] = static_cast<bwtint_t>(r.size());

            bwtintv_t ok[4]{};
            bwt_extend(bwt, &in, ok, 1);

            for (int c = 0; c < 4; ++c) {
                if (ok[c].x[2] == 0)
                    continue;

                SA_Range child = SA_Range::bidirectional(
                    &index,
                    Interval{
                        static_cast<uint64_t>(ok[c].x[0]),
                        static_cast<uint64_t>(ok[c].x[0] + ok[c].x[2])
                    },
                    Interval{
                        static_cast<uint64_t>(ok[c].x[1]),
                        static_cast<uint64_t>(ok[c].x[1] + ok[c].x[2])
                    });

                next.push_back(child);
            }
        }

        frontier.swap(next);
    }

    return result;
}

static void check_singleton(
    const BwaFMDIndex& index,
    const bwt_t* bwt,
    const SA_Range& range)
{
    for (uint8_t c = 0; c < 4; ++c) {
        bwtintv_t in{};
        const auto p = range.primary_interval();
        const auto q = range.companion_interval();

        in.x[0] = static_cast<bwtint_t>(p.l);
        in.x[1] = static_cast<bwtint_t>(q.l);
        in.x[2] = 1;
        in.info = 0;

        bwtintv_t left_ok[4]{};
        bwtintv_t right_ok[4]{};

        bwt_extend(bwt, &in, left_ok, 1);
        bwt_extend(bwt, &in, right_ok, 0);

        SA_Range left_out;
        SA_Range right_out;

        const bool left_success =
            index.extend_left_singleton(range, c, left_out);

        const bool right_success =
            index.extend_right_singleton(range, c, right_out);

        const bool expected_left = left_ok[c].x[2] == 1;

        const uint8_t right_bwa_c =
            static_cast<uint8_t>(3 - c);

        const bool expected_right =
            right_ok[right_bwa_c].x[2] == 1;

        if (left_success != expected_left) {
            std::cerr << "LEFT singleton mismatch: base=" << int(c)
                      << " range=" << interval_key(range)
                      << " expected_ok=" << expected_left
                      << " actual_ok=" << left_success << '\n';
            print_bwa("expected-left", left_ok[c]);
            return;
        }

        if (right_success != expected_right) {
            std::cerr << "RIGHT singleton mismatch: base=" << int(c)
                      << " range=" << interval_key(range)
                      << " expected_ok=" << expected_right
                      << " actual_ok=" << right_success << '\n';
            print_bwa("expected-right", right_ok[right_bwa_c]);
            return;
        }

        if (left_success)
            assert_same(left_out, left_ok[c]);

        if (right_success)
            assert_same(right_out, right_ok[right_bwa_c]);
    }
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " <bwa-index.bwt>\n";
        return 2;
    }

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

        const auto singleton_states =
            collect_singletons(index, bwt, 12);

        std::cout
            << "\nChecking singleton fast paths for "
            << singleton_states.size()
            << " singleton states...\n";

        for (const auto& s : singleton_states)
            check_singleton(index, bwt, s);

        std::cout << "Singleton fast-path checks passed.\n";
    }

    bwt_destroy(bwt);

    std::cout
        << "\nAll BWA/SA_Range diagnostic checks passed.\n";

    return 0;
}
