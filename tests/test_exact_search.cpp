#include "bwa_fmd_index.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

uint8_t base_code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
        throw std::invalid_argument("invalid DNA base");
    }
}

char complement(char c)
{
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default:
        throw std::invalid_argument("invalid DNA base");
    }
}

std::string reverse_complement(const std::string& s)
{
    std::string rc;
    rc.reserve(s.size());

    for (auto it = s.rbegin(); it != s.rend(); ++it)
        rc.push_back(complement(*it));

    return rc;
}

/*
 * Exact search from left to right.
 *
 * Once the current interval becomes a singleton, use the singleton
 * fast path for all subsequent extensions. This deliberately exercises
 * both the general and singleton extension implementations.
 */
SA_Range exact_search(
    const BwaFMDIndex& index,
    const std::string& pattern)
{
    if (pattern.empty())
        throw std::invalid_argument("empty pattern");

    SA_Range range =
        index.initial_range(base_code(pattern[0]));

    for (size_t i = 1; i < pattern.size(); ++i) {
        const uint8_t c = base_code(pattern[i]);

        if (range.empty())
            return range;

        if (range.size() == 1) {
            SA_Range next;

            if (!index.extend_right_singleton(range, c, next))
                return SA_Range::bidirectional({}, {});

            range = next;
        } else {
            range = index.extend_right(range, c);
        }
    }

    return range;
}

/*
 * Naive exact search against the actual indexed text.
 */
std::vector<uint64_t> naive_locations(
    const std::string& text,
    const std::string& pattern)
{
    std::vector<uint64_t> result;

    if (pattern.empty())
        return result;

    for (size_t p = text.find(pattern);
         p != std::string::npos;
         p = text.find(pattern, p + 1))
    {
        result.push_back(static_cast<uint64_t>(p));
    }

    return result;
}

void check_query(
    const BwaFMDIndex& index,
    const std::string& indexed_text,
    const std::string& query)
{
    const SA_Range range = exact_search(index, query);

    const auto expected =
        naive_locations(indexed_text, query);

    std::vector<uint64_t> actual;

    if (!range.empty())
        actual = index.locate(range);

    /*
     * locate() returns positions in suffix-array order.
     * The naive search returns positions in text order.
     */
    auto expected_sorted = expected;
    auto actual_sorted = actual;

    std::sort(expected_sorted.begin(), expected_sorted.end());
    std::sort(actual_sorted.begin(), actual_sorted.end());

    if (actual_sorted != expected_sorted) {
        std::cerr << "Query failed: " << query << '\n';

        std::cerr << "  expected:";
        for (uint64_t p : expected_sorted)
            std::cerr << ' ' << p;

        std::cerr << "\n  actual:";
        for (uint64_t p : actual_sorted)
            std::cerr << ' ' << p;

        std::cerr << '\n';
        std::abort();
    }

    /*
     * Independently verify every reported position against the
     * actual indexed text. This protects against an internally
     * self-consistent but incorrect FM-index calculation.
     */
    for (uint64_t p : actual) {
        assert(p + query.size() <= indexed_text.size());

        assert(
            indexed_text.compare(
                static_cast<size_t>(p),
                query.size(),
                query) == 0);
    }

    std::cout
        << "  " << query
        << ": " << actual.size()
        << " exact occurrence(s)\n";
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr
            << "Usage: " << argv[0]
            << " <index.bwt> <index.sa>\n";

        return 2;
    }

    /*
     * Reference used by the test.
     *
     * BWA's FMD index is constructed over:
     *
     *     R + RC(R)
     *
     * Therefore locate() coordinates are interpreted against that
     * complete indexed text.
     */
    const std::string reference =
        "ACGTTGCAACG";

    const std::string indexed_text =
        reference + reverse_complement(reference);

    bwt_t* bwt =
        bwt_restore_bwt(argv[1]);

    if (!bwt) {
        std::cerr
            << "Failed to restore BWT: "
            << argv[1] << '\n';

        return 1;
    }

    bwt_restore_sa(argv[2], bwt);

    if (!bwt->sa) {
        std::cerr
            << "Failed to restore SA: "
            << argv[2] << '\n';

        bwt_destroy(bwt);
        return 1;
    }

    assert(
        bwt->seq_len ==
        indexed_text.size());

    BwaFMDIndex index(bwt);

    std::cout
        << "Reference:    "
        << reference << '\n';

    std::cout
        << "Indexed text: "
        << indexed_text << '\n';

    std::cout
        << "seq_len:      "
        << bwt->seq_len << '\n';

    /*
     * Multiple exact hits.
     */
    check_query(
        index,
        indexed_text,
        "ACG");

    /*
     * Becomes singleton at "TT", then exercises
     * extend_right_singleton() for the final G.
     */
    check_query(
        index,
        indexed_text,
        "TTG");

    /*
     * Becomes singleton at "GC", then exercises
     * extend_right_singleton() for the final A.
     */
    check_query(
        index,
        indexed_text,
        "GCA");

    /*
     * Singleton result without another extension.
     */
    check_query(
        index,
        indexed_text,
        "TT");

    /*
     * Absent sequence.
     */
    check_query(
        index,
        indexed_text,
        "AAA");

    /*
     * Full reference.
     */
    check_query(
        index,
        indexed_text,
        reference);

    /*
     * Reverse complement is also present because BWA indexes
     * R + RC(R).
     */
    check_query(
        index,
        indexed_text,
        reverse_complement(reference));

    bwt_destroy(bwt);

    std::cout
        << "All exact-search checks passed.\n";

    return 0;
}
