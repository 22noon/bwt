#include "bwa_fmd_index.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

constexpr uint8_t A = 0;
constexpr uint8_t C = 1;
constexpr uint8_t G = 2;
constexpr uint8_t T = 3;

uint8_t code(char c)
{
    switch (c) {
    case 'A': return A;
    case 'C': return C;
    case 'G': return G;
    case 'T': return T;
    default:
        throw std::invalid_argument("invalid DNA base");
    }
}

SA_Range search(
    const BwaFMDIndex& index,
    const std::string& pattern)
{
    if (pattern.empty())
        throw std::invalid_argument("empty pattern");

    SA_Range range =
        index.initial_range(code(pattern[0]));

    for (size_t i = 1; i < pattern.size(); ++i)
        range = index.extend_right(
            range,
            code(pattern[i]));

    return range;
}

struct ExpectedHit {
    uint32_t ref_id;
    std::string ref_name;
    uint64_t position;
    Strand strand;
};

bool same_hit(
    const LocatedHit& actual,
    const ExpectedHit& expected)
{
    return actual.ref_id == expected.ref_id &&
           actual.ref_name == expected.ref_name &&
           actual.position == expected.position &&
           actual.strand == expected.strand;
}

void assert_expected_hits(
    const std::vector<LocatedHit>& actual,
    std::vector<ExpectedHit> expected)
{
    assert(actual.size() == expected.size());

    /*
     * locate() returns hits in SA-row order.  That order is not part
     * of the metadata-aware locate() API, so compare as an unordered
     * set of (reference, position, strand) tuples.
     */
    std::vector<bool> matched(actual.size(), false);

    for (const auto& e : expected) {
        bool found = false;

        for (size_t i = 0; i < actual.size(); ++i) {
            if (!matched[i] && same_hit(actual[i], e)) {
                matched[i] = true;
                found = true;
                break;
            }
        }

        if (!found) {
            std::cerr
                << "Expected hit not found:"
                << " ref_id=" << e.ref_id
                << " ref_name=" << e.ref_name
                << " position=" << e.position
                << " strand="
                << (e.strand == Strand::Forward
                        ? "Forward"
                        : "Reverse")
                << '\n';

            std::cerr << "Actual hits:\n";

            for (const auto& hit : actual) {
                std::cerr
                    << "  ref_id=" << hit.ref_id
                    << " ref_name=" << hit.ref_name
                    << " position=" << hit.position
                    << " strand="
                    << (hit.strand == Strand::Forward
                            ? "Forward"
                            : "Reverse")
                    << '\n';
            }

            assert(false);
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " <bwa-index-prefix>\n";
        return 2;
    }

    const std::string prefix = argv[1];

    const std::string bwt_file = prefix + ".bwt";
    const std::string sa_file  = prefix + ".sa";

    /*
     * This test assumes the BWA index was made from:
     *
     *   >ref1
     *   ACGTTGCAACG
     *
     *   >ref2
     *   TTACCGGTAA
     *
     * Therefore:
     *
     *   ref1 length = 11
     *   ref2 length = 10
     *   l_pac       = 21
     */

    bwt_t* bwt = bwt_restore_bwt(bwt_file.c_str());

    if (!bwt) {
        std::cerr
            << "Failed to restore BWT: "
            << bwt_file << '\n';
        return 1;
    }

    bwt_restore_sa(sa_file.c_str(), bwt);

    if (!bwt->sa) {
        std::cerr
            << "Failed to restore SA: "
            << sa_file << '\n';

        bwt_destroy(bwt);
        return 1;
    }

    /*
     * Load BWA reference metadata.
     *
     * bns_restore() loads .ann, .amb and .pac.
     */
    bntseq_t* bns = bns_restore(prefix.c_str());

    if (!bns) {
        std::cerr
            << "Failed to restore BWA reference metadata for prefix: "
            << prefix << '\n';

        bwt_destroy(bwt);
        return 1;
    }

    assert(bns->n_seqs == 2);
    assert(bns->l_pac == 21);

    BwaFMDIndex index(bwt, bns);

    /*
     * ------------------------------------------------------------
     * 1. ACG
     * ------------------------------------------------------------
     *
     * ref1 = A C G T T G C A A C G
     *        ^^^
     *         0
     *
     * ACG also occurs at position 8.
     *
     * RC(ACG) = CGT, which occurs at position 1.
     *
     * Therefore:
     *
     *   (ref1, 0, Forward)
     *   (ref1, 1, Reverse)
     *   (ref1, 8, Forward)
     */
    {
        const std::string pattern = "ACG";

        const SA_Range range =
            search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        assert_expected_hits(
            hits,
            {
                {0, "ref1", 0, Strand::Forward},
                {0, "ref1", 1, Strand::Reverse},
                {0, "ref1", 8, Strand::Forward}
            });
    }

    /*
     * ------------------------------------------------------------
     * 2. TGC
     * ------------------------------------------------------------
     *
     * ref1 = A C G T T G C A A C G
     *              ^^^
     *              4
     *
     * TGC occurs forward at position 4.
     *
     * RC(TGC) = GCA, which occurs at position 5:
     *
     * ref1 = A C G T T G C A A C G
     *                ^^^
     *                5
     *
     * Therefore:
     *
     *   (ref1, 4, Forward)
     *   (ref1, 5, Reverse)
     */
    {
        const std::string pattern = "TGC";

        const SA_Range range =
            search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        assert_expected_hits(
            hits,
            {
                {0, "ref1", 4, Strand::Forward},
                {0, "ref1", 5, Strand::Reverse}
            });
    }

    /*
     * ------------------------------------------------------------
     * 3. TTACC on ref2
     * ------------------------------------------------------------
     *
     * ref2 = T T A C C G G T A A
     *        ^^^^^
     *        0
     *
     * TTACC occurs only in the forward orientation.
     */
    {
        const std::string pattern = "TTACC";

        const SA_Range range =
            search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

            assert_expected_hits( hits, { {1, "ref2", 0, Strand::Forward}, {1, "ref2", 5, Strand::Reverse} });
    }

    /*
     * ------------------------------------------------------------
     * 4. CGG on ref2
     * ------------------------------------------------------------
     *
     * ref2 = T T A C C G G T A A
     *              ^^^
     *              4
     *
     * CGG occurs forward at position 4.
     *
     * RC(CGG) = CCG, which occurs at position 3.
     *
     * Therefore:
     *
     *   (ref2, 4, Forward)
     *   (ref2, 3, Reverse)
     */
    {
        const std::string pattern = "CGG";

        const SA_Range range =
            search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        assert_expected_hits(
            hits,
            {
                {1, "ref2", 4, Strand::Forward},
                {1, "ref2", 3, Strand::Reverse}
            });
    }

    /*
     * ------------------------------------------------------------
     * 5. Empty pattern must fail
     * ------------------------------------------------------------
     */
    {
        const SA_Range range =
            index.initial_range(A);

        bool threw = false;

        try {
            (void)index.locate(range, 0);
        }
        catch (const std::invalid_argument&) {
            threw = true;
        }

        assert(threw);
    }

    /*
     * ------------------------------------------------------------
     * 6. Metadata-less index must fail
     * ------------------------------------------------------------
     */
    {
        BwaFMDIndex no_metadata(bwt);

        const SA_Range range =
            search(no_metadata, "ACG");

        bool threw = false;

        try {
            (void)no_metadata.locate(
                range,
                3);
        }
        catch (const std::logic_error&) {
            threw = true;
        }

        assert(threw);
    }

    /*
     * ------------------------------------------------------------
     * 7. Basic raw locate remains unchanged
     * ------------------------------------------------------------
     */
    {
        BwaFMDIndex no_metadata(bwt);

        const SA_Range range =
            search(no_metadata, "ACG");

        const auto positions =
            no_metadata.locate(range);

        assert(!positions.empty());
    }

    bns_destroy(bns);
    bwt_destroy(bwt);

    std::cout
        << "All metadata-aware locate tests passed.\n";

    return 0;
}
