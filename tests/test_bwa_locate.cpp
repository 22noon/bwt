#include "bwa_fmd_index.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
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

const LocatedHit& only_hit(
    const std::vector<LocatedHit>& hits)
{
    assert(hits.size() == 1);
    return hits.front();
}

void check_hit(
    const LocatedHit& hit,
    uint32_t ref_id,
    const std::string& name,
    uint64_t position,
    Strand strand)
{
    assert(hit.ref_id == ref_id);
    assert(hit.ref_name == name);
    assert(hit.position == position);
    assert(hit.strand == strand);
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
    const std::string ann_file = prefix + ".ann";
    const std::string amb_file = prefix + ".amb";
    const std::string pac_file = prefix + ".pac";

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
     *
     * and the natural forward reference coordinate space is:
     *
     *   [0, 21)
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
     * Load the BWA reference metadata.
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
     * 1. Forward-strand hit
     * ------------------------------------------------------------
     *
     * ref1 = ACGTTGCAACG
     *          ^^^
     * ACG occurs at position 0.
     */
    {
        const std::string pattern = "ACG";
        const SA_Range range = search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        assert(hits.size() == 3);

        bool found_forward_0 = false;
        bool found_reverse_1 = false;
        bool found_forward_8 = false;

        for (const auto& hit : hits) {
            assert(hit.ref_id == 0);
            assert(hit.ref_name == "ref1");

            if (hit.position == 0 &&
                hit.strand == Strand::Forward)
            {
                found_forward_0 = true;
            }
            else if (hit.position == 1 &&
                    hit.strand == Strand::Reverse)
            {
                found_reverse_1 = true;
            }
            else if (hit.position == 8 &&
                    hit.strand == Strand::Forward)
            {
                found_forward_8 = true;
            }
            else {
                std::cerr
                    << "Unexpected hit: position="
                    << hit.position
                    << " strand="
                    << (hit.strand == Strand::Forward
                            ? "Forward"
                            : "Reverse")
                    << '\n';

                assert(false);
            }
        }

        assert(found_forward_0);
        assert(found_reverse_1);
        assert(found_forward_8);
    }
    /*
     * ------------------------------------------------------------
     * 2. Reverse-strand hit
     * ------------------------------------------------------------
     *
     * ref1 contains:
     *
     *   ...TGCA...
     *
     * The reverse-complement of "GCA" is "TGC".
     *
     * Searching TGC therefore gives a reverse-strand occurrence
     * corresponding to ref1 position 4.
     */
    {
        const std::string pattern = "TGC";
        const SA_Range range = search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        bool found_reverse = false;

        for (const auto& hit : hits) {
            if (hit.ref_id == 0 &&
                hit.position == 4 &&
                hit.strand == Strand::Reverse)
            {
                found_reverse = true;
                assert(hit.ref_name == "ref1");
            }
        }

        assert(found_reverse);
    }

    /*
     * ------------------------------------------------------------
     * 3. Second reference
     * ------------------------------------------------------------
     */
    {
        const std::string pattern = "TTACC";
        const SA_Range range = search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        assert(hits.size() == 1);

        check_hit(
            only_hit(hits),
            1,
            "ref2",
            0,
            Strand::Forward);
    }

    /*
     * ------------------------------------------------------------
     * 4. Reverse-strand hit on ref2
     * ------------------------------------------------------------
     *
     * ref2 = TTACCGGTAA
     *
     * Reverse complement of "CCG" is "CGG".
     *
     * CCG occurs at ref2 position 3, so CGG should locate to:
     *
     *   ref2 position 3, Reverse
     */
    {
        const std::string pattern = "CGG";
        const SA_Range range = search(index, pattern);

        const auto hits =
            index.locate(range, pattern.size());

        bool found_reverse = false;

        for (const auto& hit : hits) {
            if (hit.ref_id == 1 &&
                hit.position == 3 &&
                hit.strand == Strand::Reverse)
            {
                found_reverse = true;
                assert(hit.ref_name == "ref2");
            }
        }

        assert(found_reverse);
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
