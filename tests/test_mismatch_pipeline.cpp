#include "exact_search.hpp"
#include "mismatch_pipeline.hpp"
#include "bwt.h"

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool same_interval(const Interval& a, const Interval& b)
{
    return a.l == b.l && a.r == b.r;
}

void check_read(
    const BwaFMDIndex& index,
    const std::string& read,
    ZeroMismatchMode mode)
{
    MismatchPipelineState state;

    Run_Mismatch_Pipeline(
        read, &index, state, 0, false, mode);

    if (mode == ZeroMismatchMode::Fast) {
        const SA_Range expected = search_exact_unidirectional(
            index, read, Direction::Right);

        assert(state.hits.size() == (expected.empty() ? 0u : 1u));

        if (!expected.empty()) {
            assert(state.hits.front().is_single());
            assert(state.hits.front().orientation() ==
                   SA_Range::Orientation::Right);
            assert(same_interval(
                state.hits.front().interval(),
                expected.interval()));
        }

        for (const SA_Range& reject : state.rejects_after_0mis) {
            assert(reject.is_single());
            assert(reject.orientation() ==
                   SA_Range::Orientation::Right);
            assert(reject.has_forbidden_base());
        }
    } else {
        const SA_Range expected = search_exact_bidirectional(
            index, read, Direction::Right);

        assert(state.hits.size() == (expected.empty() ? 0u : 1u));

        if (!expected.empty()) {
            assert(state.hits.front().is_bidirectional());
            assert(same_interval(
                state.hits.front().primary_interval(),
                expected.primary_interval()));
            assert(same_interval(
                state.hits.front().companion_interval(),
                expected.companion_interval()));
        }

        for (const SA_Range& reject : state.rejects_after_0mis)
            assert(reject.is_bidirectional());
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <index.bwt> <index.sa>\n";
        return 2;
    }

    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) {
        std::cerr << "Could not restore BWT\n";
        return 1;
    }

    bwt_restore_sa(argv[2], bwt);
    if (!bwt->sa) {
        std::cerr << "Could not restore suffix array\n";
        bwt_destroy(bwt);
        return 1;
    }

    BwaFMDIndex index(bwt);

    const std::vector<std::string> reads = {
        "A",
        "AC",
        "ACG",
        "TTG",
        "AAAA",
        "TTTTTT",
        "ACGTTGCAACG"
    };

    for (const std::string& read : reads) {
        check_read(index, read, ZeroMismatchMode::Fast);
        check_read(index, read, ZeroMismatchMode::Bidirectional);
    }

    bwt_destroy(bwt);

    std::cout << "Zero-mismatch pipeline tests passed\n";
    return 0;
}
