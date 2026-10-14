#include "exact_search.hpp"
#include "bwt.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

char complement(char c)
{
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    default: std::abort();
    }
}

std::string reverse_complement(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (auto it = s.rbegin(); it != s.rend(); ++it)
        out.push_back(complement(*it));
    return out;
}

std::vector<uint64_t> occurrences(const std::string& text,
                                  const std::string& pattern)
{
    std::vector<uint64_t> out;
    if (pattern.size() > text.size()) return out;
    for (size_t i = 0; i + pattern.size() <= text.size(); ++i)
        if (text.compare(i, pattern.size(), pattern) == 0)
            out.push_back(static_cast<uint64_t>(i));
    return out;
}

void enumerate(const BwaFMDIndex& index, const std::string& text,
               size_t target_length, std::string& pattern,
               uint64_t& tested)
{
    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    if (pattern.size() == target_length) {
        const SA_Range right = search_exact_unidirectional(
            index, pattern, Direction::Right);
        const SA_Range left = search_exact_unidirectional(
            index, pattern, Direction::Left);

        const SA_Range bidirectional_right =
            search_exact_bidirectional_seed(index, pattern, Direction::Right);
        const SA_Range bidirectional_left =
            search_exact_bidirectional_seed(index, pattern, Direction::Left);
        const Interval right_expected = bidirectional_right.companion_interval();
        const Interval left_expected = bidirectional_left.primary_interval();

        const bool interval_mismatch =
            !right.is_single() ||
            right.orientation() != SA_Range::Orientation::Right ||
            !left.is_single() ||
            left.orientation() != SA_Range::Orientation::Left ||
            right.size() != bidirectional_right.size() ||
            left.size() != bidirectional_left.size() ||
            (!right.empty() &&
                (right.interval().l != right_expected.l ||
                 right.interval().r != right_expected.r)) ||
            (!left.empty() &&
                (left.interval().l != left_expected.l ||
                 left.interval().r != left_expected.r));

        if (interval_mismatch) {
            std::cerr << "Exact-search mismatch for pattern: " << pattern << '\n'
                      << "  right size: " << right.size()
                      << " expected: " << bidirectional_right.size() << '\n'
                      << "  left size:  " << left.size()
                      << " expected: " << bidirectional_left.size() << '\n';
            std::abort();
        }

        ++tested;
        return;
    }

    for (char c : bases) {
        pattern.push_back(c);
        enumerate(index, text, target_length, pattern, tested);
        pattern.pop_back();
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <index.bwt> <index.sa>\n";
        return 2;
    }

    const std::string reference = "ACGTTGCAACG";
    const std::string indexed_text = reference + reverse_complement(reference);

    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) {
        std::cerr << "Could not restore BWT: " << argv[1] << '\n';
        return 1;
    }
    bwt_restore_sa(argv[2], bwt);
    if (!bwt->sa) {
        std::cerr << "Could not restore suffix array: " << argv[2] << '\n';
        bwt_destroy(bwt);
        return 1;
    }

    if (bwt->seq_len != indexed_text.size()) {
        std::cerr << "Test fixture does not match index length\n";
        bwt_destroy(bwt);
        return 1;
    }

    BwaFMDIndex index(bwt);
    uint64_t tested = 0;
    for (size_t length = 1; length <= 6; ++length) {
        std::string pattern;
        pattern.reserve(length);
        enumerate(index, indexed_text, length, pattern, tested);
        std::cout << "Length " << length << ": passed\n";
    }

    bwt_destroy(bwt);
    std::cout << "Fast exact search passed for " << tested
              << " patterns in both directions.\n";
    return 0;
}
