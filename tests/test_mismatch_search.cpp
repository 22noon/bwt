#include "mismatch_search.hpp"

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

uint8_t code(char c)
{
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:
        std::cerr << "Invalid DNA base: " << c << '\n';
        std::abort();
    }
}

SA_Range exact_range(
    const BwaFMDIndex& index,
    const std::string& pattern,
    Direction direction)
{
    assert(!pattern.empty());

    SA_Range range;
    if (direction == Direction::Right) {
        range = index.initial_range(code(pattern.front()));
        for (size_t i = 1; i < pattern.size(); ++i)
            range = index.extend_right(range, code(pattern[i]));
    } else {
        range = index.initial_range(code(pattern.back()));
        for (size_t i = pattern.size() - 1; i-- > 0;)
            range = index.extend_left(range, code(pattern[i]));
    }
    return range;
}

void generate_patterns(
    size_t length,
    std::string& current,
    std::vector<std::string>& output)
{
    if (current.size() == length) {
        output.push_back(current);
        return;
    }

    static constexpr char bases[] = {'A', 'C', 'G', 'T'};
    for (char c : bases) {
        current.push_back(c);
        generate_patterns(length, current, output);
        current.pop_back();
    }
}

unsigned hamming_distance(
    const std::string& a,
    const std::string& b)
{
    assert(a.size() == b.size());
    unsigned distance = 0;
    for (size_t i = 0; i < a.size(); ++i)
        distance += (a[i] != b[i]);
    return distance;
}


uint64_t expected_count(
    const BwaFMDIndex& index,
    const std::string& query,
    unsigned mismatches,
    Direction direction,
    const std::vector<std::string>& patterns)
{
    uint64_t count = 0;

    // The first consumed base is already fixed by the seed.
    const size_t seed_pos =
        direction == Direction::Right ? 0 : query.size() - 1;

    for (const auto& candidate : patterns) {
        if (candidate[seed_pos] != query[seed_pos])
            continue;

        if (hamming_distance(query, candidate) == mismatches)
            count += exact_range(index, candidate, direction).size();
    }

    return count;
}


uint64_t observed_count(
    const BwaFMDIndex& index,
    const std::string& query,
    unsigned mismatches,
    Direction direction)
{
    // The seed represents the first character consumed in the selected
    // search direction. position() is the zero-based consumed position.
    const char first = direction == Direction::Right
        ? query.front()
        : query.back();

    SA_Range seed = index.initial_range(code(first));
    seed.set_position(0);
    seed.set_valid(true);

    std::vector<SA_Range> hits;
    std::vector<SA_Range> rejects;
    search_unidirectional(
        query, seed, static_cast<int>(mismatches), 0,
        hits, rejects, direction);

    uint64_t count = 0;
    for (const auto& hit : hits) {
        assert(!hit.empty());
        assert(hit.mismatch_count() == mismatches);
        count += hit.size();
    }
    return count;
}

void test_query(
    const BwaFMDIndex& index,
    const std::string& query,
    const std::vector<std::string>& patterns)
{
    for (Direction direction : {Direction::Right, Direction::Left}) {
        for (unsigned mismatches = 0; mismatches <= 1; ++mismatches) {
            const uint64_t expected = expected_count(
                index, query, mismatches, direction, patterns);
            const uint64_t observed = observed_count(
                index, query, mismatches, direction);

            if (expected != observed) {
                std::cerr
                    << "Mismatch-search count failure\n"
                    << "  query: " << query << '\n'
                    << "  direction: "
                    << (direction == Direction::Right ? "right" : "left")
                    << '\n'
                    << "  mismatches: " << mismatches << '\n'
                    << "  expected total interval size: " << expected << '\n'
                    << "  observed total interval size: " << observed << '\n';
                std::abort();
            }
        }
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <index.bwt> <index.sa>\n";
        return 2;
    }

    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) {
        std::cerr << "Failed to restore BWT: " << argv[1] << '\n';
        return 1;
    }

    bwt_restore_sa(argv[2], bwt);
    if (!bwt->sa) {
        std::cerr << "Failed to restore SA: " << argv[2] << '\n';
        bwt_destroy(bwt);
        return 1;
    }

    BwaFMDIndex index(bwt);
    size_t queries_tested = 0;

    // Length 1 is deliberately omitted: the API seed is already one
    // consumed character, so there is no position at which to introduce
    // a mismatch for a one-character query.
    for (size_t length = 2; length <= 4; ++length) {
        std::vector<std::string> patterns;
        std::string current;
        generate_patterns(length, current, patterns);

        for (const auto& query : patterns) {
            test_query(index, query, patterns);
            ++queries_tested;
        }

        std::cout << "Mismatch-search length " << length << ": passed\n";
    }

    bwt_destroy(bwt);
    std::cout << "Mismatch-search tests passed; queries tested: "
              << queries_tested << '\n';
    return 0;
}
