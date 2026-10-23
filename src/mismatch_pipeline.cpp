#include "mismatch_pipeline.hpp"
#include "exact_search.hpp"
#include "mismatch_search.hpp"

#include <cctype>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace {

std::string normalize_read(const std::string& read) {
    std::string result = read;
    for (char& c : result)
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return result;
}

void bisect_read(const std::string& read, std::string& left, std::string& right) {
    const size_t mid = read.size() / 2;
    left = read.substr(0, mid);
    right = read.substr(mid);
}

SA_Range exact_seed(const BwaFMDIndex& index, const std::string& seed,
                    Direction direction, ZeroMismatchMode mode) {
    return mode == ZeroMismatchMode::Bidirectional
        ? search_exact_bidirectional(index, seed, direction)
        : search_exact_unidirectional(index, seed, direction);
}

void append(std::vector<SA_Range>& destination, std::vector<SA_Range>&& source) {
    destination.insert(destination.end(),
                       std::make_move_iterator(source.begin()),
                       std::make_move_iterator(source.end()));
}

// Stage 1: exact LH seed, then extend through RH with zero mismatches.
// The mismatch search routine records the frontier states needed by stage 2.
void run_zero_mismatch(const std::string& read, const BwaFMDIndex& index,
                       MismatchPipelineState& state, ZeroMismatchMode mode) {
    if (read.empty()) return;

    std::string left, right;
    bisect_read(read, left, right);
    if (left.empty()) {
        // A one-base read has no non-empty LH seed. Search it directly.
        SA_Range hit = exact_seed(index, read, Direction::Right, mode);
        if (!hit.empty()) state.hits.push_back(std::move(hit));
        return;
    }

    SA_Range lh_seed = exact_seed(index, left, Direction::Right, mode);
    if (lh_seed.empty()) return;

    std::vector<SA_Range> hits, rejects;
    extend(read, lh_seed, 0, 0, hits, rejects, Direction::Right);
    append(state.hits, std::move(hits));
    append(state.rejects_after_0mis, std::move(rejects));
}

// Stage 2a: continue the zero-mismatch reject frontier through the right
// half with one mismatch allowed. The rejects emitted here are the deferred
// frontier for a future two-mismatch pass.
void run_one_mismatch_right(const std::string& read,
                            MismatchPipelineState& state, bool all_hits) {
    const int max_hits = all_hits ? 0 : 1;
    for (const SA_Range& seed : state.rejects_after_0mis) {
        std::vector<SA_Range> hits, rejects;
        extend(read, seed, 1, max_hits, hits, rejects, Direction::Right, true);
        append(state.hits, std::move(hits));
        append(state.rejects_after_1mis, std::move(rejects));
        if (!all_hits && !state.hits.empty()) return;
    }
}

// Stage 2b: independently seed from an exact RH match, then extend left
// through LH with one mismatch allowed. This covers mismatches in LH without
// depending on the rightward zero-mismatch reject frontier.
void run_one_mismatch_left(const std::string& read, const BwaFMDIndex& index,
                           MismatchPipelineState& state, bool all_hits,
                           ZeroMismatchMode mode) {
    std::string left, right;
    bisect_read(read, left, right);
    if (left.empty() || right.empty()) return;

    SA_Range rh_seed = exact_seed(index, right, Direction::Left, mode);
    if (rh_seed.empty()) return;

    std::vector<SA_Range> hits, rejects;
    const int max_hits = all_hits ? 0 : 1;
    extend(read, rh_seed, 1, max_hits, hits, rejects, Direction::Left, true);
    append(state.hits, std::move(hits));
    append(state.rejects_after_1mis, std::move(rejects));
}

} // namespace

void Run_Zero_Mismatch(const std::string& read, const BwaFMDIndex* index,
                       MismatchPipelineState& state, ZeroMismatchMode mode) {
    if (!index) throw std::invalid_argument("Run_Zero_Mismatch: null index");
    const std::string normalized = normalize_read(read);
    if (normalized.empty()) return;
    run_zero_mismatch(normalized, *index, state, mode);
}

void Run_Mismatch_Pipeline(const std::string& read, const BwaFMDIndex* index,
                           MismatchPipelineState& state, int max_mismatches,
                           bool do_all_hits, ZeroMismatchMode zero_mode) {
    state.clear();
    if (!index) throw std::invalid_argument("Run_Mismatch_Pipeline: null index");
    if (max_mismatches < 0 || max_mismatches > 1)
        throw std::invalid_argument(
            "Run_Mismatch_Pipeline currently supports max_mismatches 0 or 1");

    const std::string normalized = normalize_read(read);
    if (normalized.empty()) return;

    run_zero_mismatch(normalized, *index, state, zero_mode);
    if (max_mismatches == 0 || (!do_all_hits && !state.hits.empty())) return;

    run_one_mismatch_right(normalized, state, do_all_hits);
    if (!do_all_hits && !state.hits.empty()) return;
    run_one_mismatch_left(normalized, *index, state, do_all_hits, zero_mode);
}
