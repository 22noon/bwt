#pragma once
#include "bwa_fmd_index.hpp"
#include <vector>

enum class ZeroMismatchMode { Fast, Bidirectional };

struct MismatchPipelineState {
    std::vector<SA_Range> hits;
    std::vector<SA_Range> rejects_after_0mis;
    std::vector<SA_Range> rejects_after_1mis;
    void clear() { hits.clear(); rejects_after_0mis.clear(); rejects_after_1mis.clear(); }
};

void Run_Zero_Mismatch(const std::string& read, const BwaFMDIndex* index,
                       MismatchPipelineState& state, ZeroMismatchMode mode);
void Run_Mismatch_Pipeline(const std::string& read, const BwaFMDIndex* index,
                           MismatchPipelineState& state, int max_mismatches,
                           bool do_all_hits, ZeroMismatchMode zero_mode);
