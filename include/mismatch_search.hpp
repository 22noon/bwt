#pragma once

#include "bwa_fmd_index.hpp"
#include <string>
#include <vector>

// Search from a seed. seed.position() is the zero-based read position
// already consumed; Direction::Right scans left-to-right, Direction::Left
// scans right-to-left. Report hits with exactly Count mismatches.
void search_unidirectional(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases = false);
