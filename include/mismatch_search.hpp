#pragma once

#include "bwa_fmd_index.hpp"
#include <string>
#include <vector>

// Search from a seed. seed.position() is the zero-based read position
// already consumed; Direction::Right scans left-to-right, Direction::Left
// scans right-to-left. Report hits with exactly Count mismatches.
//
// Rejects is an ordered collection of conservative deferred/rejected search
// states, not a promise that each state will be explored or produce a hit.
// The search may stop with states unused (for example, at the end of the
// scanned read or when MAXHITS is reached); callers may use them in a later
// pass, but this function does not automatically resume them. The collection
// can mix distinct state types: at the mismatch threshold it may contain a
// copy of the parent range with the expected base forbidden, preserving
// possible alternative branches; it may also contain an actual child range
// whose mismatch count exceeds Count. These entries retain the interval,
// paired-range, position, validity and mismatch metadata appropriate to their
// construction, and forbidden_base is meaningful for the threshold parent
// state rather than a general marker on every reject. Conservative entries
// may be redundant (including around the sentinel). Hit_List and Rejects are
// appended in traversal order; their ordering is part of the observable
// result and callers should not assume Rejects is a deduplicated list.
void search_unidirectional(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases = false);

void extend(
    const std::string& read,
    const SA_Range& seed,
    int Count,
    int MAXHITS,
    std::vector<SA_Range>& Hit_List,
    std::vector<SA_Range>& Rejects,
    Direction direction,
    bool capture_reference_bases = false);