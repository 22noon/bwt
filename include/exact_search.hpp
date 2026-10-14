#pragma once

#include "bwa_fmd_index.hpp"
#include <string>

/*
 * Exact unidirectional search using a single SA interval.
 *
 * Direction::Right scans read[0] .. read[n-1] and stores the interval
 * on the reverse-complement side (Orientation::Right).
 * Direction::Left scans read[n-1] .. read[0] and stores the primary
 * interval (Orientation::Left).
 *
 * position() is the number of consumed characters minus one in scan
 * coordinates. No companion interval or branches are constructed.
 */
SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);

/*
 * Exact search that returns a paired FMD range suitable for the current
 * mismatch-search/branching pipeline. Use this when the result will be
 * passed to search_unidirectional().
 */
SA_Range search_exact_bidirectional_seed(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);
