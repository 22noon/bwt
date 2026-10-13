#pragma once

#include "bwa_fmd_index.hpp"
#include <string>

/*
 * Exact unidirectional search of a read or seed substring.
 *
 * Direction::Right scans read[0] .. read[n-1].
 * Direction::Left scans read[n-1] .. read[0].
 *
 * The returned SA_Range stores position() == read.size() - 1 for a
 * completed non-empty search, so it can be passed to the mismatch
 * extension routine together with the full read.  position() is the
 * number of consumed characters minus one in scan coordinates, not a
 * natural left-to-right read coordinate.
 *
 * Once an interval becomes a singleton, subsequent extensions use the
 * dedicated singleton fast path.  No branching or reject ranges are
 * generated here.
 */
SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);
