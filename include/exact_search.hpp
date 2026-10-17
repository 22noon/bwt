#pragma once

#include "bwa_fmd_index.hpp"
#include <string>

// Pure exact search returning one SA interval.
SA_Range search_exact_unidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);

// Pure exact search returning a paired FMD interval.
SA_Range search_exact_bidirectional(
    const BwaFMDIndex& index,
    const std::string& read,
    Direction direction);

// Extend an existing single SA interval by the supplied sequence.
// Direction::Right consumes extension left-to-right.
// Direction::Left consumes extension right-to-left.
SA_Range extend_exact_unidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction);

// Extend an existing paired FMD interval by the supplied sequence.
// Uses singleton-cardinality and sentinel-aware single-path
// optimizations before falling back to general FMD extension.
SA_Range extend_exact_bidirectional(
    const BwaFMDIndex& index,
    const SA_Range& seed,
    const std::string& extension,
    Direction direction);

