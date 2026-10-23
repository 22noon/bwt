#pragma once

#include <array>
#include <cstdint>
#include <vector>

extern "C" {
#include "bwt.h"
}

/*
 * BatMIS-style rank support built in memory from a decoded BWA bwt_t.
 * The BWA serialized BWT is never reinterpreted as a BatMIS file or struct.
 * Packed symbols use BatMIS's high-to-low 2-bit ordering and occurrence
 * checkpoints are placed every 256 symbols, matching its rank convention.
 */
class BatMisRankBackend {
public:
    explicit BatMisRankBackend(const bwt_t* bwt);

    // Equivalent to BatMIS BWTOccValue(bwt, index, c): count symbols before
    // the logical BWT row `index`, with the omitted primary row accounted for.
    uint64_t occ(uint64_t index, uint8_t c) const noexcept;

    // Return the packed BWT symbol at a logical row. The primary row is
    // omitted from the packed sequence and therefore returns false.
    bool symbol_at(uint64_t row, uint8_t& c) const noexcept;

    uint64_t primary() const noexcept { return primary_; }
    uint64_t packed_length() const noexcept { return packed_length_; }

private:
    uint64_t primary_ = 0;
    uint64_t packed_length_ = 0;
    std::vector<uint32_t> packed_;
    std::vector<std::array<uint32_t, 4>> checkpoints_;

    uint32_t decode(uint32_t index1, uint32_t index2, uint8_t c) const noexcept;
};
