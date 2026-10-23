#include "batmis_rank_backend.hpp"

#include <algorithm>
#include <limits>
#if defined(__SSE2__)
#include <emmintrin.h>
#endif
#include <stdexcept>

namespace {
constexpr uint64_t kCheckpoint = 256;
} // namespace

BatMisRankBackend::BatMisRankBackend(const bwt_t* bwt)
{
    if (!bwt || bwt->seq_len == 0 || bwt->primary > bwt->seq_len ||
        bwt->seq_len > std::numeric_limits<uint32_t>::max())
        throw std::invalid_argument("BatMIS rank backend: invalid BWA BWT");

    primary_ = bwt->primary;
    packed_length_ = bwt->seq_len; // BWA omits the primary row but stores N symbols.
    // BatMIS's SIMD decoder reads 64-symbol blocks around an endpoint.
    // Keep four zero-padded words beyond the packed payload for terminal reads.
    packed_.assign(static_cast<size_t>((packed_length_ + 63) / 16 + 4), 0);

    const uint64_t checkpoint_count = packed_length_ / kCheckpoint +
        (packed_length_ % kCheckpoint != 0 ? 1 : 0) + 1;
    checkpoints_.resize(static_cast<size_t>(checkpoint_count));
    std::array<uint32_t, 4> counts{};
    checkpoints_[0] = counts;

    for (uint64_t i = 0; i < packed_length_; ++i) {
        const uint8_t c = static_cast<uint8_t>(bwt_B0(bwt, i));
        packed_[static_cast<size_t>(i / 16)] |=
            static_cast<uint32_t>(c) << (30 - 2 * (i % 16));
        ++counts[c];
        if ((i + 1) % kCheckpoint == 0)
            checkpoints_[static_cast<size_t>((i + 1) / kCheckpoint)] = counts;
    }
    if (packed_length_ % kCheckpoint != 0)
        checkpoints_.back() = counts;
}

uint32_t BatMisRankBackend::decode(
    uint32_t index1,
    uint32_t index2,
    uint8_t character) const noexcept
{
#if defined(__SSE2__)
    uint32_t numChar1, numChar2, minIndex, maxIndex, minIndex128, maxIndex128;
    uint32_t r;
    alignas(16) static const uint32_t partition_one_1[4] = {47, 31, 15, 0};
    alignas(16) static const uint32_t partition_one_2[4] = {0, 15, 31, 47};
    alignas(16) static const uint32_t partition_zero_1[4] = {63, 47, 31, 15};
    alignas(16) static const uint32_t partition_zero_2[4] = {15, 31, 47, 63};

    r = (index1 - index2) & -(index1 < index2);
    minIndex = index2 + r;
    maxIndex = index1 - r;
    minIndex128 = minIndex & ~uint32_t{63};
    maxIndex128 = (maxIndex - (maxIndex - minIndex > 64)) & ~uint32_t{63};
    numChar1 = maxIndex128 - minIndex;
    numChar2 = maxIndex - maxIndex128;

    const __m128i r1e = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        packed_.data() + minIndex128 / 16));
    const __m128i r2e = _mm_loadu_si128(reinterpret_cast<const __m128i*>(
        packed_.data() + maxIndex128 / 16));
    __m128i m0 = _mm_set1_epi32(0xffffffffu + (character & 1u));
    __m128i m1 = _mm_set1_epi32(0xffffffffu + (character >> 1));
    __m128i mcl = _mm_set1_epi32(0x55555555);
    __m128i r1a = _mm_set1_epi32(static_cast<int>(numChar1));
    __m128i r2a = _mm_set1_epi32(static_cast<int>(numChar2));
    __m128i r1b = _mm_load_si128(reinterpret_cast<const __m128i*>(partition_one_1));
    __m128i r2b = _mm_load_si128(reinterpret_cast<const __m128i*>(partition_one_2));
    __m128i r1c = _mm_load_si128(reinterpret_cast<const __m128i*>(partition_zero_1));
    __m128i r2c = _mm_load_si128(reinterpret_cast<const __m128i*>(partition_zero_2));
    r1b = _mm_cmpgt_epi32(r1a, r1b);
    r2b = _mm_cmpgt_epi32(r2a, r2b);
    r1c = _mm_cmpgt_epi32(r1a, r1c);
    r2c = _mm_cmpgt_epi32(r2a, r2c);
    r1b = _mm_srl_epi32(r1b, _mm_cvtsi32_si128(static_cast<int>((16 - numChar1 % 16) * 2)));
    r2b = _mm_sll_epi32(r2b, _mm_cvtsi32_si128(static_cast<int>((16 - numChar2 % 16) * 2)));
    r1c = _mm_and_si128(_mm_or_si128(r1b, r1c), mcl);
    r2c = _mm_and_si128(_mm_or_si128(r2b, r2c), mcl);
    r1b = _mm_xor_si128(_mm_srli_epi32(r1e, 1), m1);
    r2b = _mm_xor_si128(_mm_srli_epi32(r2e, 1), m1);
    r1a = _mm_and_si128(_mm_xor_si128(r1e, m0), r1b);
    r2a = _mm_and_si128(_mm_xor_si128(r2e, m0), r2b);
    r1a = _mm_add_epi32(_mm_and_si128(r1a, r1c),
                         _mm_and_si128(r2a, r2c));
    mcl = _mm_set1_epi32(0x33333333);
    r1b = _mm_and_si128(_mm_srli_epi32(r1a, 2), mcl);
    r1a = _mm_add_epi32(_mm_and_si128(r1a, mcl), r1b);
    r1a = _mm_and_si128(_mm_add_epi32(r1a, _mm_srli_epi32(r1a, 4)),
                         _mm_set1_epi32(0x0f0f0f0f));
    r1a = _mm_sad_epu8(r1a, _mm_setzero_si128());
    return static_cast<uint32_t>(_mm_extract_epi16(r1a, 0) +
                                 _mm_extract_epi16(r1a, 4));
#else
    const uint32_t begin = std::min(index1, index2);
    const uint32_t end = std::max(index1, index2);
    uint32_t count = 0;
    for (uint32_t i = begin; i < end; ++i) {
        const uint32_t word = packed_[i / 16];
        count += ((word >> (30 - 2 * (i % 16))) & 3u) == c;
    }
    return count;
#endif
}

uint64_t BatMisRankBackend::occ(uint64_t index, uint8_t c) const noexcept
{
    if (index > primary_)
        --index;
    if (index > packed_length_)
        index = packed_length_;

    // BatMIS rounds to the nearest 256-symbol checkpoint, choosing the lower
    // checkpoint at an exact midpoint, then counts forward or backward.
    const uint64_t checkpoint = (index + 127) / kCheckpoint;
    const uint64_t checkpoint_pos = std::min(
        checkpoint * kCheckpoint, packed_length_);
    const uint64_t base = checkpoints_[static_cast<size_t>(checkpoint)][c];
    if (checkpoint_pos == index)
        return base;
    const uint32_t decoded = decode(
        static_cast<uint32_t>(checkpoint_pos),
        static_cast<uint32_t>(index), c);
    if (checkpoint_pos < index)
        return base + decoded;
    return base - decoded;
}

bool BatMisRankBackend::symbol_at(uint64_t row, uint8_t& c) const noexcept
{
    if (row > packed_length_ || row == primary_) return false;
    const uint64_t packed_row = row - (row > primary_);
    if (packed_row >= packed_length_) return false;
    const uint32_t word = packed_[static_cast<size_t>(packed_row / 16)];
    c = static_cast<uint8_t>((word >> (30 - 2 * (packed_row % 16))) & 3u);
    return true;
}
