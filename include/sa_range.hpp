#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>
#include <string>

class BwaFMDIndex;
class BranchSet;

enum class Direction {
    Left,
    Right
};

enum class Strand {
    Forward,
    Reverse
};

struct LocatedHit {
    uint32_t ref_id;       // Zero-based index into bntseq_t::anns
    std::string ref_name;  // Reference/contig name
    uint64_t position;     // Zero-based leftmost position within reference
    Strand strand;
};

struct Mismatch {
    uint32_t Pos = 0;
    char Char = 0;
};

struct Interval {
    uint64_t l = 0;
    uint64_t r = 0;  // half-open [l,r)

    uint64_t size() const noexcept { return r - l; }
    bool empty() const noexcept { return l == r; }
};

class SA_Range {
public:
    enum class Mode { Single, Bidirectional };
    enum class Orientation { Left, Right };

private:
    const BwaFMDIndex* index_ = nullptr;

    Mode mode_ = Mode::Single;
    Orientation orientation_ = Orientation::Left;

    Interval interval_;
    Interval primary_;
    Interval companion_;

    uint32_t Pos = 0;
    bool Valid = false;

    uint32_t Mismatches = 0;
    std::vector<Mismatch> MM;

    char forbidden_char = '\0';
    bool unidirectional_search = false;

#ifdef DEBUG
    std::string Scanned_Sequence;
#endif

    void copy_search_state_to(SA_Range& destination) const;

public:
    SA_Range() = default;

    static SA_Range single(
        const BwaFMDIndex* index,
        Interval interval,
        Orientation orientation = Orientation::Left);

    static SA_Range bidirectional(
        const BwaFMDIndex* index,
        Interval primary,
        Interval companion);

    // Accessors
    const BwaFMDIndex* index() const noexcept { return index_; }

    bool valid() const noexcept { return Valid; }
    void set_valid(bool value) noexcept { Valid = value; }

    uint32_t position() const noexcept { return Pos; }
    void set_position(uint32_t pos) noexcept { Pos = pos; }
    void advance_position() noexcept { ++Pos; }

    uint32_t mismatch_count() const noexcept { return Mismatches; }

    const std::vector<Mismatch>& mismatches() const noexcept { return MM; }
    std::vector<Mismatch>& mismatches() noexcept { return MM; }

    void add_mismatch(uint32_t pos, char query_base)
    {
        MM.push_back({pos, query_base});
        Mismatches = static_cast<uint32_t>(MM.size());
    }

    char forbidden_base() const noexcept { return forbidden_char; }
    void set_forbidden_base(char base) noexcept { forbidden_char = base; }
    bool has_forbidden_base() const noexcept { return forbidden_char != '\0'; }

    bool is_unidirectional() const noexcept { return unidirectional_search; }
    void set_unidirectional(bool value) noexcept { unidirectional_search = value; }

#ifdef DEBUG
    const std::string& scanned_sequence() const noexcept { return Scanned_Sequence; }
    std::string& scanned_sequence() noexcept { return Scanned_Sequence; }
#endif

    bool is_single() const noexcept;
    bool is_bidirectional() const noexcept;
    Orientation orientation() const noexcept;

    Interval interval() const;
    Interval primary_interval() const;
    Interval companion_interval() const;

    uint64_t size() const;
    bool empty() const;

    SA_Range to_single(Orientation orientation) const;

    // Search methods
    SA_Range extend_left(uint8_t c) const;
    SA_Range extend_right(uint8_t c) const;
    BranchSet branch(Direction direction) const;
    std::vector<LocatedHit> locate() const;
};
