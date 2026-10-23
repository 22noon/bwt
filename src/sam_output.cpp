#include "sam_output.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <tuple>
#include <vector>

namespace {
char complement_base(char base) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(base)))) {
    case 'A': return 'T'; case 'C': return 'G'; case 'G': return 'C'; case 'T': return 'A';
    default: return 'N';
    }
}
std::string reverse_complement(const std::string& seq) {
    std::string out(seq.rbegin(), seq.rend());
    for (char& c : out) c = complement_base(c);
    return out;
}
std::string reverse_string(const std::string& s) { return std::string(s.rbegin(), s.rend()); }
std::string build_md_tag(size_t read_len, const std::vector<Mismatch>& mismatches) {
    if (mismatches.empty()) return std::to_string(read_len);
    auto sorted = mismatches;
    std::sort(sorted.begin(), sorted.end(), [](const Mismatch& a, const Mismatch& b) { return a.Pos < b.Pos; });
    std::string md;
    size_t last = 0;
    for (const auto& mm : sorted) {
        if (mm.Pos == 0 || mm.Pos > read_len) continue;
        const size_t pos = static_cast<size_t>(mm.Pos - 1);
        if (pos < last) continue;
        md += std::to_string(pos - last); // MD needs a zero count before a mismatch at the first/consecutive base.
        char ref = static_cast<char>(std::toupper(static_cast<unsigned char>(mm.Char)));
        if (ref != 'A' && ref != 'C' && ref != 'G' && ref != 'T') ref = 'N';
        md.push_back(ref);
        last = pos + 1;
    }
    md += std::to_string(read_len - last);
    return md;
}
bool write_record(samFile* out, sam_hdr_t* hdr, const std::string& qname,
                  uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
                  const std::vector<uint32_t>& cigar, const std::string& seq,
                  const std::string& qual, const std::vector<Mismatch>* mismatches,
                  size_t nh, size_t hi) {
    bam1_t* rec = bam_init1();
    if (!rec) return false;
    std::string encoded_qual;
    if (!qual.empty()) {
        encoded_qual.reserve(qual.size());
        for (unsigned char q : qual)
            encoded_qual.push_back(static_cast<char>(q >= 33 ? q - 33 : 0));
    }
    const int rc = bam_set1(rec, qname.size(), qname.c_str(), flag, tid, pos, mapq,
                            cigar.size(), cigar.empty() ? nullptr : cigar.data(),
                            -1, -1, 0, seq.size(), seq.c_str(), encoded_qual.empty() ? nullptr : encoded_qual.c_str(), 0);
    if (rc < 0) { bam_destroy1(rec); return false; }
    if (mismatches) {
        const char nm[2] = {'N','M'}, mdtag[2] = {'M','D'}, nhtag[2] = {'N','H'}, hitag[2] = {'H','I'};
        const std::string md = build_md_tag(seq.size(), *mismatches);
        if (bam_aux_update_int(rec, nm, static_cast<int64_t>(mismatches->size())) < 0 ||
            bam_aux_update_str(rec, mdtag, -1, md.c_str()) < 0 ||
            bam_aux_update_int(rec, nhtag, static_cast<int64_t>(nh)) < 0 ||
            bam_aux_update_int(rec, hitag, static_cast<int64_t>(hi)) < 0) {
            bam_destroy1(rec); return false;
        }
    }
    const int wr = sam_write1(out, hdr, rec);
    bam_destroy1(rec);
    return wr >= 0;
}
struct Placement { LocatedHit location; const SA_Range* range; };
} // namespace

bool build_reference_header(sam_hdr_t* header, const BwaFMDIndex* index) {
    if (!header || !index || !index->reference_metadata()) return false;
    const bntseq_t* ref = index->reference_metadata();
    if (sam_hdr_add_line(header, "HD", "VN", "1.6", "SO", "unsorted", NULL) != 0) return false;
    for (int32_t i = 0; i < ref->n_seqs; ++i) {
        const bntann1_t& ann = ref->anns[i];
        if (!ann.name || ann.len < 0) return false;
        const std::string length = std::to_string(ann.len);
        if (sam_hdr_add_line(header, "SQ", "SN", ann.name, "LN", length.c_str(), NULL) != 0) return false;
    }
    return true;
}

bool write_sam_unmapped(samFile* out, sam_hdr_t* hdr, const std::string& name,
                        const std::string& seq, const std::string& qual) {
    return write_record(out, hdr, name, BAM_FUNMAP, -1, -1, 0, {}, seq, qual, nullptr, 0, 0);
}

bool write_sam_hits(samFile* out, sam_hdr_t* hdr, const BwaFMDIndex* index,
                    const std::string& name, const std::string& read_seq,
                    const std::string& read_qual, const std::vector<SA_Range>& hits) {
    if (!out || !hdr || !index) return false;
    std::vector<Placement> placements;
    for (const auto& hit : hits) {
        if (hit.empty()) continue;
        for (const auto& loc : index->locate(hit, read_seq.size())) placements.push_back({loc, &hit});
    }
    std::sort(placements.begin(), placements.end(), [](const Placement& a, const Placement& b) {
        return std::make_tuple(a.location.ref_id, a.location.position, a.location.strand, a.range->mismatch_count()) <
               std::make_tuple(b.location.ref_id, b.location.position, b.location.strand, b.range->mismatch_count());
    });
    placements.erase(std::unique(placements.begin(), placements.end(), [](const Placement& a, const Placement& b) {
        return a.location.ref_id == b.location.ref_id && a.location.position == b.location.position &&
               a.location.strand == b.location.strand && a.range->mismatch_count() == b.range->mismatch_count();
    }), placements.end());
    if (placements.empty()) return write_sam_unmapped(out, hdr, name, read_seq, read_qual);

    const size_t nh = placements.size();
    for (size_t i = 0; i < placements.size(); ++i) {
        const auto& p = placements[i];
        const bool reverse = p.location.strand == Strand::Reverse;
        const std::string seq = reverse ? reverse_complement(read_seq) : read_seq;
        const std::string qual = read_qual.empty() ? std::string() : (reverse ? reverse_string(read_qual) : read_qual);
        std::vector<Mismatch> mm = p.range->mismatches();
        if (reverse) {
            for (auto& m : mm) {
                if (m.Pos > 0 && m.Pos <= read_seq.size()) m.Pos = static_cast<uint32_t>(read_seq.size() - m.Pos + 1);
                // On the reverse-complement half of the FMD index, mismatch chars are in reverse-complement text orientation.
                m.Char = complement_base(m.Char);
            }
        }
        const int tid = sam_hdr_name2tid(hdr, p.location.ref_name.c_str());
        if (tid < 0) continue;
        const std::vector<uint32_t> cigar{bam_cigar_gen(static_cast<uint32_t>(read_seq.size()), BAM_CMATCH)};
        const uint16_t flag = static_cast<uint16_t>((reverse ? BAM_FREVERSE : 0) | (i ? BAM_FSECONDARY : 0));
        const uint8_t mapq = nh == 1 ? 60 : 0;
        if (!write_record(out, hdr, name, flag, tid, static_cast<hts_pos_t>(p.location.position),
                          mapq, cigar, seq, qual, &mm, nh, i + 1)) {
            std::fprintf(stderr, "Failed to write SAM record for read '%s'\n", name.c_str());
            return false;
        }
    }
    return true;
}

bool write_sam_forward_exact_hits(samFile* out, sam_hdr_t* hdr,
                                  const BwaFMDIndex* index,
                                  const std::string& name,
                                  const std::string& read_seq,
                                  const std::string& read_qual,
                                  const std::vector<uint64_t>& positions) {
    if (!out || !hdr || !index || !index->reference_metadata()) return false;
    if (positions.empty()) return write_sam_unmapped(out, hdr, name, read_seq, read_qual);
    const bntseq_t* ref = index->reference_metadata();
    if (ref->n_seqs != 1 || ref->l_pac != static_cast<int64_t>(index->bwt()->seq_len) ||
        ref->anns[0].offset != 0) {
        std::fprintf(stderr, "Forward-only SAM output currently requires one contig spanning the indexed text\n");
        return false;
    }
    const uint64_t ref_length = static_cast<uint64_t>(ref->anns[0].len);
    const uint64_t nh = positions.size();
    const std::vector<uint32_t> cigar{
        bam_cigar_gen(static_cast<uint32_t>(read_seq.size()), BAM_CMATCH)};
    const std::vector<Mismatch> no_mismatches;
    for (size_t i = 0; i < positions.size(); ++i) {
        const uint64_t pos = positions[i];
        if (pos > ref_length || read_seq.size() > ref_length - pos) {
            std::fprintf(stderr, "Located exact-match coordinate is outside the forward reference\n");
            return false;
        }
        const uint16_t flag = i == 0 ? 0 : BAM_FSECONDARY;
        const uint8_t mapq = nh == 1 ? 60 : 0;
        if (!write_record(out, hdr, name, flag, 0,
                          static_cast<hts_pos_t>(pos), mapq, cigar,
                          read_seq, read_qual, &no_mismatches, nh, i + 1)) {
            std::fprintf(stderr, "Failed to write SAM record for read '%s'\n", name.c_str());
            return false;
        }
    }
    return true;
}
