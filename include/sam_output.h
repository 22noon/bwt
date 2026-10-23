#ifndef FMSEARCH_SAM_OUTPUT_H
#define FMSEARCH_SAM_OUTPUT_H
#include "bwa_fmd_index.hpp"
#include <string>
#include <vector>
#include <htslib/sam.h>

bool build_reference_header(sam_hdr_t* header, const BwaFMDIndex* index);
bool write_sam_hits(samFile* out, sam_hdr_t* header, const BwaFMDIndex* index,
                    const std::string& read_name, const std::string& read_sequence,
                    const std::string& read_quality, const std::vector<SA_Range>& hits);
bool write_sam_forward_exact_hits(samFile* out, sam_hdr_t* header,
                                  const BwaFMDIndex* index,
                                  const std::string& read_name,
                                  const std::string& read_sequence,
                                  const std::string& read_quality,
                                  const std::vector<uint64_t>& positions);
bool write_sam_unmapped(samFile* out, sam_hdr_t* header, const std::string& read_name,
                        const std::string& read_sequence, const std::string& read_quality);
#endif
