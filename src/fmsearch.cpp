#include "mismatch_pipeline.hpp"
#include "sam_output.h"
#include "exact_search.hpp"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <unistd.h>
#include <zlib.h>
#include <htslib/kseq.h>

extern "C" {
#include "bwt.h"
#include "bntseq.h"
}
KSEQ_INIT(gzFile, gzread)

namespace {
struct Options {
    std::string index_prefix, reads_path, output_path = "-";
    int max_mismatches = 1;
    ZeroMismatchMode zero_mode = ZeroMismatchMode::Fast;
    bool all_hits = false, progress = true, load_sa = false, locate = false;
    bool sam_output = false, timing = false, seeded_6mer = false;
    bool batmis_singleton_fast_path = false;
    RankBackend rank_backend = RankBackend::Bwa;
};
struct Progress {
    bool enabled = false;
    uint64_t total_bytes = 0, target_bytes = 0, reads_seen = 0;
    uint64_t reads_since_emit = 0, last_emit_bytes = 0;
};
constexpr uint64_t kWarmupReads = 200, kWarmupEmitEvery = 50;
constexpr uint64_t kMinTargetBytes = 256ULL * 1024ULL, kMinBatchReads = 50, kMaxBatchReads = 1000;
uint64_t clamp_u64(uint64_t v, uint64_t lo, uint64_t hi) { return std::max(lo, std::min(v, hi)); }
uint64_t adaptive_batch(const Progress& p, uint64_t consumed) {
    if (p.reads_seen <= kWarmupReads) return kWarmupEmitEvery;
    const uint64_t avg = std::max<uint64_t>(1, consumed / p.reads_seen);
    return clamp_u64((p.target_bytes + avg / 2) / avg, kMinBatchReads, kMaxBatchReads);
}
void show_progress(uint64_t consumed, uint64_t total, uint64_t reads) {
    if (!total) return;
    consumed = std::min(consumed, total);
    constexpr int width = 30;
    const double ratio = static_cast<double>(consumed) / static_cast<double>(total);
    const int filled = std::max(0, std::min(width, static_cast<int>(ratio * width)));
    char bar[width + 1];
    for (int i = 0; i < width; ++i) bar[i] = i < filled ? '#' : '-';
    bar[width] = '\0';
    std::fprintf(stderr, "\rAlign progress [%s] %6.2f%% (%llu/%llu bytes, %llu reads)", bar,
        ratio * 100.0, static_cast<unsigned long long>(consumed),
        static_cast<unsigned long long>(total), static_cast<unsigned long long>(reads));
    std::fflush(stderr);
}
void usage(const char* program) {
    std::fprintf(stderr,
        "Usage: %s [align] -i <bwa-index-prefix> -q <reads.fa[.gz]|reads.fq[.gz]> [-o output.sam --sam-output] [-n 0|1] [-z fast|bidirectional] [--all-hits] [--load-sa|--locate] [--seeded-6mer] [--rank-backend bwa|batmis] [--batmis-singleton-fastpath] [--timing] [--no-progress]\n"
        "  -i PREFIX  Existing BWA index prefix (.bwt/.sa/.ann/.amb/.pac)\n"
        "  -q FILE    FASTA/FASTQ, plain or gzip; - reads stdin\n"
        "  -o FILE    SAM output path (used with --sam-output; default stdout)\n"
        "  -n INT     Maximum mismatches 0 or 1 (default 1)\n"
        "  -z MODE    fast or bidirectional zero-mismatch mode\n"
        "  --all-hits Also report one-mismatch hits when exact hits exist\n"
        "  --load-sa  Load sampled SA without locating hits\n"
        "  --locate   Load SA and locate all exact-match hits (forward-only index)\n"
        "  --sam-output  Write exact-match SAM records; implies --locate and requires -n 0\n"
        "  --seeded-6mer Opt in to 6-mer seeded exact search (requires -n 0)\n"
        "  --rank-backend Select rank backend: bwa (default) or batmis (requires -n 0)\n"
        "  --batmis-singleton-fastpath Use the experimental BatMIS singleton extension path\n"
        "  --timing   Report search, locating, SAM-output and total elapsed times on stderr\n"
        "  --no-progress Disable progress bar (stderr)\n", program);
}
bool parse_n(const char* s, int& out) {
    char* end = nullptr; errno = 0;
    long v = std::strtol(s, &end, 10);
    if (errno || !end || *end || v < 0 || v > 1) return false;
    out = static_cast<int>(v); return true;
}
bool parse_options(int argc, char** argv, Options& o) {
    int i = 1;
    if (i < argc && std::string(argv[i]) == "align") ++i;
    while (i < argc) {
        const std::string arg = argv[i++];
        auto value = [&]() -> const char* {
            if (i >= argc) { std::fprintf(stderr, "Missing value for %s\n", arg.c_str()); return nullptr; }
            return argv[i++];
        };
        if (arg == "-h" || arg == "--help") { usage(argv[0]); std::exit(0); }
        if (arg == "-i" || arg == "--index") { const char* v=value(); if(!v)return false; o.index_prefix=v; }
        else if (arg == "-q" || arg == "--reads") { const char* v=value(); if(!v)return false; o.reads_path=v; }
        else if (arg == "-o" || arg == "--output") { const char* v=value(); if(!v)return false; o.output_path=v; }
        else if (arg == "-n" || arg == "--max-mismatches") { const char* v=value(); if(!v || !parse_n(v,o.max_mismatches)) { std::fprintf(stderr,"-n currently accepts only 0 or 1\n"); return false; } }
        else if (arg == "-z" || arg == "--zero-mode") {
            const char* v=value(); if(!v)return false; const std::string mode=v;
            if (mode == "fast" || mode == "Fast") o.zero_mode=ZeroMismatchMode::Fast;
            else if (mode == "bidirectional" || mode == "bi" || mode == "Bidirectional" || mode == "BiDirectional") o.zero_mode=ZeroMismatchMode::Bidirectional;
            else { std::fprintf(stderr,"Unknown zero mode: %s\n",v); return false; }
        } else if (arg == "--all-hits") o.all_hits=true;
        else if (arg == "--load-sa") o.load_sa=true;
        else if (arg == "--locate") { o.locate=true; o.load_sa=true; }
        else if (arg == "--sam-output" || arg == "--sam") { o.sam_output=true; o.locate=true; o.load_sa=true; }
        else if (arg == "--seeded-6mer") o.seeded_6mer=true;
        else if (arg == "--batmis-singleton-fastpath") o.batmis_singleton_fast_path=true;
        else if (arg == "--rank-backend") {
            const char* v=value(); if(!v)return false; const std::string backend=v;
            if (backend == "bwa") o.rank_backend=RankBackend::Bwa;
            else if (backend == "batmis") o.rank_backend=RankBackend::BatMIS;
            else { std::fprintf(stderr,"Unknown rank backend: %s\n",v); return false; }
        }
        else if (arg == "--timing") o.timing=true;
        else if (arg == "--no-progress") o.progress=false;
        else { std::fprintf(stderr,"Unknown option: %s\n",arg.c_str()); return false; }
    }
    if (o.index_prefix.empty() || o.reads_path.empty()) { usage(argv[0]); return false; }
    if (o.sam_output && o.max_mismatches != 0) {
        std::fprintf(stderr, "--sam-output currently supports exact matches only; specify -n 0\n");
        return false;
    }
    if (o.rank_backend == RankBackend::BatMIS && o.max_mismatches != 0) {
        std::fprintf(stderr, "--rank-backend batmis currently supports exact search only; specify -n 0\n");
        return false;
    }
    if (o.batmis_singleton_fast_path && o.rank_backend != RankBackend::BatMIS) {
        std::fprintf(stderr, "--batmis-singleton-fastpath requires --rank-backend batmis\n");
        return false;
    }
    if (o.seeded_6mer && o.max_mismatches != 0) {
        std::fprintf(stderr, "--seeded-6mer currently supports exact search only; specify -n 0\n");
        return false;
    }
    return true;
}
uint64_t file_size(const std::string& path) {
    if (path == "-") return 0;
    std::error_code ec; auto n = std::filesystem::file_size(path, ec);
    return ec ? 0 : static_cast<uint64_t>(n);
}
bool load_sa_checked(const std::string& path, bwt_t* bwt, std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) { error = "Cannot open sampled suffix array: " + path; return false; }
    uint64_t header[7]{};
    in.read(reinterpret_cast<char*>(header), sizeof(header));
    if (!in) { error = "Truncated or invalid sampled suffix array: " + path; return false; }
    if (header[0] != bwt->primary || header[6] != bwt->seq_len) {
        error = "Sampled suffix array is incompatible with BWT (primary/sequence length mismatch): " + path;
        return false;
    }
    for (int i = 0; i < 4; ++i) {
        if (header[i + 1] != bwt->L2[i + 1]) {
            error = "Sampled suffix array is incompatible with BWT (cumulative counts mismatch): " + path;
            return false;
        }
    }
    const uint64_t intv = header[5];
    if (intv == 0 || intv > static_cast<uint64_t>(std::numeric_limits<int>::max()) || (intv & (intv - 1))) {
        error = "Invalid sampled-SA interval in " + path;
        return false;
    }
    const uint64_t samples = (bwt->seq_len + intv) / intv;
    std::error_code ec;
    const uint64_t bytes = std::filesystem::file_size(path, ec);
    if (ec || samples == 0 || bytes != sizeof(header) + (samples - 1) * sizeof(uint64_t)) {
        error = "Sampled suffix array has an invalid size or is truncated: " + path;
        return false;
    }
    bwt_restore_sa(path.c_str(), bwt);
    if (!bwt->sa || bwt->sa_intv != static_cast<int>(intv)) {
        error = "Failed to restore sampled suffix array: " + path;
        return false;
    }
    for (uint64_t i = 1; i < bwt->n_sa; ++i) {
        if (bwt->sa[i] > bwt->seq_len) {
            error = "Sampled suffix array contains an out-of-range coordinate: " + path;
            return false;
        }
    }
    return true;
}
gzFile open_reads(const std::string& path) {
    if (path != "-") return gzopen(path.c_str(), "r");
    int fd = dup(::fileno(stdin));
    if (fd < 0) return nullptr;
    gzFile fp = gzdopen(fd, "r");
    if (!fp) close(fd);
    return fp;
}
bool acgt_only(const std::string& s) {
    if (s.empty()) return false;
    for (unsigned char c : s) {
        switch (static_cast<char>(std::toupper(c))) {
        case 'A': case 'C': case 'G': case 'T': break;
        default: return false;
        }
    }
    return true;
}
} // namespace

int main(int argc, char** argv) {
    Options options;
    if (!parse_options(argc, argv, options)) return 2;
    const auto total_start = std::chrono::steady_clock::now();
    std::chrono::duration<double> bwt_load_elapsed{0}, reference_load_elapsed{0};
    std::chrono::duration<double> input_elapsed{0}, sa_load_elapsed{0}, sam_output_elapsed{0};
    std::chrono::duration<double> seed_table_elapsed{0}, rank_backend_setup_elapsed{0};
    std::chrono::duration<double> search_elapsed{0}, locate_elapsed{0}, progress_elapsed{0};
    const std::string bwt_path = options.index_prefix + ".bwt";
    const std::string sa_path = options.index_prefix + ".sa";
    const auto bwt_load_start = std::chrono::steady_clock::now();
    bwt_t* bwt = bwt_restore_bwt(bwt_path.c_str());
    bwt_load_elapsed = std::chrono::steady_clock::now() - bwt_load_start;
    if (!bwt) { std::fprintf(stderr,"Failed to load BWT index: %s\n",bwt_path.c_str()); return 1; }
    if (options.load_sa) {
        const auto sa_start = std::chrono::steady_clock::now();
        std::string error;
        if (!load_sa_checked(sa_path, bwt, error)) {
            std::fprintf(stderr, "%s\n", error.c_str());
            bwt_destroy(bwt);
            return 1;
        }
        sa_load_elapsed += std::chrono::steady_clock::now() - sa_start;
    }
    const auto reference_load_start = std::chrono::steady_clock::now();
    bntseq_t* reference = bns_restore(options.index_prefix.c_str());
    reference_load_elapsed = std::chrono::steady_clock::now() - reference_load_start;
    if (!reference) {
        std::fprintf(stderr,"Failed to load reference metadata for prefix '%s' (.ann/.amb/.pac required)\n",options.index_prefix.c_str());
        bwt_destroy(bwt); return 1;
    }
    if (options.sam_output &&
        (reference->n_seqs != 1 || reference->n_holes != 0 ||
         reference->l_pac != static_cast<int64_t>(bwt->seq_len) ||
         reference->anns[0].offset != 0 ||
         reference->anns[0].len < 0 ||
         static_cast<uint64_t>(reference->anns[0].len) != bwt->seq_len)) {
        std::fprintf(stderr, "--sam-output currently requires a single-contig A/C/G/T forward-only index\n");
        bns_destroy(reference); bwt_destroy(bwt); return 1;
    }
    const auto backend_setup_start = std::chrono::steady_clock::now();
    BwaFMDIndex index(bwt, reference, options.rank_backend,
                      options.batmis_singleton_fast_path);
    rank_backend_setup_elapsed = std::chrono::steady_clock::now() - backend_setup_start;
    std::unique_ptr<ExactSearch6merTable> seeded_search;
    if (options.seeded_6mer) {
        const auto table_start = std::chrono::steady_clock::now();
        seeded_search = std::make_unique<ExactSearch6merTable>(index);
        seed_table_elapsed = std::chrono::steady_clock::now() - table_start;
    }
    samFile* output = nullptr;
    sam_hdr_t* header = nullptr;
    if (options.sam_output) {
        const auto output_start = std::chrono::steady_clock::now();
        output = sam_open(options.output_path.empty() || options.output_path == "-" ? "-" : options.output_path.c_str(), "w");
        if (!output) {
            std::fprintf(stderr,"Could not open SAM output: %s\n",options.output_path.c_str()); bns_destroy(reference); bwt_destroy(bwt); return 1;
        }
        header = sam_hdr_init();
        if (!header || !build_reference_header(header, &index) || sam_hdr_write(output, header) < 0) {
            std::fprintf(stderr,"Failed to write SAM header\n"); if(header)sam_hdr_destroy(header); sam_close(output); bns_destroy(reference); bwt_destroy(bwt); return 1;
        }
        sam_output_elapsed += std::chrono::steady_clock::now() - output_start;
    }
    gzFile input = open_reads(options.reads_path);
    if (!input) {
        std::fprintf(stderr,"Could not open reads: %s\n",options.reads_path.c_str()); if(header)sam_hdr_destroy(header); if(output)sam_close(output); bns_destroy(reference); bwt_destroy(bwt); return 1;
    }
    kseq_t* seq = kseq_init(input);
    if (!seq) {
        std::fprintf(stderr,"Could not initialize kseq parser\n"); gzclose(input); if(header)sam_hdr_destroy(header); if(output)sam_close(output); bns_destroy(reference); bwt_destroy(bwt); return 1;
    }
    Progress progress;
    progress.total_bytes = file_size(options.reads_path);
    progress.enabled = options.progress && progress.total_bytes > 0;
    progress.target_bytes = std::max<uint64_t>(kMinTargetBytes, progress.total_bytes / 200);
    int read_status = 0, exit_status = 0;
    uint64_t located_hits = 0, exact_hits = 0;
    while (true) {
        const auto input_start = options.timing
            ? std::chrono::steady_clock::now()
            : std::chrono::steady_clock::time_point{};
        read_status = kseq_read(seq);
        if (read_status < 0) {
            if (options.timing) input_elapsed += std::chrono::steady_clock::now() - input_start;
            break;
        }
        const std::string name(seq->name.s ? seq->name.s : "", seq->name.l);
        const std::string input_sequence(seq->seq.s ? seq->seq.s : "", seq->seq.l);
        std::string read = input_sequence;
        std::string qual;
        if (seq->qual.l == seq->seq.l && seq->qual.l > 0) qual.assign(seq->qual.s, seq->qual.l);
        for (char& c : read) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (options.timing) input_elapsed += std::chrono::steady_clock::now() - input_start;
        bool ok = false;
        if (!acgt_only(read)) {
            if (options.sam_output) {
                const auto output_start = std::chrono::steady_clock::now();
                ok = write_sam_unmapped(output, header, name, input_sequence, qual);
                sam_output_elapsed += std::chrono::steady_clock::now() - output_start;
            } else ok = true;
        } else {
            try {
                SA_Range exact;
                std::vector<uint64_t> positions;
                if (options.sam_output || options.locate || options.max_mismatches == 0) {
                    const auto search_start = std::chrono::steady_clock::now();
                    exact = seeded_search
                        ? seeded_search->search(index, read)
                        : search_exact_unidirectional(index, read, Direction::Left);
                    search_elapsed += std::chrono::steady_clock::now() - search_start;
                    exact_hits += exact.size();
                } else {
                    MismatchPipelineState state;
                    const auto search_start = std::chrono::steady_clock::now();
                    Run_Mismatch_Pipeline(read, &index, state, options.max_mismatches, options.all_hits, options.zero_mode);
                    search_elapsed += std::chrono::steady_clock::now() - search_start;
                    if (options.max_mismatches == 0)
                        for (const auto& hit : state.hits) exact_hits += hit.size();
                }
                if (options.locate) {
                    const auto one_locate_start = std::chrono::steady_clock::now();
                    positions = index.locate_forward(exact);
                    located_hits += positions.size();
                    locate_elapsed += std::chrono::steady_clock::now() - one_locate_start;
                }
                if (options.sam_output) {
                    const auto output_start = std::chrono::steady_clock::now();
                    ok = write_sam_forward_exact_hits(output, header, &index, name,
                                                      input_sequence, qual, positions);
                    sam_output_elapsed += std::chrono::steady_clock::now() - output_start;
                } else ok = true;
            } catch (const std::exception& e) {
                std::fprintf(stderr,"\nRead '%s' failed: %s\n",name.c_str(),e.what());
                ok = false;
            }
        }
        if (!ok) { exit_status=1; break; }
        if (progress.enabled) {
            const auto progress_start = options.timing
                ? std::chrono::steady_clock::now()
                : std::chrono::steady_clock::time_point{};
            ++progress.reads_seen; ++progress.reads_since_emit;
            const z_off_t offset = gzoffset(input);
            if (offset >= 0) {
                const uint64_t consumed = std::min<uint64_t>(static_cast<uint64_t>(offset),progress.total_bytes);
                const uint64_t batch = adaptive_batch(progress,consumed);
                bool emit = progress.reads_since_emit >= batch;
                if (!emit && consumed >= progress.last_emit_bytes) emit = consumed-progress.last_emit_bytes >= progress.target_bytes;
                if (emit) { show_progress(consumed,progress.total_bytes,progress.reads_seen); progress.reads_since_emit=0; progress.last_emit_bytes=consumed; }
            }
            if (options.timing) progress_elapsed += std::chrono::steady_clock::now() - progress_start;
        }
    }
    if (read_status == -2) { std::fprintf(stderr,"\nMalformed FASTQ: truncated or mismatched quality string\n"); exit_status=1; }
    if (progress.enabled) {
        const auto progress_start = options.timing
            ? std::chrono::steady_clock::now()
            : std::chrono::steady_clock::time_point{};
        show_progress(progress.total_bytes,progress.total_bytes,progress.reads_seen);
        std::fprintf(stderr,"\n");
        if (options.timing) progress_elapsed += std::chrono::steady_clock::now() - progress_start;
    }
    kseq_destroy(seq); gzclose(input); if(header)sam_hdr_destroy(header);
    if (output) {
        const auto output_close_start = std::chrono::steady_clock::now();
        if (sam_close(output) < 0) exit_status=1;
        sam_output_elapsed += std::chrono::steady_clock::now() - output_close_start;
    }
    bns_destroy(reference); bwt_destroy(bwt);
    const auto total_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start);
    if (options.max_mismatches == 0)
        std::fprintf(stderr, "fmsearch: exact hits=%llu located=%llu\n",
                     static_cast<unsigned long long>(exact_hits),
                     static_cast<unsigned long long>(located_hits));
    if (options.timing) {
        const double classified = bwt_load_elapsed.count() + reference_load_elapsed.count() +
            input_elapsed.count() + search_elapsed.count() +
            seed_table_elapsed.count() + rank_backend_setup_elapsed.count() + sa_load_elapsed.count() + locate_elapsed.count() + sam_output_elapsed.count() +
            progress_elapsed.count();
        const double other = std::max(0.0, total_elapsed.count() - classified);
        std::fprintf(stderr, "fmsearch timing (steady_clock): bwt_load=%.6fs reference_load=%.6fs input_parse=%.6fs search=%.6fs seed_table_build=%.6fs rank_backend_setup=%.6fs sa_load=%.6fs locate=%.6fs sam_output=%.6fs progress=%.6fs other=%.6fs total=%.6fs\n",
                     bwt_load_elapsed.count(), reference_load_elapsed.count(), input_elapsed.count(),
                     search_elapsed.count(), seed_table_elapsed.count(), rank_backend_setup_elapsed.count(), sa_load_elapsed.count(),
                     locate_elapsed.count(), sam_output_elapsed.count(), progress_elapsed.count(),
                     other, total_elapsed.count());
    }
    return exit_status;
}
