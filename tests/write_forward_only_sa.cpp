#include "bwt.h"

#include <cstdio>
#include <cstdint>

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::fprintf(stderr, "Usage: write_forward_only_sa <index.bwt> <index.sa>\n");
        return 2;
    }
    bwt_t* bwt = bwt_restore_bwt(argv[1]);
    if (!bwt) return 1;
    bwt_cal_sa(bwt, 16);
    FILE* out = std::fopen(argv[2], "wb");
    if (!out) { bwt_destroy(bwt); return 1; }
    bool ok = std::fwrite(&bwt->primary, sizeof(uint64_t), 1, out) == 1 &&
              std::fwrite(bwt->L2 + 1, sizeof(uint64_t), 4, out) == 4;
    const uint64_t interval = static_cast<uint64_t>(bwt->sa_intv);
    ok = ok && std::fwrite(&interval, sizeof(interval), 1, out) == 1 &&
         std::fwrite(&bwt->seq_len, sizeof(uint64_t), 1, out) == 1 &&
         std::fwrite(bwt->sa + 1, sizeof(uint64_t), bwt->n_sa - 1, out) == bwt->n_sa - 1;
    if (std::fclose(out) != 0) ok = false;
    bwt_destroy(bwt);
    return ok ? 0 : 1;
}
