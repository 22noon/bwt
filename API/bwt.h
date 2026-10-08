/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef API_BWT_H
#define API_BWT_H

#include <stdint.h>
#include <stddef.h>

#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#ifndef BWA_UBYTE
#define BWA_UBYTE
typedef unsigned char ubyte_t;
#endif

typedef uint64_t bwtint_t;

typedef struct {
	bwtint_t primary;
	bwtint_t L2[5];
	bwtint_t seq_len;
	bwtint_t bwt_size;
	uint32_t *bwt;
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	bwtint_t x[3], info;
} bwtintv_t;

typedef struct { size_t n, m; bwtintv_t *a; } bwtintv_v;

#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))

#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

#define bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)

#ifdef __cplusplus
extern "C" {
#endif

	void bwt_gen_cnt_table(bwt_t *bwt);
	bwt_t *bwt_restore_bwt(const char *fn);
	void bwt_restore_sa(const char *fn, bwt_t *bwt);
	void bwt_destroy(bwt_t *bwt);
	void bwt_cal_sa(bwt_t *bwt, int intv);

	bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c);
	void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4]);
	bwtint_t bwt_sa(const bwt_t *bwt, bwtint_t k);

	void bwt_2occ(const bwt_t *bwt, bwtint_t k, bwtint_t l, ubyte_t c, bwtint_t *ok, bwtint_t *ol);
	void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4]);

	void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back);

#ifdef __cplusplus
}
#endif

#endif
