READ=tests/data/0mis/100/reads.fq

/usr/bin/time -f 'tool\torig_batmis_forward\nreal_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.orig_batmis_forward.time.tsv \
  /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/Original/BatMis-3.00/bin/batman \
  -g ../batlib_port/benchmark/dataset_5m_100bp_n4/ref.fa \
  -q "$READ" \
  -o benchmark/orig_batmis.forward.sam \
  -n 0 -mall -c

/usr/bin/time -f 'tool\tfmsearch_forward\nreal_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.fmsearch_forward.time.tsv \
  ./build/fmsearch \
  -i tests/index/ref.forward \
  -q "$READ" \
  -n 0
/usr/bin/time -f 'tool\tfmsearch_forward_loadsa\nreal_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.fmsearch_forwardloadsa.time.tsv \
  ./build/fmsearch \
  -i tests/index/ref.forward \
  -q tests/data/0mis/100/reads.fq \
  -n 0 --load-sa

/usr/bin/time -f 'tool\tfmsearch_forward_locate\nreal_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.fmsearch_forward_locate.time.tsv \
  ./build/fmsearch \
  -i tests/index/ref.forward \
  -q tests/data/0mis/100/reads.fq \
  -n 0 --locate


/usr/bin/time -f 'real_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.search.time.tsv \
  ./build/fmsearch \
  -i tests/index/ref.forward \
  -q "$READ" \
  -n 0 --timing

/usr/bin/time -f 'real_sec\t%e\nuser_sec\t%U\nsys_sec\t%S' \
  -o benchmark/100.sam.time.tsv \
  ./build/fmsearch \
  -i tests/index/ref.forward \
  -q "$READ" \
  -n 0 --sam-output \
  -o benchmark/fmsearch.forward.sam \
  --timing
