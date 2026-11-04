#!/bin/bash
# P3-1 gate S driver: fresh R-body runs (4 FORK workers each, 4 runs at a time) and exact-engine runs, nice'd, then the
# comparison. Seeds 3120xxx (R body) and 3121xxx (exact) are not used anywhere else in dev/perf.
cd ~/research/evolveTrial-merge
LIB=dev/perf/lib/p31_dev; OUT=dev/perf/phase3/p31_gateS; mkdir -p $OUT
jobs() {
  # frame  r_seeds(first count)  r_n  exact_first exact_count exact_n
  echo "B16_R01_unset_thr 3120011 4 5000 3121011 4 50000"
  echo "B02_R01_off_cut   3120021 4 5000 3121021 4 50000"
  echo "B07_R01_off_drop  3120031 4 5000 3121031 4 50000"
  echo "B15_R01_unset     3120041 2 5000 3121041 2 50000"
  echo "B20_R23_off       3120051 8 5000 3121051 4 50000"
}
jobs | while read f rs rc rn es ec en; do
  for i in $(seq 0 $((ec - 1))); do s=$((es + i)); echo "$f exact $s $en 4 $OUT/${f}_exact_$s.rds"; done
  for i in $(seq 0 $((rc - 1))); do s=$((rs + i)); echo "$f r $s $rn 4 $OUT/${f}_r_$s.rds"; done
done > $OUT/runs.txt
cat $OUT/runs.txt | xargs -P 4 -L 1 sh -c 'test -f "$5" || nice -n 10 Rscript dev/perf/phase3/p31_gateS_run.R '"$LIB"' "$0" "$1" "$2" "$3" "$4" "$5"'
Rscript dev/perf/phase3/p31_gateS_compare.R $OUT dev/perf/phase3/p31_gateS_report.tsv
