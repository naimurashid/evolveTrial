#!/bin/bash
# P3-1 gate S confirmation (PI ruling 2026-11-04): sizes fixed in advance at the extended readout's final sizes, new seeds
# 3140xxx (R body) and 3141xxx (exact), installed 0.0.0.9004 (commit 325c664). Pass/fail is this run alone, no pooling.
cd ~/research/evolveTrial-merge
LIB=$HOME/R/library; OUT=dev/perf/phase3/p31c_gateS; mkdir -p $OUT
jobs() {
  # frame  r_first r_count r_n  exact_first exact_count exact_n   (R body total / exact total)
  echo "B02_R01_off_cut   3140021  8 5000 3141021  8 50000"   # 40,000 / 400,000
  echo "B07_R01_off_drop  3140031  4 5000 3141031  4 50000"   # 20,000 / 200,000
  echo "B15_R01_unset     3140041  2 5000 3141041  2 50000"   # 10,000 / 100,000
  echo "B16_R01_unset_thr 3140111 20 5000 3141111 16 50000"   # 100,000 / 800,000
  echo "B20_R23_off       3140051 16 5000 3141051  8 50000"   # 80,000 / 400,000
}
jobs | while read f rs rc rn es ec en; do
  for i in $(seq 0 $((ec - 1))); do s=$((es + i)); echo "$f exact $s $en 4 $OUT/${f}_exact_$s.rds"; done
  for i in $(seq 0 $((rc - 1))); do s=$((rs + i)); echo "$f r $s $rn 4 $OUT/${f}_r_$s.rds"; done
done > $OUT/runs.txt
cat $OUT/runs.txt | xargs -P 4 -L 1 sh -c 'test -f "$5" || nice -n 10 Rscript dev/perf/phase3/p31_gateS_run.R '"$LIB"' "$0" "$1" "$2" "$3" "$4" "$5"'
Rscript dev/perf/phase3/p31_gateS_compare.R $OUT dev/perf/phase3/p31c_gateS_report.tsv
