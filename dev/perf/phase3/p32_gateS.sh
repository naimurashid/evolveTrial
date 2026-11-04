#!/bin/bash
# P3-2 gate S (plan dev/plans/2026-11-04-p3-2-plan.md, Section 5). Exact engine against fresh R-body runs on the four
# bidirectional frames of p32_frames.R, each under its null (scenario 1) and alternative (scenario 2): eight frame x
# scenario cells. Sizes fixed in advance from the pilot (p32_gateS_pilot.sh, p32_gateS_size.R, p32_gateS_sizing.tsv) for
# test power >= 0.95 per statistic at a true difference of 0 under the phase 2 margins; see P3-2_REPORT.md Section 3.
# New seeds 3230xxx (R body) and 3231xxx (exact), not used by the pilot or elsewhere in dev/perf. Pass or fail is this
# run alone: no extension after the readout.
# LIB defaults to the build of the commit under test; set LIB=... to override.
cd ~/research/evolveTrial-merge
LIB=${LIB:-dev/perf/lib/p32_final}; OUT=dev/perf/phase3/p32_gateS; mkdir -p $OUT
jobs() {
  # frame                 r_first r_count r_n  exact_first exact_count exact_n   (R body total / exact total)
  echo "C01_R01_on_eff        3230101 17 5000 3231101 11 50000"   #  85,000 /   550,000
  echo "C02_R01_off_eff       3230201 15 5000 3231201 10 50000"   #  75,000 /   500,000
  echo "D01_R23_off_current   3230301 26 5000 3231301 34 50000"   # 130,000 / 1,700,000
  echo "D03_R23_off_symmetric 3230401 17 5000 3231401 24 50000"   #  85,000 / 1,200,000
}
jobs | while read f rs rc rn es ec en; do
  for i in $(seq 0 $((rc - 1))); do s=$((rs + i)); echo "$f r $s $rn 4 $OUT/${f}_r_$s.rds"; done
  for i in $(seq 0 $((ec - 1))); do s=$((es + i)); echo "$f exact $s $en 4 $OUT/${f}_exact_$s.rds"; done
done > $OUT/runs.txt
cat $OUT/runs.txt | xargs -P 4 -L 1 sh -c 'test -f "$5" || nice -n 10 Rscript dev/perf/phase3/p32_gateS_run.R '"$LIB"' "$0" "$1" "$2" "$3" "$4" "$5"'
Rscript dev/perf/phase3/p32_gateS_compare.R $OUT dev/perf/phase3/p32_gateS_report.tsv
