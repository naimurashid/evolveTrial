#!/bin/bash
# P3-2 gate S pilot: estimates the spread of each statistic (type I error, power, EN under the null and the alternative)
# and the cost per trial of each engine on the four gate S frames, to size the full run in advance (p32_gateS_size.R).
# Seeds 3220xxx (R body) and 3221xxx (exact), never used by the full run. The pilot is not part of the gate readout.
cd ~/research/evolveTrial-merge
LIB=${LIB:-dev/perf/lib/p32_dev}; OUT=dev/perf/phase3/p32_gateS_pilot; mkdir -p $OUT
{
  i=0
  for f in C01_R01_on_eff C02_R01_off_eff D01_R23_off_current D03_R23_off_symmetric; do
    i=$((i + 1))
    echo "$f exact $((3221000 + i)) 20000 4 $OUT/${f}_exact_$((3221000 + i)).rds"
    echo "$f r $((3220000 + i)) 2000 4 $OUT/${f}_r_$((3220000 + i)).rds"
  done
} > $OUT/runs.txt
cat $OUT/runs.txt | xargs -P 4 -L 1 sh -c 'test -f "$5" || nice -n 10 Rscript dev/perf/phase3/p32_gateS_run.R '"$LIB"' "$0" "$1" "$2" "$3" "$4" "$5"'
Rscript dev/perf/phase3/p32_gateS_compare.R $OUT dev/perf/phase3/p32_gateS_pilot_report.tsv
