#!/bin/bash
# P3-1 gate S, B16 powered confirmation (PI ruling 2026-11-04 after p31c: B16 EN_alt failed at test power 0.72).
# Size fixed before running: sigma(EN_alt) ~ 5.46 from p31c; R body 200,000 (40 x 5,000), exact 1,600,000 (32 x 50,000)
# gives se ~0.0130 against a margin of ~0.0498, test power ~0.97. New seeds 3150xxx (R body), 3151xxx (exact).
# Library dev/perf/lib/p31_commit built from 325c664. Pass or fail is this run alone, no pooling with p31_gateS or p31c.
cd ~/research/evolveTrial-merge
LIB=dev/perf/lib/p31_commit; OUT=dev/perf/phase3/p31d_gateS; mkdir -p $OUT
f=B16_R01_unset_thr
{ for i in $(seq 0 31); do s=$((3151001 + i)); echo "$f exact $s 50000 4 $OUT/${f}_exact_$s.rds"; done
  for i in $(seq 0 39); do s=$((3150001 + i)); echo "$f r $s 5000 4 $OUT/${f}_r_$s.rds"; done; } > $OUT/runs.txt
cat $OUT/runs.txt | xargs -P 4 -L 1 sh -c 'test -f "$5" || nice -n 10 Rscript dev/perf/phase3/p31_gateS_run.R '"$LIB"' "$0" "$1" "$2" "$3" "$4" "$5"'
Rscript dev/perf/phase3/p31_gateS_compare.R $OUT dev/perf/phase3/p31d_gateS_report.tsv
