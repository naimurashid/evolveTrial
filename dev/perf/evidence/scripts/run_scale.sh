#!/bin/bash
cd /tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf
for D in R01 R23; do for sc in 1 2; do for w in 1 4 8 16; do for m in internal external; do
  echo "== $D s$sc w$w $m $(date +%T)"
  R_LIBS=/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf/Rlib timeout 300 Rscript scale.R $D $sc 2400 $w $m 2>&1 | tail -4
done; done; done; done
echo ALLDONE
