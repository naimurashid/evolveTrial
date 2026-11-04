#!/bin/bash
S=$(dirname "$0"); P=/home/naimrashid/research/evolveTrial-perf/dev/perf/lib
for D in R01 R23; do for w in 1 16; do for L in ref dev_l2; do
  echo "== $D w$w $L $(date +%T) load $(cut -d' ' -f1 /proc/loadavg)"
  R_LIBS=$P/$L:$HOME/R/library R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null Rscript $S/t7.R $P/$L $D $w $S/t7_${D}_w${w}_${L}.rds 2>&1 | tail -2
done; done; done; echo ALLDONE
