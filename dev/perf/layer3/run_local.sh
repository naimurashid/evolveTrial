#!/bin/bash
# run_local.sh <mode> [tag]  -- local gate runs of the three candidates (sequential, 16 workers)
M=$1; O=dev/perf/layer3/out/local/${2:-$M}; mkdir -p $O
for c in "PB-baton-1001 1" "PB-baton-1003 95" "PB-ga-1001 274"; do set -- $c
  env -i HOME="$HOME" PATH="$PATH" LANG=C.UTF-8 TZ=America/New_York R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null \
    SLURM_CPUS_PER_TASK=16 R_LIBS="$HOME/R/library" L3_ENGINE_SHA=1c918711f48c411fecd1c1dad7e1bdd799171fd9 ${L3_FORCE_R:+L3_FORCE_R=$L3_FORCE_R} \
    Rscript dev/perf/layer3/layer3_search.R $1 $2 $M $O/$1_$2.rds 2>&1 | tee $O/$1_$2.log | cut -c1-400
done
