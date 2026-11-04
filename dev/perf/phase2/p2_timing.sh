#!/bin/bash
# p2_timing.sh -- gate G6 (dev/plans/2026-11-04-phase2-plan.md Section 3). Times R01 and R23 at 1 and 16 workers for the
# phase 1 library (dev_l2, engine option "cpp") and the phase 2 library (dev_p2, engine "exact"), 1,000 trials each,
# production path, as in evidence/step7_timing/scripts. Run it on a quiet machine (the script logs the load average before
# each run and refuses to start above P2_MAX_LOAD, default 1.0). Output CSV: evidence/phase2_gate/timing/timing_p2.csv
# (columns design, lib, engine, workers, elapsed, per_1000_s, cpu, digest, sha, load1, rep); per-run rds beside it.
#   bash p2_timing.sh [reps]     reps default 1 (use 3 for a gate-quality table; rows are kept, medians are in the report)
set -uo pipefail
unset CDPATH; H=$(cd "$(dirname "$0")" >/dev/null && pwd); WT=$(cd "$H/../../.." >/dev/null && pwd); P=$WT/dev/perf/lib
O=$WT/dev/perf/evidence/phase2_gate/timing; mkdir -p "$O"; REPS=${1:-1}; MAXLOAD=${P2_MAX_LOAD:-1.0}
PL=${P2_TIMING_LIB:-dev_p2}; [ -d "$P/$PL" ] || { echo "$PL missing" >&2; exit 2; }
load=$(cut -d' ' -f1 /proc/loadavg); awk -v l="$load" -v m="$MAXLOAD" 'BEGIN{exit !(l>m)}' && { echo "load $load > $MAXLOAD; machine not quiet" >&2; exit 2; }
for rep in $(seq 1 "$REPS"); do for D in R01 R23; do for w in 1 16; do for LE in dev_l2:cpp $PL:exact; do
  L=${LE%:*}; E=${LE#*:}; f=$O/t_${D}_w${w}_${L}_rep${rep}.rds
  [ -f "$f" ] && continue
  echo "== $D w$w $L $E rep$rep $(date +%T) load $(cut -d' ' -f1 /proc/loadavg)"
  R_LIBS=$P/$L:$HOME/R/library R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null Rscript "$H/p2_t7.R" "$P/$L" "$D" "$w" "$E" "$f" 2>&1 | tail -2
done; done; done; done
Rscript -e 'f <- list.files("'"$O"'", "^t_.*rds$", full.names = TRUE)
  d <- do.call(rbind, lapply(f, function(z) { r <- readRDS(z); r$rep <- as.integer(sub(".*_rep([0-9]+)\\.rds$", "\\1", z)); r }))
  d$per_1000_s <- d$elapsed; d <- d[order(d$design, d$workers, d$lib, d$rep), ]; write.csv(d, "'"$O"'/timing_p2.csv", row.names = FALSE)
  m <- aggregate(per_1000_s ~ design + workers + lib, d, median); w <- reshape(m, idvar = c("design", "lib"), timevar = "workers", direction = "wide")
  w$efficiency16 <- w$per_1000_s.1 / (16 * w$per_1000_s.16); print(w, row.names = FALSE)'
echo "(efficiency16 = T(1 worker) / (16 x T(16 workers)); G6 target: R23 < 8 s per 1,000 trials at one worker for dev_p2, and efficiency16 above dev_l2.)"
