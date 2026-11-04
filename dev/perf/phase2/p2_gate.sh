#!/bin/bash
# p2_gate.sh -- gate G5 runner (dev/plans/2026-11-04-phase2-plan.md Section 3). Runs every (design, seed) row of
# acceptance_set_reference_rows.csv (75 rows, 10 designs) through p2_cert.R in mode p2 at R = 10,000 with 19 workers,
# one row at a time, nice'd, resumable (an existing rows/<design>_<seed>.rds is skipped).
#   bash p2_gate.sh [--dry-run]
# Outputs under dev/perf/evidence/phase2_gate/: rows/<design>_<seed>.rds, rows/<design>_<seed>.log, run.log.
# Environment: P2_ENGINE_SHA (defaults to the worktree HEAD; recorded in every rds), P2_LIB (default lib/dev_p2),
# P2_GATE_R (default 10000; for smoke tests only, a value other than 10000 is written to run.log as SMOKE).
# Exit 0 when all 75 rows exist and ran clean in this or an earlier invocation, 1 otherwise.
set -uo pipefail
unset CDPATH; H=$(cd "$(dirname "$0")" >/dev/null && pwd); WT=$(cd "$H/../../.." >/dev/null && pwd)
EV=$WT/dev/perf/evidence; G=$EV/phase2_gate; ROWS=$G/rows; REF=$EV/acceptance_set_reference_rows.csv
R=${P2_GATE_R:-10000}
mkdir -p "$ROWS"
export P2_ENGINE_SHA=${P2_ENGINE_SHA:-$(git -C "$WT" rev-parse HEAD)}
[ -d "${P2_LIB:-$WT/dev/perf/lib/dev_p2}" ] || { echo "p2 library missing: ${P2_LIB:-$WT/dev/perf/lib/dev_p2}" >&2; exit 2; }
exec 9>"$G/.lock"; flock -n 9 || { echo "another p2_gate.sh is running" >&2; exit 2; }
LOG=$G/run.log
log() { echo "$(date +%FT%T) $*" | tee -a "$LOG"; }
[ "$R" = 10000 ] || log "SMOKE: R=$R (not a gate run)"
log "start: engine $P2_ENGINE_SHA, R $R, host $(hostname), load $(cut -d' ' -f1 /proc/loadavg)"
fail=0; done_n=0; skip_n=0
while IFS=, read -r design seed; do
  out=$ROWS/${design}_${seed}.rds
  if [ -f "$out" ]; then skip_n=$((skip_n+1)); continue; fi
  if [ "${1:-}" = "--dry-run" ]; then log "would run $design $seed"; continue; fi
  rm -f "$out.tmp"; t0=$(date +%s)
  nice -n 10 Rscript "$H/p2_cert.R" "$design" "$seed" "$R" p2 "$out.tmp" > "$ROWS/${design}_${seed}.log" 2>&1; st=$?
  if [ $st -eq 0 ] && [ -f "$out.tmp" ]; then
    mv "$out.tmp" "$out"; done_n=$((done_n+1)); log "ok   $design $seed $(( $(date +%s) - t0 )) s: $(grep -h '^\[layer2\] R' "$ROWS/${design}_${seed}.log" | cut -c1-200)"
  else
    fail=$((fail+1)); rm -f "$out.tmp"; log "FAIL $design $seed exit $st after $(( $(date +%s) - t0 )) s (see rows/${design}_${seed}.log)"
  fi
done < <(awk -F, 'NR>1 {print $1","$7}' "$REF")
have=$(ls "$ROWS"/*.rds 2>/dev/null | wc -l)
log "finished: ran $done_n, skipped $skip_n, failed $fail, rds present $have of 75"
[ "$fail" -eq 0 ] && [ "$have" -eq 75 ]
