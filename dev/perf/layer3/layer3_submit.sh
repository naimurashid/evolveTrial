#!/bin/bash
# layer3_submit.sh -- submit layer 3 on Longleaf (plan step 6, dev/plans/2026-11-03-step6-layer3.md). Run on the login node
# after layer3_stage.sh, once the PI has the core-hour estimate (about 0.3 core-h, at most 0.5).
#   bash /proj/rashidlab/nur2/e2_ll/engine_cpp_l3/harness/layer3_submit.sh <tag>
# Submits the 3-task array and one dependent 1-CPU job (layer3_gate.sh) that runs layer3_compare.R on the outputs and writes
# the array's sacct accounting. Results: $L3/out/<tag>/compare.csv, compare.txt (gate line, PASS or FAIL), sacct.tsv, core_h.txt.
set -euo pipefail
. "${L3_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l3/harness/layer3_env.sh}"
TAG=${1:?tag}
[[ "$TAG" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "bad tag" >&2; exit 2; }
[ ! -e "$L3/out/$TAG" ] && [ ! -e "$L3/logs/$TAG" ] || { echo "refusing: tag $TAG used" >&2; exit 2; }
[ -f "$L3_LIB/evolveTrial/DESCRIPTION" ] && grep -qx "EngineCppSHA: $L3_ENGINE_SHA" "$L3_LIB/evolveTrial/DESCRIPTION" ||
  { echo "engine library $L3_LIB is not engine $L3_ENGINE_SHA" >&2; exit 2; }
N=${#L3_CANDS[@]}; [ "$N" = 3 ] || { echo "candidates: $N, expected 3" >&2; exit 2; }
for c in "${L3_CANDS[@]}"; do read -r r e <<< "$c"; [ -f "$L3_RESULTS/$r/sim_ledger.csv" ] && [ -f "$L3_RESULTS/$(l3_rds_path "$r")" ] ||
  { echo "results not staged for $r (run layer3_stage.sh)" >&2; exit 2; }; done
mkdir -p "$L3/logs/$TAG"
A=$(sbatch --parsable --array=1-"$N" --output="$L3/logs/$TAG/%A_%a.log" "$L3_HARNESS/layer3.sbatch" "$TAG")
C=$(sbatch --parsable --dependency=afterany:"$A" --job-name=l3cmp --account="$L3_ACCOUNT" --partition="$L3_PARTITION" \
  --cpus-per-task=1 --mem=4G --time=00:15:00 --output="$L3/logs/$TAG/compare_%j.log" "$L3_HARNESS/layer3_gate.sh" "$TAG" "$A")
echo "$A $C $(date -Is)" > "$L3/logs/$TAG/jobs.txt"
echo "submitted array $A (${N} tasks x $L3_CPUS CPUs) and compare job $C, tag $TAG"
