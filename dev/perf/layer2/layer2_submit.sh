#!/bin/bash
# layer2_submit.sh -- submit layer 2 on Longleaf (plan step 5, dev/plans/2026-11-03-step5-layer2.md, Section 3). Run on
# the login node after layer2_stage.sh, once the PI has the core-hour estimate.
#   bash /proj/rashidlab/nur2/e2_ll/engine_cpp_l2/harness/layer2_submit.sh <tag>
# Submits the 75-task array and one dependent 1-CPU job (layer2_gate.sh) that runs layer2_compare.R on the outputs and
# writes the array's sacct accounting. Results: $L2/out/<tag>/compare.csv, compare.txt (gate line, PASS or FAIL), sacct.tsv, core_h.txt.
set -euo pipefail
. "${L2_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l2/harness/layer2_env.sh}"
TAG=${1:?tag}
[[ "$TAG" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "bad tag" >&2; exit 2; }
[ ! -e "$L2/out/$TAG" ] && [ ! -e "$L2/logs/$TAG" ] || { echo "refusing: tag $TAG used" >&2; exit 2; }
[ -f "$L2_LIB/evolveTrial/DESCRIPTION" ] && grep -qx "EngineCppSHA: $L2_ENGINE_SHA" "$L2_LIB/evolveTrial/DESCRIPTION" ||
  { echo "engine library not staged (run layer2_stage.sh)" >&2; exit 2; }
N=$(($(wc -l < "$L2_REF_ROWS") - 1)); [ "$N" = 75 ] || { echo "reference rows: $N, expected 75" >&2; exit 2; }
mkdir -p "$L2/logs/$TAG"
A=$(sbatch --parsable --array=1-"$N" --output="$L2/logs/$TAG/%A_%a.log" "$L2_HARNESS/layer2.sbatch" "$TAG")
C=$(sbatch --parsable --dependency=afterany:"$A" --job-name=l2cmp --account="$L2_ACCOUNT" --partition="$L2_PARTITION" \
  --cpus-per-task=1 --mem=4G --time=00:15:00 --output="$L2/logs/$TAG/compare_%j.log" "$L2_HARNESS/layer2_gate.sh" "$TAG" "$A")
echo "$A $C $(date -Is)" > "$L2/logs/$TAG/jobs.txt"
echo "submitted array $A (${N} tasks x $L2_CPUS CPUs) and compare job $C, tag $TAG"
