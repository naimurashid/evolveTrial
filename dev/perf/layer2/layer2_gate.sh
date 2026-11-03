#!/bin/bash
# layer2_gate.sh -- the dependent 1-CPU job of layer2_submit.sh (plan step 5, dev/plans/2026-11-03-step5-layer2.md).
#   sbatch --dependency=afterany:<array> ... layer2_gate.sh <tag> <array job id>
# Runs layer2_compare.R on $L2/out/<tag> and writes the array's sacct accounting. Outputs in $L2/out/<tag>: compare.csv,
# compare.txt (gate line and "compare exit <status>"), sacct.tsv, core_h.txt.
set -uo pipefail
. "${L2_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l2/harness/layer2_env.sh}"
TAG=${1:?tag}; A=${2:?array job id}; O=$L2/out/$TAG; mkdir -p "$O"
l2_modules || exit 2
L2_ENGINE_SHA=$L2_ENGINE_SHA Rscript "$L2_HARNESS/layer2_compare.R" "$O" "$L2_REF_ROWS" "$O/compare.csv" > "$O/compare.txt" 2>&1
echo "compare exit $?" >> "$O/compare.txt"
sacct -j "$A" -X -n -P -o JobID,Elapsed,AllocCPUS,State,NodeList > "$O/sacct.tsv"
awk -F'|' '{ split($2, t, ":"); c += (t[1] * 3600 + t[2] * 60 + t[3]) * $3; n++ } END { printf "%d tasks, %.2f core-h\n", n, c / 3600 }' \
  "$O/sacct.tsv" > "$O/core_h.txt"
cat "$O/compare.txt" "$O/core_h.txt"
