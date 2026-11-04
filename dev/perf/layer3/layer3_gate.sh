#!/bin/bash
# layer3_gate.sh -- the dependent 1-CPU job of layer3_submit.sh (plan step 6, dev/plans/2026-11-03-step6-layer3.md).
#   sbatch --dependency=afterany:<array> ... layer3_gate.sh <tag> <array job id>
# Runs layer3_compare.R on $L3/out/<tag> against the staged results and writes the array's sacct accounting. Outputs in
# $L3/out/<tag>: compare.csv, compare.txt (gate line and "compare exit <status>"), sacct.tsv, core_h.txt.
set -uo pipefail
. "${L3_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l3/harness/layer3_env.sh}"
TAG=${1:?tag}; A=${2:?array job id}; O=$L3/out/$TAG; mkdir -p "$O"
l3_modules || exit 2
R_LIBS="$L3_DEPS_LIB" L3_ENGINE_SHA=$L3_ENGINE_SHA Rscript "$L3_HARNESS/layer3_compare.R" "$O" "$L3_RESULTS" "$O/compare.csv" dev > "$O/compare.txt" 2>&1
RC=$?; echo "compare exit $RC" >> "$O/compare.txt"
sacct -j "$A" -X -n -P -o JobID,Elapsed,AllocCPUS,State,NodeList > "$O/sacct.tsv"
awk -F'|' '{ split($2, t, ":"); c += (t[1] * 3600 + t[2] * 60 + t[3]) * $3; n++ } END { printf "%d tasks, %.2f core-h\n", n, c / 3600 }' \
  "$O/sacct.tsv" > "$O/core_h.txt"
cat "$O/compare.txt" "$O/core_h.txt"
exit $RC   # a FAIL shows as a failed gate job in sacct
