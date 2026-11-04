#!/bin/bash
#SBATCH --job-name=l3cpp
#SBATCH --account=rc_nur2_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --hint=nomultithread
#SBATCH --constraint=amd
# layer3.sbatch -- one array task per candidate (plan step 6, dev/plans/2026-11-03-step6-layer3.md).
#   sbatch --array=1-3 --output=$L3/logs/<tag>/%A_%a.log layer3.sbatch <tag>
# Task i runs candidate i of L3_CANDS through layer3_search.R with the dev engine at R = 3,000 and 16 workers. The run block
# copies layer2.sbatch (clean environment, fixed FORK cluster port from a per-host registry, R_ENVIRON_USER and
# R_PROFILE_USER disabled). The search script sets the E2c environment block itself. Writes only under $L3.
set -euo pipefail
. "${L3_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l3/harness/layer3_env.sh}"
TAG=${1:?tag}; I=${SLURM_ARRAY_TASK_ID:?array task}; JID=${SLURM_JOB_ID:?job}; CPUS=${SLURM_CPUS_PER_TASK:?cpus}
[ "$CPUS" = "$L3_CPUS" ] || { echo "cpus $CPUS != $L3_CPUS" >&2; exit 2; }
[ "$I" -ge 1 ] && [ "$I" -le "${#L3_CANDS[@]}" ] || { echo "no candidate $I" >&2; exit 2; }
read -r RUN EVI <<< "${L3_CANDS[$((I - 1))]}"
grep -qx "EngineCppSHA: $L3_ENGINE_SHA" "$L3_LIB/evolveTrial/DESCRIPTION" || { echo "library is not engine $L3_ENGINE_SHA" >&2; exit 2; }
OUTD=$L3/out/$TAG; mkdir -p "$OUTD" "$L3/ports"; OUTF=$OUTD/${RUN}_${EVI}.rds
[ ! -e "$OUTF" ] || { echo "refusing: $OUTF exists" >&2; exit 2; }
l3_modules
TMPD=$(mktemp -d "/tmp/l3cpp_${JID}.XXXXXX")
HOST=$(hostname -s); REG=$L3/ports/$HOST.reg; touch "$REG"; PORT=""
exec 9>"$REG.lock"; flock 9
for ((k = 0; k < 500; k++)); do p=$((L3_PORTS + (JID + I + k) % 500))
  if [ -z "$(ss -Htan "( sport = :$p )" 2>/dev/null)" ] && ! awk -v p="$p" '$1==p && $3=="active"{f=1} END{exit !f}' "$REG"; then PORT=$p; break; fi; done
[ -n "$PORT" ] && echo "$PORT ${JID}_$I active $HOST $(date -Is)" >> "$REG"
flock -u 9
[ -n "$PORT" ] || { echo "no free port" >&2; exit 2; }
trap '( flock 9; sed -i "s/^$PORT ${JID}_$I active /$PORT ${JID}_$I released /" "$REG" ) 9>"$REG.lock"; rm -rf "$TMPD"' EXIT
echo "$(date -Is) job ${JID}_$I host $HOST cpus $CPUS port $PORT candidate $I $RUN eval $EVI engine $L3_ENGINE_SHA out $OUTF"
env -i HOME="$HOME" PATH="$PATH" LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" LANG=C.UTF-8 TZ=America/New_York TMPDIR="$TMPD" \
  R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null SLURM_CPUS_PER_TASK="$CPUS" SLURM_JOB_ID="$JID" R_PARALLEL_PORT="$PORT" \
  R_LIBS="$L3_LIB:$L3_DEPS_LIB" R_LIBS_USER="$L3_LIB" L3_LIB="$L3_LIB" L3_BATON_LIB="$L3_DEPS_LIB" L3_F0="$L3_F0" \
  L3_AL="$L3_AL" L3_RESULTS="$L3_RESULTS" L3_ENGINE_SHA="$L3_ENGINE_SHA" \
  Rscript "$L3_HARNESS/layer3_search.R" "$RUN" "$EVI" dev "$OUTF"
echo "$(date -Is) done ${JID}_$I"
