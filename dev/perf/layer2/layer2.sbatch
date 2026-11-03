#!/bin/bash
#SBATCH --job-name=l2cpp
#SBATCH --account=rc_nur2_pi
#SBATCH --partition=general
#SBATCH --cpus-per-task=19
#SBATCH --mem=32G
#SBATCH --time=00:15:00
#SBATCH --no-requeue
#SBATCH --hint=nomultithread
# layer2.sbatch -- one array task per reference row (plan step 5, dev/plans/2026-11-03-step5-layer2.md, Section 3).
#   sbatch --array=1-75 --output=$L2/logs/<tag>/%A_%a.log layer2.sbatch <tag>
# Task i runs row i of acceptance_set_reference_rows.csv through layer2_cert.R with the dev engine at R = 10,000 and 19
# workers. The run block copies retune.sbatch's run() (clean environment, fixed FORK cluster port from a per-host
# registry, R_ENVIRON_USER and R_PROFILE_USER disabled) with this program's paths. Writes only under $L2.
set -euo pipefail
. "${L2_ENV:-/proj/rashidlab/nur2/e2_ll/engine_cpp_l2/harness/layer2_env.sh}"
TAG=${1:?tag}; I=${SLURM_ARRAY_TASK_ID:?array task}; JID=${SLURM_JOB_ID:?job}; CPUS=${SLURM_CPUS_PER_TASK:?cpus}
[ "$CPUS" = "$L2_CPUS" ] || { echo "cpus $CPUS != $L2_CPUS" >&2; exit 2; }
m() { md5sum < "$1" | cut -c1-32; }
[ "$(m "$L2_REF_ROWS")" = "$L2_REF_ROWS_MD5" ] || { echo "reference rows md5 mismatch" >&2; exit 2; }
[ "$(m "$L2_HARNESS/layer2_cert.R")" = "$L2_CERT_MD5" ] || { echo "layer2_cert.R md5 mismatch" >&2; exit 2; }
read -r DES SEED < <(awk -F, -v i="$I" 'NR == 1 { if ($1 != "design" || $7 != "seed") exit 3 } NR == i + 1 { print $1, $7 }' "$L2_REF_ROWS")
[ -n "${DES:-}" ] && [ -n "${SEED:-}" ] || { echo "no reference row $I" >&2; exit 2; }
OUTD=$L2/out/$TAG; mkdir -p "$OUTD" "$L2/ports"; OUTF=$OUTD/${DES}_${SEED}.rds
[ ! -e "$OUTF" ] || { echo "refusing: $OUTF exists" >&2; exit 2; }
l2_modules
TMPD=$(mktemp -d "/tmp/l2cpp_${JID}.XXXXXX")
HOST=$(hostname -s); REG=$L2/ports/$HOST.reg; touch "$REG"; PORT=""
exec 9>"$REG.lock"; flock 9
for ((k = 0; k < 500; k++)); do p=$((L2_PORTS + (JID + I + k) % 500))
  if [ -z "$(ss -Htan "( sport = :$p )" 2>/dev/null)" ] && ! awk -v p="$p" '$1==p && $3=="active"{f=1} END{exit !f}' "$REG"; then PORT=$p; break; fi; done
[ -n "$PORT" ] && echo "$PORT ${JID}_$I active $HOST $(date -Is)" >> "$REG"
flock -u 9
[ -n "$PORT" ] || { echo "no free port" >&2; exit 2; }
trap '( flock 9; sed -i "s/^$PORT ${JID}_$I active /$PORT ${JID}_$I released /" "$REG" ) 9>"$REG.lock"; rm -rf "$TMPD"' EXIT
echo "$(date -Is) job ${JID}_$I host $HOST cpus $CPUS port $PORT row $I $DES $SEED engine $L2_ENGINE_SHA out $OUTF"
env -i HOME="$HOME" PATH="$PATH" LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}" LANG=C.UTF-8 TZ=America/New_York TMPDIR="$TMPD" \
  R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null SLURM_CPUS_PER_TASK="$CPUS" SLURM_JOB_ID="$JID" R_PARALLEL_PORT="$PORT" \
  R_LIBS="$L2_LIB:$L2_DEPS_LIB" R_LIBS_USER="$L2_LIB" L2_LIB="$L2_LIB" L2_BATON_LIB="$L2_DEPS_LIB" L2_F0="$L2_F0" \
  L2_AL="$L2_AL" L2_ENGINE_SHA="$L2_ENGINE_SHA" \
  Rscript "$L2_HARNESS/layer2_cert.R" "$DES" "$SEED" 10000 dev "$OUTF"
echo "$(date -Is) done ${JID}_$I"
