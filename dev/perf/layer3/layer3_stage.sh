#!/bin/bash
# layer3_stage.sh -- stage layer 3 on Longleaf from the desktop (plan step 6, dev/plans/2026-11-03-step6-layer3.md). Run from
# the worktree root. Submits nothing and installs nothing.
#   bash dev/perf/layer3/layer3_stage.sh
# 1. Copies the harness files to $L3_HARNESS and checks their md5 after the copy.
# 2. On the login node: asserts the layer 2 library reports EngineCppSHA = L3_ENGINE_SHA, copies each candidate's ledger and rds
#    from $L3_E2C to $L3_RESULTS by path, and compares their sha256 with the desktop copies in ~/pilot_r1/e2c/results_live.
set -euo pipefail
WT=$(pwd); [ -f "$WT/dev/perf/layer3/layer3_env.sh" ] || { echo "run from the worktree root" >&2; exit 2; }
. "$WT/dev/perf/layer3/layer3_env.sh"
H=longleaf.unc.edu; LOC=${E2C_LOCAL:-$HOME/pilot_r1/e2c/results_live}
[ -z "$(git status --porcelain -- dev/perf/layer3 ':!dev/perf/layer3/out' ':!dev/perf/layer3/BUILD_REPORT.md')" ] ||
  { echo "layer 3 files are uncommitted or differ from HEAD: commit first" >&2; exit 2; }
ssh -o BatchMode=yes -o ConnectTimeout=20 "$H" true || { echo "cannot reach $H by ssh (VPN?)" >&2; exit 2; }
ssh -o BatchMode=yes "$H" "[ ! -e '$L3_HARNESS' ] && [ ! -e '$L3_RESULTS' ]" || { echo "refusing: $L3_HARNESS or $L3_RESULTS exists on Longleaf" >&2; exit 2; }
HF=(dev/perf/layer3/layer3_search.R dev/perf/layer3/layer3_compare.R dev/perf/layer3/layer3_env.sh dev/perf/layer3/layer3.sbatch
    dev/perf/layer3/layer3_submit.sh dev/perf/layer3/layer3_gate.sh)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
(cd "$WT" && md5sum "${HF[@]}" | sed 's#  dev/perf/layer3/#  #') > "$T/harness.md5"
: > "$T/results.sha256"
for c in "${L3_CANDS[@]}"; do read -r r e <<< "$c"
  for f in "$r/sim_ledger.csv" "$(l3_rds_path "$r")"; do (cd "$LOC" && sha256sum "$f") >> "$T/results.sha256"; done; done
ssh -o BatchMode=yes "$H" "mkdir -p '$L3_HARNESS' '$L3/logs' '$L3/out' '$L3/ports' '$L3_RESULTS'"
rsync -a "${HF[@]}" "$H:$L3_HARNESS/"
rsync -a "$T/harness.md5" "$T/results.sha256" "$H:$L3/"
ssh -o BatchMode=yes "$H" bash -s <<EOS
set -euo pipefail
. '$L3_HARNESS/layer3_env.sh'
(cd "\$L3_HARNESS" && md5sum -c --quiet "\$L3/harness.md5") || { echo "harness md5 mismatch after copy" >&2; exit 3; }
grep -qx "EngineCppSHA: \$L3_ENGINE_SHA" "\$L3_LIB/evolveTrial/DESCRIPTION" || { echo "library \$L3_LIB is not engine \$L3_ENGINE_SHA" >&2; exit 3; }
while read -r sha f; do mkdir -p "\$L3_RESULTS/\$(dirname "\$f")"; cp -p "\$L3_E2C/\$f" "\$L3_RESULTS/\$f"; done < "\$L3/results.sha256"
(cd "\$L3_RESULTS" && sha256sum -c --quiet "\$L3/results.sha256") || { echo "results sha256 mismatch after copy" >&2; exit 3; }
for p in "\$L3_AL/pilot/p6recert/R/p6r_protocol.R" "\$L3_AL/pilot/verify/R/verify_m4_switch.R" "\$L3_AL/pilot/retune/staged/f0_files.md5" "\$L3_F0/scenarios/cohortB_multiarm_production.csv"; do
  [ -f "\$p" ] || { echo "missing \$p" >&2; exit 3; }; done
echo "staged ok: \$(wc -l < "\$L3/harness.md5") harness files, \$(wc -l < "\$L3/results.sha256") result files, library engine \$L3_ENGINE_SHA"
EOS
echo "staged under $H:$L3 (nothing submitted)"
