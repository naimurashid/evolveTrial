#!/bin/bash
# layer2_stage.sh -- stage layer 2 on Longleaf from the desktop (plan step 5, dev/plans/2026-11-03-step5-layer2.md,
# Section 3). Run from the worktree root. Submits nothing.
#   bash dev/perf/layer2/layer2_stage.sh
# 1. Builds the engine tree locally with build_dev_tree.sh (the commit in layer2_env.sh) and copies it to $L2/src.
# 2. Copies the atbo-longleaf files layer2_cert.R reads into $L2/al under their repo paths, and the harness files.
# 3. On the login node: compares every copied file's md5 with the desktop's, installs the tree into $L2/lib with the
#    certify modules, and checks the installed library (both DESCRIPTION fields, the C++ entry point, BATON 0.4.0 from
#    the deps library, evolveTrial loaded from $L2/lib).
set -euo pipefail
WT=$(pwd); [ -f "$WT/dev/perf/layer2/layer2_env.sh" ] || { echo "run from the worktree root" >&2; exit 2; }
. "$WT/dev/perf/layer2/layer2_env.sh"
H=longleaf.unc.edu; AL=${AL_LOCAL:-$HOME/Downloads/atbo-longleaf}
[ -z "$(git status --porcelain -- dev/perf/layer2 dev/perf/evidence/acceptance_set_reference_rows.csv)" ] ||
  { echo "layer 2 files are uncommitted or differ from HEAD: commit first" >&2; exit 2; }
[ "$(md5sum < dev/perf/evidence/acceptance_set_reference_rows.csv | cut -c1-32)" = "$L2_REF_ROWS_MD5" ] || { echo "reference rows md5 != env" >&2; exit 2; }
[ "$(md5sum < dev/perf/layer2/layer2_cert.R | cut -c1-32)" = "$L2_CERT_MD5" ] || { echo "layer2_cert.R md5 != env" >&2; exit 2; }
ssh -o BatchMode=yes -o ConnectTimeout=20 "$H" true || { echo "cannot reach $H by ssh (VPN?)" >&2; exit 2; }
ssh -o BatchMode=yes -o ConnectTimeout=20 "$H" "[ ! -e '$L2_LIB/evolveTrial' ] && [ ! -e '$L2_TREE' ]" ||
  { echo "refusing: $L2_LIB/evolveTrial or $L2_TREE exists on Longleaf" >&2; exit 2; }
ALF=(pilot/retune/R/retune_run.R pilot/retune2a/R/retune2a_run.R pilot/p6recert/R/p6r_protocol.R pilot/verify/R/verify_m4_switch.R
     results/retune_t1/cert_design.csv results/retune_t2a/cert_design.csv)
mapfile -t STG < <(cd "$AL" && find pilot/retune/staged pilot/retune2a/staged -type f | sort)
HF=(dev/perf/layer2/layer2_cert.R dev/perf/layer2/layer2_compare.R dev/perf/layer2/layer2_env.sh dev/perf/layer2/layer2.sbatch
    dev/perf/layer2/layer2_submit.sh dev/perf/layer2/layer2_gate.sh dev/perf/evidence/acceptance_set_reference_rows.csv)
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
bash dev/perf/layer2/build_dev_tree.sh "$WT" "$L2_ENGINE_SHA" "$T/$(basename "$L2_TREE")"
(cd "$AL" && md5sum "${ALF[@]}" "${STG[@]}") > "$T/al.md5"
(cd "$WT" && md5sum "${HF[@]}" | sed 's#  dev/perf/[a-z0-9]*/#  #') > "$T/harness.md5"
(cd "$T" && find "$(basename "$L2_TREE")" -type f | sort | xargs md5sum) > "$T/tree.md5"
ssh -o BatchMode=yes "$H" "mkdir -p '$L2/src' '$L2/al' '$L2_HARNESS' '$L2_LIB' '$L2/logs' '$L2/out' '$L2/ports'"
rsync -a "$T/$(basename "$L2_TREE")" "$H:$L2/src/"
(cd "$AL" && rsync -aR "${ALF[@]}" "${STG[@]}" "$H:$L2_AL/")
rsync -a "${HF[@]}" "$H:$L2_HARNESS/"
rsync -a "$T/al.md5" "$T/harness.md5" "$T/tree.md5" "$H:$L2/"
ssh -o BatchMode=yes "$H" bash -s <<EOF
set -euo pipefail
. '$L2_HARNESS/layer2_env.sh'
(cd "\$L2_AL" && md5sum -c --quiet "\$L2/al.md5") && (cd "\$L2_HARNESS" && md5sum -c --quiet "\$L2/harness.md5") &&
  (cd "\$L2/src" && md5sum -c --quiet "\$L2/tree.md5") || { echo "md5 mismatch after copy" >&2; exit 3; }
echo "md5 ok: \$(wc -l < "\$L2/al.md5") al, \$(wc -l < "\$L2/harness.md5") harness, \$(wc -l < "\$L2/tree.md5") tree files"
l2_modules
R_LIBS="\$L2_LIB:\$L2_DEPS_LIB" R_LIBS_USER="\$L2_LIB" R CMD INSTALL --no-multiarch -l "\$L2_LIB" "\$L2_TREE" > "\$L2/logs/install.log" 2>&1 || { tail -20 "\$L2/logs/install.log"; exit 3; }
export L2_LIB L2_DEPS_LIB L2_ENGINE_SHA
R_LIBS="\$L2_LIB:\$L2_DEPS_LIB" R_LIBS_USER="\$L2_LIB" Rscript -e '
  lib <- normalizePath(Sys.getenv("L2_LIB")); loadNamespace("evolveTrial", lib.loc = lib)
  d <- utils::packageDescription("evolveTrial", lib.loc = lib)
  stopifnot(identical(d\$EngineCppSHA, Sys.getenv("L2_ENGINE_SHA")), identical(d\$P6RecertEvolveTrialSHA, "7565a17c625573d2a697d0d4a85efb2fa6c58ad9"),
            is.null(d\$VerifyM4PatchMD5), exists(".simulate_chunk_cpp", envir = asNamespace("evolveTrial"), inherits = FALSE))
  suppressPackageStartupMessages(library(BATON, lib.loc = Sys.getenv("L2_DEPS_LIB")))
  stopifnot(identical(as.character(utils::packageVersion("BATON")), "0.4.0"),
            dirname(normalizePath(getNamespaceInfo("evolveTrial", "path"))) == lib)
  cat("library ok:", R.version.string, "evolveTrial", d\$EngineCppSHA, "BATON", as.character(utils::packageVersion("BATON")), "\n")'
EOF
echo "staged under $H:$L2 (nothing submitted)"
