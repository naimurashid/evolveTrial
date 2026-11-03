# layer2_env.sh -- Longleaf paths of layer 2 (plan step 5, dev/plans/2026-11-03-step5-layer2.md). Sourced by
# layer2.sbatch and layer2_submit.sh. Read-only inputs: the M4 library (BATON 0.4.0 and every dependency, the reference
# evolveTrial is never loaded from it here) and the F0 worktree. Writes only under $L2.
L2=/proj/rashidlab/nur2/e2_ll/engine_cpp_l2
L2_ENGINE_SHA=1c918711f48c411fecd1c1dad7e1bdd799171fd9   # step 4 commit; R/, src/, DESCRIPTION identical at 830fc73
L2_TREE=$L2/src/l2_tree_${L2_ENGINE_SHA:0:12}          # build_dev_tree.sh output, copied from the desktop
L2_LIB=$L2/lib                                           # evolveTrial (the tree above) and nothing else
L2_DEPS_LIB=/proj/rashidlab/nur2/e2_ll/R_libs_verify_m4  # RETUNE_LIB of the certify runs: BATON 0.4.0, tidyverse
L2_F0=/proj/rashidlab/nur2/pilot_ll/f0                   # RETUNE_F0_DIR of the certify runs (git 1863cecb)
L2_AL=$L2/al                                             # md5-checked copy of the atbo-longleaf files layer2_cert.R reads
L2_HARNESS=$L2/harness                                   # layer2_cert.R, layer2_compare.R, the reference rows
L2_REF_ROWS=$L2_HARNESS/acceptance_set_reference_rows.csv
L2_REF_ROWS_MD5=be671e759e60f2b0e0de400c428cd26f
L2_CERT_MD5=548fd36ad2bf5e4c87fe3b5a28c5f3b1
L2_ACCOUNT=rc_nur2_pi
L2_PARTITION=general
L2_CPUS=19                                               # certify worker count (plan I8)
L2_PORTS=16000                                           # FORK cluster port range 16000-16499; 11000-15499 are taken by other harnesses
l2_modules() {
  if ! type module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    . /etc/profile.d/modules.sh >/dev/null 2>&1 || . /usr/share/lmod/lmod/init/bash
  fi
  module load gcc/12.2.0 r/4.5.0 >/dev/null 2>&1
  command -v Rscript >/dev/null || { echo "Rscript not on PATH after module load" >&2; return 1; }
}
