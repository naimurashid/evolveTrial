# layer3_env.sh -- Longleaf paths of layer 3 (plan step 6, dev/plans/2026-11-03-step6-layer3.md). Sourced by layer3.sbatch,
# layer3_submit.sh and layer3_gate.sh. Read-only inputs: the layer 2 staged library and atbo files, the F0 worktree, the E2c
# results. Writes only under $L3.
L3=/proj/rashidlab/nur2/e2_ll/engine_cpp_l3
L2=/proj/rashidlab/nur2/e2_ll/engine_cpp_l2              # layer 2 staging, read-only here
L3_ENGINE_SHA=1c918711f48c411fecd1c1dad7e1bdd799171fd9   # as layer 2; the staged library's EngineCppSHA is asserted to equal it
L3_LIB=$L2/lib                                           # the library layer 2 passed with (evolveTrial only)
L3_DEPS_LIB=/proj/rashidlab/nur2/e2_ll/R_libs_verify_m4  # BATON 0.4.0, tidyverse
L3_F0=/proj/rashidlab/nur2/pilot_ll/f0                   # F0 worktree (git 1863cecb)
L3_AL=$L2/al                                             # md5-checked atbo-longleaf files staged by layer 2 (p6r, m4s, f0_files.md5)
L3_E2C=/proj/rashidlab/nur2/e2_ll/e2c/results            # the stored E2c rows (source of the staged copies)
L3_RESULTS=$L3/results                                   # md5-checked copies of the three candidates' ledger and rds
L3_HARNESS=$L3/harness                                   # layer3_search.R, layer3_compare.R and these scripts
L3_ACCOUNT=rc_nur2_pi
L3_PARTITION=general
L3_CPUS=16                                               # E2c search worker count
L3_PORTS=16000                                           # FORK cluster port range 16000-16499 (live ports are skipped by ss)
# candidates: task index, run id, eval index
L3_CANDS=("PB-baton-1001 1" "PB-baton-1003 95" "PB-ga-1001 274")
l3_rds_path() { case "$1" in PB-baton-*) echo "$1/f0_out/doublet_triplet_admissible_fut_warmstart_result.rds";; *) echo "$1/run_state.rds";; esac; }
l3_modules() {
  if ! type module >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    . /etc/profile.d/modules.sh >/dev/null 2>&1 || . /usr/share/lmod/lmod/init/bash
  fi
  module load gcc/12.2.0 r/4.5.0 >/dev/null 2>&1
  command -v Rscript >/dev/null || { echo "Rscript not on PATH after module load" >&2; return 1; }
}
