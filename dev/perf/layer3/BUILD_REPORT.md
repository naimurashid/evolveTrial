# Layer 3 build report (2026-11-03)

Built in `dev/perf/layer3/`, nothing committed, nothing run on Longleaf, nothing installed. Record Section 5 filled.

## Result
- Local `ref`: 3 of 3 pass (counts equal ledger, objective/E0[N]/E1[N] identical() to rds doubles).
- Local `dev` (lib/dev_l2, 1c91871): 3 of 3 pass.
- Negative control (`L3_FORCE_R=1`, C2 dev): C++ probe lines 0, trials 0 against 6,000 expected, engine check failed, exit 3, no rds.
- Walls (s): ref 29.1, 15.8, 16.9. dev 1.5, 1.1, 1.2. Desktop was loaded by other jobs.

## Files
layer3_search.R, layer3_compare.R, layer3_env.sh, layer3.sbatch, layer3_stage.sh, layer3_submit.sh, layer3_gate.sh (bash -n clean, untested on Longleaf), run_local.sh (local driver, not for Longleaf).

## Conflicts and deviations from the brief, record or recipe
1. **lib/dev has no EngineCppSHA.** Record decision 5 and the brief named `lib/dev`. Coordinator ruled `lib/dev_l2`, as layer2_cert.R does. Record corrected by the coordinator.
2. **Library path for tidyverse.** Recon does not say where tidyverse comes from locally. `jasa_r1_v040` holds only BATON and evolveTrial, so `run_local.sh` sets `R_LIBS=$HOME/R/library` (read-only use). On Longleaf `L3_DEPS_LIB` (R_libs_verify_m4) supplies it, as in layer 2.
3. **rds read order.** Reading the BATON rds loads DiceKriging, which trips the "no packages loaded before the env block" assertion. The rds is therefore read after the libraries load (the ledger csv is still read first). The assertion is kept unchanged.
4. **Packages loaded per row type.** BATON rows load tidyverse and BATON 0.4.0, GA rows load neither explicitly (scenario_parser.R loads tidyverse itself), following recon L3-c.
5. **Dropped from layer 2:** runner and cert_design md5 checks, STAGED.md5, task table and weight tables (they belong to the certify runners, not E2c). Kept: P6R and M4S md5, f0_files.md5, F0 git SHA.
6. **Objective source.** Objective is taken from the wrapper result field `weighted_EN_N`. It was identical() to the rds objective on all rows.
7. **Compare takes an optional 4th argument `ref|dev`** (default dev), so the same script gates the local ref run. The Longleaf gate calls it with `dev`. It also checks sha256 of the ledger and rds against the values recorded in each output.
8. **Stage script** requires layer3 files committed (as layer 2 does), excluding `out/` and this report. It copies by path on Longleaf from `$L3_E2C`, not rsync from the desktop, then checks sha256 against the desktop copies in `~/pilot_r1/e2c/results_live`.
9. **Port range** 16000-16499 is shared with layer 2, with a separate registry under `$L3/ports`. Live ports are skipped by `ss`. Not an issue unless layer 2 jobs run at the same time on one host.
10. `out/local/` is under `dev/perf/layer3/` and not gitignored.
