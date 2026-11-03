# Step 2a capture report (agent A, 2026-11-03)

Plan: `dev/plans/2026-11-03-engine-cpp-plan.md`, Step 2a (inventory items I1 to I6). Nothing committed, no package source edited.

## Result

All 10 designs reproduce their reference row on all four checks (power count, type I count, `EN_null` and `EN_alt` at 4 decimals), on local R 4.6.0 through `lib/ref` at R = 10,000, 19 workers, first cert seed. Item I6 (R 4.5.0 vs 4.6.0) did not arise. Every capture also passed the in-script assertions. Cutpoint count equals theta `n_intervals` for all 10, `random_dropout` is FALSE, `min_follow_up_at_final` is 6, the `looks_stop_at_final` option is TRUE (both at the hook and after the call), `num_workers` is 19, `num_simulations` and `n_rep` are 10,000. The wrapper made exactly one `run_scenarios` call, and the runner's engine-echo check returned `engine_ok TRUE` with the same echo string as the Longleaf `certify_16/evals.csv` (`p6[calls=1 F=6 opt=TRUE workers=38 ...] dropout[calls=1 switch=off passed=FALSE workers=38 seen=FALSE]`).

Object counts: 10 RDS files in `dev/perf/evidence/args/`, and 10 rows in `capture_check.csv` with all four `_ok` columns TRUE.

| design | task | protocol | cert seed | K theta | K cutpoints | power count ref / cap | type I count ref / cap | E0[N] cap (= ref) | E1[N] cap (= ref) | wall s |
|---|---|---|---|---|---|---|---|---|---|---|
| R01 | T04 | xphil | 900001 | 8 | 8 | 8245 / 8245 | 150 / 150 | 11.2470 | 28.2910 | 23.2 |
| R02 | T07 | xphil | 900001 | 8 | 8 | 8608 / 8608 | 145 / 145 | 13.2455 | 27.6548 | 15.8 |
| R03 | T03 | rws | 900001 | 9 | 9 | 9585 / 9585 | 111 / 111 | 21.4618 | 29.8526 | 16.0 |
| R04 | U01 | rws | 910001 | 4 | 4 | 8308 / 8308 | 57 / 57 | 13.6633 | 28.1529 | 9.9 |
| R23 | T11 | bws | 900001 | 15 | 15 | 8002 / 8002 | 826 / 826 | 58.6885 | 90.9617 | 165.6 |
| R24 | T17 | bws | 900001 | 9 | 9 | 8085 / 8085 | 892 / 892 | 63.9323 | 83.9193 | 56.3 |
| R25 | T22 | bws | 900001 | 12 | 12 | 8139 / 8139 | 939 / 939 | 70.5507 | 81.3908 | 45.9 |
| R26 | U13 | rws | 910001 | 13 | 13 | 8139 / 8139 | 853 / 853 | 60.9857 | 86.0771 | 100.5 |
| R27 | T23 | rws | 900001 | 14 | 14 | 8174 / 8174 | 958 / 958 | 74.2830 | 82.1041 | 37.9 |
| R28 | T24 | rws | 900001 | 15 | 15 | 8172 / 8172 | 976 / 976 | 66.7974 | 79.7318 | 52.2 |

Wall time is the wrapper call. Whole-process times (load included) were 11 to 166 s, 8.9 minutes for all ten, run one at a time on this 20-core machine, which was otherwise idle (load average 0.9 before the start).

## Commands

```bash
cd ~/research/evolveTrial-perf
mkdir -p dev/perf/lib/ref dev/perf/build    # .gitignore += dev/perf/lib/, dev/perf/build/, dev/perf/harness/cache/
bash ~/Downloads/atbo-longleaf/pilot/verify/engine/build_verify_source.sh ~/research/evolveTrial dev/perf/build
#   -> built p6 (7565a17), m4 (7565a17 + patch 83de2c8b2de1902630eddf00c334b508), m1 (422eb40 + patch d12a0e20901ab56f23b675ea91ed56e9) under dev/perf/build
R CMD INSTALL --no-multiarch -l dev/perf/lib/ref dev/perf/build/evolveTrial_m4    # log: dev/perf/build/install_ref.log
md5sum dev/perf/build/evolveTrial_m4/R/simulation_driver.R                      # 54763b77b65793392c016d07c14acd7d (same as R/simulation_driver.R on this branch)
# installed DESCRIPTION: P6RecertEvolveTrialSHA 7565a17c625573d2a697d0d4a85efb2fa6c58ad9, VerifyM4PatchMD5 83de2c8b2de1902630eddf00c334b508,
#                        VerifyM4DriverMD5 54763b77b65793392c016d07c14acd7d, Version 0.0.0.9001
bash dev/perf/build/run_captures.sh    # for d in R01 R02 R03 R04 R23 R24 R25 R26 R27 R28:
                                       #   Rscript dev/perf/evidence/scripts/capture_cert.R $d 10000 dev/perf/evidence/args/args_$d.rds
                                       # logs: dev/perf/build/capture_logs/<design>.log, summary.txt (all exit 0)
```

`capture_check.csv` was built with the script below (kept outside the repo, since it is not a listed output; reproduced here so the table can be rebuilt).

```r
setwd("~/research/evolveTrial-perf/dev/perf/evidence")
ref <- read.csv("acceptance_set_reference_rows.csv", stringsAsFactors = FALSE)
# per design: x <- readRDS(args/args_<d>.rds); r <- ref row with design == d and seed == x$cert_seed (exactly one)
# power_count_ok = round(power * 10000) == r$power_count; type1_count_ok likewise
# EN_null_ok = round(EN_null, 4) == round(r$EN_null, 4); EN_alt_ok likewise; plus K, P6 and dropout fields
```

`capture_cert.R` md5 at capture time: `1e508f09fc7d5ef6e1939cb4ac63566c`. Each RDS holds `design, task, scenario_id, trial_type, protocol, csv, theta, cert_seed, base_args, scens, dots` (`parallel FALSE, seed <cert seed>, return_variance TRUE`), `seed_before_run_scenarios` (`.Random.seed` at hook entry), `wrapper_res, counts, R, tier, runner, runner_md5, checks, evals_row` (the runner's evals.csv row with the engine echo), `wall_s, cpu, env, sessionInfo`.

Input checks the script makes before running. The runner files have the md5 values recorded in the Longleaf certify manifests (`retune_run.R` 6e2a3905, `retune2a_run.R` b388c79f, `p6r_protocol.R` 52d2b43c, `verify_m4_switch.R` 405f7d60). The runner's `STAGED.md5` and `f0_files.md5` lists pass. The F0 worktree is at 1863cecb. `cert_design.csv` md5 is 1e1997bb (t1) or ec1c987c (t2a). BATON is 0.4.0 from `~/R_libs_frozen/jasa_r1_v040`. RNGkind is Mersenne-Twister/Inversion/Rejection.

A smoke run at R = 200 (R23, written to the scratchpad, not kept) ran before the production captures. It checked the script end to end and nothing in the script was changed after it.

## Where the plan conflicts with the runner (the runner was followed)

1. **`BO_SEED`.** Task 2a.2 says the ENV block is set with `SEED` = the design's first cert seed. In the runner, `SEED` is the task's search seed (`retune_run.R:81`, `SEED <- as.integer(tr$seed)`), and the Longleaf manifests confirm it (`certify_1`: `BO_SEED=1060001` for R01/T04, `certify_51`: `BO_SEED=2240001` for R26/U13). The capture uses the task seed. The cert seed is passed only as the `seed` argument of `W`/`W0`, as at `retune_run.R:396`. `BO_SEED` is not read on the certify path, which acceptance_set.md Section 2.1 also notes.
2. **R04 and R26 use a different runner.** The plan says to find how `retune_run.R` resolves the t2a tasks. It does not: `retune_run.R:68` reads only `t1_tasks.tsv`. Tier 2a was certified by `pilot/retune2a/R/retune2a_run.R`, a marked copy, with task table `pilot/retune2a/staged/t2a_tasks.tsv`, staged md5 lists under `pilot/retune2a/staged/`, a literal tier 2a weight table (`:144-156`), a held/protocol check (`:86-88`), a minimax ev rule (`:162`), a staged n_init check (`:182`) and the `base_design` id (`:411-412`). The capture follows `retune2a_run.R` for R04 and R26. The task id comes from `results/retune_t2a/cert_design.csv` (`task` column: U01, U13), not from `t2a_cert_tasks.csv`, which holds only (cert_task, design, seed, status). R26's seed-910001 row comes from the logged rerun `c2a_r` (`retune2a_rerun_c2a.sh`), which only doubles the time limit and runs the same command.
3. **Protocols.** Inventory I2 names R27 and R28 as the `rws` designs. From the task tables, `rws` is R03, R04, R26, R27 and R28, `xphil` is R01 and R02, and `bws` is R23, R24 and R25. The runner's `switch(PROT, ...)` was copied verbatim, so every design was parsed the way the runner parsed it.
4. **The K-snap list in I1.** I1 lists R03, R23, R25, R26, R27 and R28 as having K outside {4, 6, 8, 12}. R25 has K = 12, which is inside that set, and R24 has K = 9, which is outside it. The correct list is R03, R23, R24, R26, R27, R28, as in acceptance_set.md. The capture sets `BATON_DISABLE_KSNAP=1` for every design, so the error in the list had no effect.
5. **Loading evolveTrial.** Section 0 says `library(evolveTrial, lib.loc = ...)`. The runner uses `.libPaths(c(LIB, ...))` and then `loadNamespace("evolveTrial", lib.loc = LIB)`, without attaching (`retune_run.R:105-107`). The capture copies the runner and adds the Section 0 path assertion. One side effect: `loadNamespace` at top level prints `<environment: namespace:evolveTrial>` in the log, as it would in the runner.
6. **`W0` versus `W`.** The plan says to call `W0(theta, fidelity = "high", seed = <cert seed>)`. The runner calls `W` (`:202-246`, which calls `W0` at `:215`) after `setwd(F0)` (`:393`). `W` adds the per-call engine-echo check. The capture copies `W` verbatim and calls it with `setwd(F0)`, so `W0` is still the function that runs. `W` consumes no random numbers.
7. **Frozen `cert_design.md5`.** The runner reads `<select dir>/cert_design.md5` (`:74`). `results/retune_t1` and `results/retune_t2a` hold no such file. The capture compares against the frozen md5 values as literals: 1e1997bb..., which equals `select_t1a/cert_design.md5` in the reproducibility build, and ec1c987c..., which equals `~/pilot_r1/retune_t2a/select_t2a/cert_design.md5` and `SELECTION.md:43`.
8. **Where BATON comes from.** On Longleaf, BATON loaded from the M4 library (`R_libs_verify_m4`, per the manifests). Locally it loads from `~/R_libs_frozen/jasa_r1_v040`, as the plan says. Both are version 0.4.0, and the runner's version and path checks pass.
9. **Not copied, because they do not apply to certification.** Smoke overrides (`:166-169`), the `pareto_front_select_best` coercion (`:175-186`, used by BO selection only), and the manifest write (`:249-260`). `R_PARALLEL_PORT` was set on Longleaf (for example 14125) and is unset locally. The runner marks it optional, and all 10 rows reproduced without it.
10. **Placement of the capture hook (I4).** This is not a conflict. Both installers capture whatever `run_scenarios` is in the namespace at install time (`verify_m4_switch.R:38`, `p6r_protocol.R:51`), so a hook installed first sits innermost. The script checks `identical(.vm4$orig_rs, hook)`. The captured `base_args` therefore already carry the P6 and dropout settings.

## Files

- `dev/perf/lib/ref/evolveTrial` (engine of record, gitignored), `dev/perf/build/` (source trees, install log, `run_captures.sh`, `capture_logs/`; gitignored)
- `dev/perf/evidence/scripts/capture_cert.R`
- `dev/perf/evidence/args/args_R01.rds` ... `args_R28.rds` (10), `dev/perf/evidence/args/capture_check.csv` (10 rows)
- `.gitignore` (three lines added)
