# evolveTrial engine C++ port, phase 1 (bit-identical): implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (or subagent-driven-development) to implement this plan task by task.

**Goal.** Replace the body of `simulate_chunk` with one C++ call per chunk that reproduces every number of the engine of record (58a4850) exactly, 10x or better on between-arm designs.

**Architecture.** The R driver (seeding, `splitIndices`, FORK workers, `aggregate_results`, post-processing) stays as it is. A C++ function `simulate_chunk_cpp` runs one chunk on R's own random number stream through libR's `R::` samplers, in the draw order recorded in `dev/perf/evidence/rng_driver_inventory.md`. A seam in `run_simulation_pure` exposes the R body and the C++ body to a test harness that compares them chunk by chunk and bisects to the first diverging trial.

**Tech stack.** R 4.6.0 locally and R 4.5.0 on Longleaf (gcc/12.2.0, r/4.5.0), Rcpp attributes, testthat 3, Slurm on Longleaf for layer 2.

**Design of record.** `dev/plans/2026-11-03-engine-cpp-design.md` (approved 2026-11-03). Evidence in `dev/perf/evidence/`. Chunk-body inventory in `dev/perf/evidence/chunk_body_inventory.md`.

---

## 0. Conventions for every step

- **Worktree.** All work is in `~/research/evolveTrial-perf` on branch `perf/engine-cpp`. Never touch `~/research/evolveTrial` (main), the other worktrees, `~/R/library`, `~/R_libs_frozen`, or any `atbo-*` repo. Those repos are read-only inputs.
- **Libraries.** Two private libraries, both gitignored.
  - `dev/perf/lib/ref`: the engine of record, built by `~/Downloads/atbo-longleaf/pilot/verify/engine/build_verify_source.sh` (tree `evolveTrial_m4`, which is 7565a17 plus `m4_nodropout.patch` plus the DESCRIPTION fields `P6RecertEvolveTrialSHA`, `VerifyM4PatchMD5`, `VerifyM4DriverMD5`). Its `R/simulation_driver.R` has md5 `54763b77b65793392c016d07c14acd7d`, the same as this branch at c2a88d5. Built once, never rebuilt.
  - `dev/perf/lib/dev`: this branch, reinstalled after every source change with `R CMD INSTALL --no-multiarch --preclean -l dev/perf/lib/dev . && git checkout -- src/symbols.rds`. Keep `--preclean`: R's package make has no header dependencies, so without it an edit to `src/r_semantics.h` leaves stale objects. `--preclean` deletes the tracked `src/symbols.rds`, hence the checkout (inventory B7).
  - Every script loads evolveTrial with `library(evolveTrial, lib.loc = normalizePath(<lib>))` (absolute, since testthat changes directory) and then asserts `dirname(getNamespaceInfo("evolveTrial", "path")) == normalizePath(<lib>)`.
- **Gate commands** (built in step 2b, all run from the worktree root).
  - Layer 0: `Rscript dev/perf/harness/layer0_seam_check.R dev/perf/lib/<ref|dev> dev/perf/evidence/args/args_<D>.rds <out.rds>` for both libraries, then `Rscript dev/perf/harness/layer0_compare.R <ref.rds> <dev.rds>`, for D in R01 and R23.
  - Layer 1: `Rscript dev/perf/harness/run_layer1.R [--designs R01,...] [--cand simulate_chunk_cpp]`. Exits 1 on any mismatch and writes `dev/perf/harness/out/layer1_<timestamp>.csv`.
  - Harness self-tests: `Rscript -e 'library(evolveTrial, lib.loc = normalizePath("dev/perf/lib/dev")); testthat::test_file("dev/perf/harness/tests/test-harness.R")'`.
  - Package tests: `Rscript dev/perf/harness/run_pkg_tests.R dev/perf/lib/dev <out.rds>`. Baseline at c2a88d5: 314 passed, 0 failed, 33 skipped.
- **BATON.** Only the capture step needs it, version 0.4.0 from `~/R_libs_frozen/jasa_r1_v040`.
- **RNG.** Every harness entry point asserts `identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))`.
- **Delegation.** Agents never commit. The coordinator verifies each agent's central claim independently, by rerunning the named command and counting objects, and makes every commit. Agent briefs say "report where this plan conflicts with the code you find, do not improvise", and every agent writes its report to a file and returns at most 150 words.
- **Commits.** One commit per step below, with the gate result in the message. Nothing in the C++ port is committed with layer 1 failing.

## 1. Breakage inventory

What can break, where, what it would do, and which check catches it. "Step" is the step that owns the guard.

### 1a. Inputs to the gate

| ID | What breaks | Where | Consequence | Guard | Step |
|---|---|---|---|---|---|
| I1 | `args_R23.rds` was captured with the K-snap on. It has 12 intervals, the certified R23 has 15. The capture script set none of the certify environment block (`BATON_DISABLE_KSNAP=1`, `BO_*`, `SLURM_CPUS_PER_TASK`). | `dev/perf/evidence/scripts/capture.R` vs `atbo-longleaf/pilot/retune/R/retune_run.R:88-93` | Layer 1 and layer 2 would test the wrong design for every design with K outside {4, 6, 8, 12} (R03 9, R23 15, R24 9, R26 13, R27 14, R28 15). | Recapture all 10 with the full environment block. Assert cutpoint count equals theta `n_intervals`. **Closed 2026-11-03**, 10 of 10. | 2a |
| I2 | `capture.R` parses every design with `parse_single_scenario`. The certify runners use `parse_warmstart_scenarios` for protocol `rws` (R03, R04, R26, R27, R28). | `retune_run.R:124-127` | Captured base_args could differ from what certification used. | Copy the runner's `spec <- switch(PROT, ...)` verbatim. **Closed.** | 2a |
| I3 | `capture.R` reads only `results/retune_t1/cert_design.csv`. R04 and R26 were certified by a different runner, `pilot/retune2a/R/retune2a_run.R` with `t2a_tasks.tsv`. | `capture.R:12-15` | R04 and R26 cannot be captured. | Follow `retune2a_run.R` for those two. **Closed.** | 2a |
| I4 | The capture hook must sit innermost, under both `verify_m4_install("off")` and `p6r_install("p6")`, to see the arguments that reach the real `run_scenarios` (`random_dropout = FALSE`, `min_follow_up_at_final = 6`). | `retune_run.R:189-192` | Captured args miss the P6 or dropout settings. | Install the hook before the two installers. Assert both fields in the capture. | 2a |
| I5 | Nothing yet shows the captured args reproduce a reference row. | none | The harness could be exact against the wrong inputs. | Run the captured call at R = 10,000, 19 workers, cert seed, through `lib/ref`, and compare with `acceptance_set_reference_rows.csv` (`power_count`, `type1_count` as integers, `EN_null`, `EN_alt` at 4 decimals). | 2a |
| I6 | Local R is 4.6.0 and the reference rows come from R 4.5.0. A mismatch in I5 is ambiguous between a capture defect and a libR difference. | design Section 4, Risk 1 | Wrong attribution. | If I5 fails, stop and report. The coordinator reruns the same capture on Longleaf. **Result 2026-11-03:** I5 passed 10 of 10 locally, so R 4.6.0 reproduces the R 4.5.0 rows at the first cert seed of every design. Layer 2 still runs on Longleaf. | 2a |
| I7 | `looks_stop_at_final` is a lazy default, `getOption("evolveTrial.looks_stop_at_final", FALSE)`. If a harness restores options before the promise is forced, the chunk runs the non-P6 branch. | `simulation_driver.R:236` | Silent test of the wrong branch, in both bodies alike, so identity still passes. | The harness sets the option for the whole session and asserts `isTRUE(env$looks_stop_at_final)` after building the frame. | 2b |
| I8 | The certify runner ran with 19 workers, searches with 16. Results depend on the worker count. | `rng_driver_inventory.md` Section 2 | A gate at the wrong worker count proves nothing about the stored rows. | Layer 1 runs both counts. Layer 2 uses 19 and asserts it. Layer 3 uses 16. | 2b, 5, 6 |
| I9 | Below `evolveTrial.sequential_threshold` (100) the driver runs one sequential chunk with `seed = NULL`, so layer 1 at R = 50 does not follow the production chunk plan. | `simulation_driver.R:965-974` | None for chunk identity, which is what layer 1 tests. The plan emulation is checked separately at R = 200. | Harness self-test H3. | 2b |

### 1b. The seam and the build

| ID | What breaks | Where | Consequence | Guard | Step |
|---|---|---|---|---|---|
| B1 | The seam (`simulate_chunk_r` rename plus the option-gated early return) changes the R path. | `simulation_driver.R:468`, `:944` | Every later comparison is against a modified reference. | Layer 0: `run_scenarios` on R01 and R23 at R = 200, 4 workers, seed 900001, `lib/ref` vs `lib/dev`, `identical()` on the returned table. Rerun after every commit that touches R code. | 2b |
| B2 | `src/Makevars` is new (none exists). Writing it can drop flags or libraries the build gets by default. RcppArmadillo normally wants `PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)`. The package links today without them. | `src/` | Link failure, or a changed numerical path. | Makevars holds `PKG_CXXFLAGS = -ffp-contract=off` and nothing else. Layer 0 covers the existing C++ gamma loop. | 2b |
| B3 | `Rcpp::compileAttributes()` regenerates `R/RcppExports.R` and `src/RcppExports.cpp`. A different local Rcpp can rewrite unrelated wrappers, including their `RNGScope`. | `src/RcppExports.cpp:61` (the gamma loop relies on the wrapper's RNGScope) | The single-arm gamma draws stop reaching R's stream state. | Review `git diff` of both files after every regeneration. Only additions are allowed. Floating-point helpers are exported with `rng = false`. `simulate_chunk_cpp` keeps the default `rng = true`. | 2b, 3 |
| B4 | `R CMD INSTALL` from the worktree leaves `.o` and `.so` files in `src/`. | `src/.gitignore` | Stray binaries committed. | Check `src/.gitignore` covers them, `git status` before every commit. | 2b |
| B5 | `-ffp-contract=off` is a no-op on x86-64 without `-march` (no FMA). It would change results on aarch64. | design Section 2 | None on Linux x86-64. Mac builds would differ from Longleaf before and after. | Recorded only. Not a supported target. | none |
| B6 | Installing into `~/R/library` would replace the evolveTrial that other projects use. | `.libPaths()[1]` | Other projects silently run the branch engine. | Every install names `-l dev/perf/lib/...`. | all |
| B7 | `R CMD INSTALL --preclean` deletes the tracked `src/symbols.rds`. Without `--preclean`, header edits leave stale objects. | `src/` | A deleted tracked file in a commit, or a library built from stale objects. | Install command in Section 0 restores the file. `git status` before every commit. | all |
| B8 | The harness runner fell back to the superseded K-snapped `evidence/args_R23.rds` when a design was missing from `args/`. | `run_layer1.R:26-28` (step 2b) | Identity tested on the wrong design. | Fallback removed. Old captures moved to `dev/perf/evidence/superseded/` with a README. Self-tests repointed to `args/`. **Closed.** | 2b |

### 1c. The port itself

These are the defect classes the port can introduce. Layer 1 catches all of them, and the bisection tool localizes them to a trial. Items marked (inv) come from `dev/perf/evidence/chunk_body_inventory.md`.

| ID | What breaks | Where | Guard | Step |
|---|---|---|---|---|
| P1 | Draw order within a trial. The arrival `rexp` comes before every break test, the look sits between the arrival draw and the arm draw, and the arrival that ends accrual has already drawn. | `simulation_driver.R:542-679` | Layer 1, `.Random.seed` after each prefix. | 3 |
| P2 | Arm assignment uses one `unif_rand` through `ProbSampleNoReplace` and `revsort`, even with one eligible arm. With equal probabilities `revsort` puts the second arm first. | `simulation_driver.R:652` | Unit test against `sample()` on 1, 2 and 3 arms with equal and unequal probabilities. | 3 |
| P3 | The dropout uniform is drawn and discarded when dropout is off. `runif(1, 0, Inf)` draws nothing. | `simulation_driver.R:657-658` | Copy `R::runif` exactly, no shortcut. Layer 1. | 3 |
| P4 | `R::rexp` takes a scale, not a rate. | `rng_driver_inventory.md` Section 4 | Layer 1 on the first trial. | 3 |
| P5 | Between-arm draws: `num_posterior_draws` (800) at the vs-reference look, not the interim count, `R::rnorm` then `R::rgamma` per sample in R's loop order. The dead `sample_vs_ref_medians_ph_cpp` (`arma::randn`) must not be reused. | `posterior_helpers.R:105-135`, `interim_logic.R:84-169` | Layer 1 on R23 to R28. Grep for `arma::randn` in the new code. | 4 |
| P6 | Discarded control-cache draws (V0) when more than one experimental arm is evaluated and the PH model is off. Not reached by the frozen designs. | `interim_logic.R:116-119` | Out of the frozen configuration. The port refuses it (see P12). | 4 |
| P7 | Long-double reductions: R's `sum` and `mean` on doubles accumulate in long double, and `mean` adds a second correction pass. `stats::median` is the `mean` of the two middle values. Type 7 quantile. `findInterval` semantics plus `pmin(idx, K)`. | `state_management.R:115-149`, `simulation_driver.R:1331-1398` | Unit tests on each helper against the R function (step 2b), then layer 1. | 2b, 3 |
| P8 | Rebalanced cutpoints once 50 or more events are observed, reset to the baseline at the start of every trial. | `simulation_driver.R:501-520`, `:572-603` | Layer 1, long chunk. | 3 |
| P9 | Per-trial accumulators use plain double `+` in trial order, then across chunks in chunk order. | `simulation_driver.R:856-863`, `:901-944` | Layer 1 compares the returned sums with `identical()`. `aggregate_results` stays in R. | 3 |
| P10 | Return value shape: names, the names attribute on per-arm vectors, list order, raw vectors only when `return_variance` or `return_percentiles` is on. | `simulation_driver.R:885-898` | Layer 1 `identical()` on the whole list. | 3 |
| P11 | State carried across trials inside a chunk (anything not reset per trial). | (inv) | Harness self-test H2 (trial-by-trial composition of the R body equals the R body) proves there is none outside the RNG. | 2b |
| P12 | Configuration branches the frozen designs never exercise. A port that silently takes a wrong branch there gives wrong numbers to a future caller. | (inv) | `simulate_chunk_cpp` supports exactly the frozen branch set. For any other configuration the dispatcher falls back to the R body, and a test asserts the fallback. | 3, 4 |
| P13 | `tick()` (progress bar) is called per trial in sequential mode. | `simulation_driver.R` (inv) | The C++ path does not tick. With `show_progress = TRUE` the dispatcher uses the R body. | 3 |
| P14 | `EVOLVETRIAL_USE_CPP` is read on every posterior call. FALSE switches the gamma loop to R's `rgamma` (same draws) and `ph_beta_mode_var` to its R fallback, which sums in long double where the C++ solver uses Armadillo double sums, so between-arm results change. | `posterior_cpp_dispatchers.R:11-13, 98-100`, inventory 2.69, 4.30 | The supported-set predicate requires the variable to be TRUE (as frozen). Otherwise the dispatcher uses the R body. Test asserts the fallback. | 3, 4 |
| P15 | Dead draws on the frozen between-arm path. The 800 x K gamma draws at `posterior_helpers.R:130` feed medians the decision never reads. Skipping them shifts the stream. | inventory 6.1, 4.33 | The port draws them in order after the 800 normals. The median computation may be skipped only if nothing downstream reads it. Layer 1 on R23. | 4 |
| P16 | `ph_beta_mode_var_cpp` uses Armadillo reductions, not R's. | inventory 4.30, 6.24 | Call the existing function unchanged. Do not rewrite it with `r_semantics.h`. | 4 |
| P17 | Look times are built by repeated addition (`next_calendar_look + beat`, `post_enroll_time + beat`), not `j * beat`. The comment at `simulation_driver.R:698` says `j * beat`. | inventory 4.8, 4.41 | Copy the repeated addition. Layer 1. | 3 |
| P18 | `mean()` of a logical vector is a long-double count divided by n, which differs from plain `c / n` in about 1 in 8,000 cases. `seq(0, 1, length.out = K + 1)` equals `i * (1 / K)`, not `i / K`. | inventory 4.12, 4.36 | Helpers `mean_logical` and `seq_probs` with unit tests. | 2b |
| P19 | In vs-reference mode the reference arm never leaves "recruiting", so `all_stopped_early` is always FALSE for R23 to R28. `is_eligible` compares `enrolled_counts` with `max_total_patients_per_arm` by position while other code indexes by name. | inventory 6.9, 6.10 | Copy as written, positional where R is positional. The supported-set predicate requires `names(max_total_patients_per_arm) == arm_names`. | 3, 4 |
| P20 | After rebalancing the cutpoints span `[0, max_follow_up_sim]` (120), not the modelling horizon (24). Person-time is clipped at the last cutpoint, but events past it still count in the last interval. | inventory 6.14, 6.15 | Copy as written. Unit test of the whole cutpoint function against `rebalance_cutpoints_from_events`. | 2b, 3 |
| P21 | Bit identity is specific to x86-64 (80-bit long double). | inventory 6.23 | Recorded. Local and Longleaf are both x86-64. | none |

## 2. Steps

### Step 2a. Reference library and recapture of the 10 designs (Opus agent A)

Output: `dev/perf/lib/ref`, `dev/perf/evidence/scripts/capture_cert.R`, `dev/perf/evidence/args/args_<design>.rds` for the 10 designs, `dev/perf/evidence/args/capture_check.csv`, report `dev/perf/evidence/capture_report.md`.

**Task 2a.1. Build `lib/ref`.**
1. `mkdir -p dev/perf/lib/ref dev/perf/build` and add `dev/perf/lib/`, `dev/perf/build/` and `dev/perf/harness/cache/` to `.gitignore`.
2. `bash ~/Downloads/atbo-longleaf/pilot/verify/engine/build_verify_source.sh ~/research/evolveTrial dev/perf/build`. Expected last line: `built p6 (7565a17), m4 (...), m1 (...)`.
3. `R CMD INSTALL --no-multiarch -l dev/perf/lib/ref dev/perf/build/evolveTrial_m4`.
4. Verify: `md5sum dev/perf/build/evolveTrial_m4/R/simulation_driver.R` gives `54763b77b65793392c016d07c14acd7d`, and the installed DESCRIPTION has `VerifyM4PatchMD5: 83de2c8b2de1902630eddf00c334b508` and `P6RecertEvolveTrialSHA: 7565a17c625573d2a697d0d4a85efb2fa6c58ad9`.

**Task 2a.2. Write `capture_cert.R`.** Usage `Rscript capture_cert.R <design> <R> <out.rds>`. It must copy, not paraphrase, the certify runner's set-up from `~/Downloads/atbo-longleaf/pilot/retune/R/retune_run.R`, citing line numbers in comments.
- The `ENV` block (lines 88-93) with `WORKERS = 19` and `SEED` = the design's first cert seed, set before any package load.
- Libraries: evolveTrial from `dev/perf/lib/ref` by `lib.loc`, BATON 0.4.0 from `~/R_libs_frozen/jasa_r1_v040`, tidyverse. Assert both versions and paths.
- Task and CSV rows: Tier 1 designs from `results/retune_t1/cert_design.csv` and `pilot/retune/staged/t1_tasks.tsv`. R04 and R26 from `results/retune_t2a/cert_design.csv` and `t2a_cert_tasks.csv` (find how `retune_run.R` resolves the t2a task's protocol and CSV and copy that).
- F0 sources in the protocol's order (lines 118-123) and the `spec <- switch(PROT, ...)` parse (lines 124-127), run in `~/Downloads/atbo-f0`.
- The assertions at lines 131-166 that apply to certification (weights, workers, fidelity 3000/5000/10000).
- The engine set-up (lines 189-192): first replace `evolveTrial:::run_scenarios` in the namespace by a capture hook that stores `base_args`, `scens` and `...` and then calls the original, then `source` the two switch files and run `verify_m4_install("off", ...)` and `p6r_install("p6", "probe", ...)` with their logs in a temporary directory. The hook therefore sits innermost.
- If `<R>` differs from 10,000, override `fidelity_levels[["high"]]` exactly as the runner does with `TEST_R` (lines 159-160).
- Call `W0(theta, fidelity = "high", seed = <cert seed>)`, then save `list(design, task, scenario_id, trial_type, protocol, csv, theta, cert_seed, base_args, scens, dots, wrapper_res, env = ENV, sessionInfo = sessionInfo())`.
- Assert in the capture: `length(base_args$interval_cutpoints_sim) - 1 == theta$n_intervals`, `identical(base_args$random_dropout, FALSE)`, `identical(base_args$min_follow_up_at_final, 6)`, `isTRUE(getOption("evolveTrial.looks_stop_at_final"))`, and `base_args$num_workers == 19`.

**Task 2a.3. Capture and check.** For each of R01, R02, R03, R04, R23, R24, R25, R26, R27, R28, run `capture_cert.R <design> 10000 dev/perf/evidence/args/args_<design>.rds`, one design at a time (19 workers each). Estimated 3 core-hours, about 20 minutes of wall time. Then write `capture_check.csv` with one row per design: the reference row for (design, cert seed) from `acceptance_set_reference_rows.csv` next to the captured `wrapper_res`, and columns `power_count_ok`, `type1_count_ok`, `EN_null_ok`, `EN_alt_ok` (counts as `round(x * 10000)` equal to the stored integer, EN values equal at 4 decimals).

Expected: 10 of 10 rows pass on all four. If any row fails, stop, report the row and the K and P6 fields of its capture, and do not change anything to make it pass (inventory I6).

**Task 2a.4. Report.** `dev/perf/evidence/capture_report.md`: commands, the table, wall times, where the plan conflicted with the runner.

### Step 2b. Seam, layer 0, layer 1 harness, bisection, floating-point helpers (Sonnet agent B)

Output: the seam in `R/simulation_driver.R`, `src/Makevars`, `src/r_semantics.h`, `src/r_semantics_exports.cpp`, regenerated `RcppExports`, `tests/testthat/test-r-semantics.R`, `dev/perf/harness/chunk_harness.R`, `dev/perf/harness/layer0_seam_check.R`, `dev/perf/harness/run_layer1.R`, `dev/perf/harness/tests/test-harness.R`, report `dev/perf/harness/HARNESS_REPORT.md`.

Agent B can start before 2a finishes. It uses the existing `dev/perf/evidence/args_R01.rds` and `args_R23.rds` for its self-tests (R23 with K = 12 is still a valid configuration for identity testing). The full 10-design run waits for 2a.

**Task 2b.1. Baseline.** Install unmodified HEAD into `dev/perf/lib/dev`. Run the package tests against it, `Rscript -e 'library(evolveTrial, lib.loc = "dev/perf/lib/dev"); testthat::test_dir("tests/testthat", package = "evolveTrial", load_package = "installed")'`, and record pass, fail, skip counts and wall time in the report. These are the baseline counts that later steps must not worsen.

**Task 2b.2. The seam.** In `run_simulation_pure`:
1. Rename the closure at line 468 from `simulate_chunk` to `simulate_chunk_r`. Body unchanged.
2. Directly after `aggregate_results` (after line 944, before `chunk_results <- list()`), add

```r
  simulate_chunk <- simulate_chunk_r

  # Test seam for the C++ engine gate (dev/perf/harness): hand the frame to a probe before any chunk runs.
  chunk_probe <- getOption("evolveTrial.chunk_probe")
  if (is.function(chunk_probe)) return(chunk_probe(environment()))
```

3. Nothing else changes. Confirm that no random draw happens in `run_simulation_pure` before line 946 (inventory), so the early return leaves the stream untouched.

**Task 2b.3. Layer 0.** Write `dev/perf/harness/layer0_seam_check.R <lib> <args.rds> <out.rds>`: loads evolveTrial from `<lib>`, sets `options(evolveTrial.looks_stop_at_final = TRUE)`, sets `base_args$num_simulations <- 200` and `base_args$num_workers <- 4`, calls `evolveTrial::run_scenarios(base_args, scens, parallel = FALSE, seed = 900001, return_variance = TRUE)` and saves the result together with `.Random.seed`. Run it for R01 and R23 against `lib/ref` and `lib/dev` (after 2b.2 is installed). Expected: `identical()` TRUE for both designs, on the table and on `.Random.seed`. This script is rerun after every later change to R code.

**Task 2b.4. Harness core**, `dev/perf/harness/chunk_harness.R`. Functions:

- `load_design(path)`: reads an args RDS and returns it.
- `chunk_frame(design, scenario, R, lib)`: builds `args_i` exactly as `run_scenarios`' `run_one` does (`utils::modifyList(base_args, scens[[scenario]], keep.null = TRUE)`, then `return_percentiles = FALSE`, `percentile_probs = c(0, 0.25, 0.5, 0.75, 0.9, 1.0)`, `return_variance = dots$return_variance`, then keep only `names(formals(run_simulation_pure))`), sets `num_simulations = R`, calls `evolveTrial:::run_simulation_pure` with `options(evolveTrial.chunk_probe = function(env) env)` and returns the frame. It asserts `isTRUE(env$looks_stop_at_final)` and that `env$simulate_chunk_r` is a function (inventory I7).
- `chunk_plan(master_seed, R, workers, n_scen = 2)`: `set.seed(master_seed)`, then for each scenario in order, `chunks <- parallel::splitIndices(R, max(1L, min(workers, R)))`, drop empty chunks, `seeds <- replicate(length(chunks), sample.int(.Machine$integer.max, 1L), simplify = TRUE)`. Returns a list per scenario of (indices, seed) pairs. This mirrors `simulation_driver.R:976-1043` and `run_scenarios:1272` for R at or above the sequential threshold.
- `run_chunk(fn, indices, seed)`: `fn(indices, seed = seed)`, returns `list(result = ..., rng = .Random.seed)`.
- `compare_chunk(env, indices, seed, ref = "simulate_chunk_r", cand = "simulate_chunk_cpp")`: runs both, returns `list(same_result, same_rng, ref, cand)`. `cand` may also be a function, for the self-tests.
- `bisect_chunk(env, indices, seed, ref, cand)`: finds the smallest k such that running `indices[1:k]` differs in `.Random.seed` or in the result. Binary search, each probe a fresh run from `seed`. Returns `list(k, sim_index = indices[k], kind = "rng" | "result" | "both", ref_trial, cand_trial)`, where the trial values are the k-th entries of the raw per-trial vectors and the per-arm sums of the k-trial prefix minus the (k-1)-trial prefix.
- A cache of reference results in `dev/perf/harness/cache/`, keyed by an md5 of `deparse(env$simulate_chunk_r)`, the serialized `args`, the indices and the seed, so that porting iterations only rerun the candidate. Use `tools::md5sum` on a serialized temporary file if `digest` is not installed.
- `run_layer1(designs, master_seeds = c(900001, 900002, 900003), workers = c(16L, 19L), R = 50L, long_chunk = 200L, cand)`: for every design, scenario, master seed and worker count, every chunk of the plan, then one extra chunk `seq_len(long_chunk)` per design and scenario at the first master seed. Returns a data.frame (design, scenario, seed, workers, chunk, n, same_result, same_rng, bisect_k, bisect_kind). On any mismatch it runs `bisect_chunk` for that chunk.

`dev/perf/harness/run_layer1.R [--designs R01,R23] [--cand simulate_chunk_cpp]` wraps `run_layer1`, writes `dev/perf/harness/out/layer1_<timestamp>.csv`, prints the count of chunks and of mismatches, and exits 1 on any mismatch.

**Task 2b.5. Harness self-tests**, `dev/perf/harness/tests/test-harness.R`, run with `Rscript -e 'testthat::test_file("dev/perf/harness/tests/test-harness.R")'` against `lib/dev`. Write each test first and watch it fail before the function exists.

- H1. Reference against itself: `compare_chunk` with `cand = env$simulate_chunk_r` gives `same_result` and `same_rng` TRUE on R01 and R23, one chunk of 20 trials.
- H2. Composition: a candidate that runs the R body one trial at a time (`fn(i, seed = NULL)` per index after one `set.seed(seed)`) and sums the per-arm vectors and the time in trial order, concatenating the raw vectors, is `identical()` to the R body on a 30-trial chunk of R01 and of R23. This proves that no state other than the RNG carries across trials inside a chunk (inventory P11). If it fails, stop and report, because the port's design assumes it.
- H3. Plan emulation: for R01 at R = 200 and 4 workers, the chunk sums from `chunk_plan` plus the R body, aggregated in chunk order, reproduce the `run_scenarios` output columns that are sums divided by R (find the formulas at `simulation_driver.R:1058-1135` and use the same arithmetic). Report any column whose formula you could not reproduce instead of loosening the test.
- H4. Bisection finds a stream fault: a candidate built like H2 that calls `runif(1)` once before trial 7 of a 20-trial chunk. Expected `k == 7`, `kind == "rng"`.
- H5. Bisection finds a value fault: a candidate built like H2 that adds 1 to the trial-12 final N. Expected `k == 12`, `kind == "result"`.
- H6. Option guard: `chunk_frame` errors if `evolveTrial.looks_stop_at_final` is unset.

**Task 2b.6. Floating-point helpers.** `src/r_semantics.h`, a header of inline functions in `namespace evt_rsem`, each copied from the R source of the running R version (cite the file and function in a comment):
- `sum_ld(const double* x, R_xlen_t n)`: R's `rsum` (`src/main/summary.c`), long double accumulator, NA handling as `sum(x)` with `na.rm = FALSE`.
- `mean_ld(const double* x, R_xlen_t n)`: R's `real_mean` (`src/main/summary.c`), long double sum divided by n, then the correction pass when the result is finite.
- `median_r(std::vector<double> x)`: `stats::median.default` with `na.rm = FALSE`, odd n the middle order statistic, even n the `mean` (by `mean_ld`) of the two middle values.
- `quantile7(std::vector<double> x, double p)`: `stats::quantile.default` type 7 for one probability, from the R source of the installed `stats`, including how `index`, `lo`, `hi` and `h` are formed.
- `find_interval(const double* breaks, int nb, double x, bool rightmost_closed, bool all_inside, bool left_open)`: R's `findInterval` (`src/appl/interv.c`), with the flag values the engine uses (read them at the call sites in `state_management.R`).
- `mean_logical(int count, R_xlen_t n)`: `mean()` of a logical vector, `(double)((long double)count / (long double)n)` (inventory 4.36). Check the R source to confirm this is what `real_mean`/`mean` does for LGLSXP.
- `seq_probs(int K)`: `seq(0, 1, length.out = K + 1)`, which the inventory found equals `i * (1.0 / K)` (inventory 4.12). Copy from `seq.default`'s source rather than from the inventory's observation, and test K from 1 to 40.
- `rebalance_cutpoints(std::vector<double> event_times, double max_follow_up, int K)`: the whole of `evolveTrial:::rebalance_cutpoints_from_events` (`simulation_driver.R:1331-1398`, inventory 4.11 to 4.16), built from the helpers above, returning an empty vector where the R function returns NULL.
- Nothing else. Slice metrics and the posterior are port work (step 3), not helpers.

Inventory Section 4 (`dev/perf/evidence/chunk_body_inventory.md:252-328`) lists every floating-point site with its R semantics. Read it before writing the helpers.

`src/r_semantics_exports.cpp` exports one thin wrapper per helper (`rsem_sum`, `rsem_mean`, `rsem_median`, `rsem_quantile7`, `rsem_find_interval`, ...) with `// [[Rcpp::export(rng = false)]]`. `src/Makevars` is the single line `PKG_CXXFLAGS = -ffp-contract=off`. Run `Rcpp::compileAttributes()` and check that the diff of `R/RcppExports.R` and `src/RcppExports.cpp` is additions only (inventory B3).

**Task 2b.7. Unit tests**, `tests/testthat/test-r-semantics.R`, written before the helpers. Each helper is compared with `identical()` to the R function on at least 2,000 random inputs per helper, seeded. Inputs cover lengths 0, 1, 2, 3, odd and even, and 10^5. Values cover exponential and Weibull samples like the engine's event times, values at 1e8 scale with small offsets (where long double matters), decimal fractions such as multiples of 0.1, ties, sorted and reverse-sorted input, and for `find_interval` points exactly on the breaks and outside the range. Also a direct long-double test: `x <- c(1, rep(1e-16, 10))` must give `sum(x)` and `mean(x)` exactly as R does, which a double accumulator would not. `mean_logical` is tested on all (count, n) with n up to 2,000, `rebalance_cutpoints` on 2,000 random event-time vectors of length 50 to 300 drawn like the engine's (Weibull shape 1.3, medians 3.8 to 9, some beyond 120), with K from 3 to 16, plus vectors with ties and with values at 0 and at the cap.

**Task 2b.8. Report.** `dev/perf/harness/HARNESS_REPORT.md`: baseline test counts, layer 0 result, self-test results, the unit test counts per helper, the R source locations each helper was copied from, and every place the plan conflicted with the code.

### Coordinator verification for step 2 (before the commit)

1. `md5sum` of the ref tree driver and the DESCRIPTION fields (2a.1).
2. Count rows in `capture_check.csv`, expect 10 with all four checks TRUE. Spot-check one between-arm capture by hand: K, `random_dropout`, `min_follow_up_at_final`, `num_workers`.
3. Rerun layer 0 for R01 and R23 myself.
4. Rerun the harness self-tests and `test-r-semantics.R` myself and count the expectations.
5. Read the seam diff in `R/simulation_driver.R` (must be the rename plus the five lines) and the RcppExports diff (additions only).
6. Run `run_layer1.R` with `cand = simulate_chunk_r` on all 10 designs. Expected zero mismatches. This is the dry run of the gate and gives the chunk count for the record.
7. Commit: "Layer 1 harness, seam, floating-point helpers; layer 0 identical on R01 and R23; captures reproduce 10/10 reference rows".

### Step 3. Single-arm port (Opus agent C)

Brief written by the coordinator at the start of the step, from this plan, the design and the inventory. Shape:

- `simulate_chunk_cpp(sim_indices, cfg, seed_is_set)` in `src/simulate_chunk.cpp`, where `cfg` is a list the R side packs once per frame from the closure's free variables (inventory Section 1). In R, `simulate_chunk_cpp <- function(sim_indices, seed = NULL, tick = function() {}) { if (!is.null(seed)) set.seed(seed); .simulate_chunk_cpp(sim_indices, cfg) }` is defined in the frame next to `simulate_chunk_r`. `set.seed` stays in R.
- `simulate_chunk` becomes a dispatcher: the C++ body when `getOption("evolveTrial.use_cpp_engine", TRUE)` is TRUE and the configuration is in the supported set (single-arm only in this step), otherwise the R body.
- Tasks in order: config packing and the supported-set predicate with a fallback test; registry arrays and the arrival loop with no looks (layer 1 on a design with looks disabled is not available, so the first gate is the RNG state after trial 1 of R01); slicing and metrics through `r_semantics.h`; the interim look with gates and the gamma posterior; rebalancing; post-enrollment looks under `looks_stop_at_final`; the final analysis; accumulators and the return list.
- Gate: `run_layer1.R --designs R01,R02,R03,R04` zero mismatches, layer 0 identical with the option FALSE, package tests no worse than baseline. Coordinator reruns the gate, then times R01 at 1 worker, 1,000 trials, against `lib/ref`, and commits with both results in the message.

### Step 4. Between-arm port (Opus agent D)

Same shape. The vs-reference look, the C++ copy of the live R PH sampler (`R::rnorm(mean, sd)` for S = `num_posterior_draws`, then for each sample `R::rgamma` over the K intervals, in the order of `posterior_helpers.R:124-133`), the existing `ph_beta_mode_var_cpp` called directly, and the final vs-reference analysis. The supported set gains `compare_arms_option = TRUE` with `use_ph_model_vs_ref = TRUE`. Gate: `run_layer1.R` on all 10 designs zero mismatches, then timing for R23. Commit with both.

### Step 5. Layer 2 on Longleaf (coordinator)

Install `lib/dev` on Longleaf with `module load gcc/12.2.0 r/4.5.0`. A Slurm array of 75 tasks, one per reference row, each with `--cpus-per-task=19`, running `capture_cert.R`'s call path at R = 10,000 with the C++ engine on. Compare with `acceptance_set_reference_rows.csv` on `power_count`, `type1_count` (integers) and `EN_null`, `EN_alt` (doubles at the stored precision). Pass means 75 of 75. Estimated 3 to 5 core-hours after the port. Before submitting, the coordinator states the core-hour estimate to the PI.

### Step 6. Layer 3, search-side (coordinator)

Two or three E2c BATON evaluations from the stored ledgers at R = 3,000 and 16 workers from their recorded seeds, through the C++ engine, compared with the ledger values. Under 1 core-hour.

### Step 7. Readout and handoff

Timing table (R01 and R23 at 1 and 16 workers, ref vs dev), gate results, `dev/HANDOFF-<date>.md`, commit, `/clear`.
