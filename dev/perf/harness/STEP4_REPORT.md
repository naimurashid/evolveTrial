# Step 4 report: between-arm C++ port of `simulate_chunk` (agent D, 2026-11-03)

Branch `perf/engine-cpp`, worktree `~/research/evolveTrial-perf`, on top of 0870256. Nothing committed or pushed.
Installed only into `dev/perf/lib/dev` with
`R CMD INSTALL --no-multiarch --preclean -l dev/perf/lib/dev . && git checkout -- src/symbols.rds`.

## Gate results

| Gate | Command | Result |
|---|---|---|
| 1. Layer 1 | `Rscript dev/perf/harness/run_layer1.R --cand simulate_chunk_cpp` (all 10 designs) | 2,120 chunks compared, 0 mismatches (0 result, 0 stream), 293 s. 212 chunks per design, 0 mismatches in each of R01 to R04 and R23 to R28. CSV `dev/perf/harness/out/layer1_20261103_141657.csv`, log `dev/perf/harness/out/layer1_step4_gate.log`. Exit 0. |
| 2. Layer 0, option FALSE | `layer0_seam_check.R dev/perf/lib/<ref,dev> args_<D>.rds <out> FALSE`, then `layer0_compare.R` | R01: table identical TRUE, `.Random.seed` identical TRUE (2 x 23). R23: TRUE, TRUE (4 x 23). |
| 2. Layer 0, option TRUE | same with `TRUE` | R01: TRUE, TRUE. R23: TRUE, TRUE, now through the C++ body (lib/dev 1.0 s against lib/ref 9.0 s at 200 trials, 4 workers). Outputs and log in `dev/perf/harness/out/l0s4/` (`layer0_step4.log`). |
| 3. C++ actually ran | `Rscript dev/perf/harness/cpp_entry_check.R layer1 R23,R24,R25,R26,R27,R28` and `... layer0 dev/perf/lib/dev args_<D>.rds <TRUE/FALSE> l0s4/dev_<D>_<mode>.rds` | Layer 1: 1,272 entries into `.simulate_chunk_cpp` for 1,272 between-arm chunks, 0 mismatches. Layer 0: R23 TRUE 8 entries (2 scenarios x 4 chunks), R23 FALSE 0, R01 TRUE 8, R01 FALSE 0. Each traced run identical to its gate output. Log `dev/perf/harness/out/cpp_entry_check_step4.log`. The script needed no change. |
| 4. Package tests | `Rscript dev/perf/harness/run_pkg_tests.R dev/perf/lib/dev dev/perf/harness/out/step4_tests.rds` | 723 passed, 0 failed, 33 skipped, 0 errors, 0 warnings. That is 495 (step 3) plus 228 in the new file. The step 3 file still has 148 expectations; one of them changed (listed below). Log `dev/perf/harness/out/step4_tests.log`. |
| 5. Harness self-tests | `testthat::test_file("dev/perf/harness/tests/test-harness.R")` against lib/dev | 26 passed, 0 failed, 0 skipped. Log `dev/perf/harness/out/harness_selftest_step4.log`. |
| 6. RcppExports (B3) | `Rcpp::compileAttributes(".")`, `git diff --numstat` | `R/RcppExports.R` 4 added, 0 deleted. `src/RcppExports.cpp` 23 added, 0 deleted. The additions are one R wrapper (`.chunk_ph_sample`), one C++ wrapper with `Rcpp::RNGScope rcpp_rngScope_gen`, and one registration line. No existing wrapper changed. |
| 7. Tree | `git status`, `git diff --quiet src/symbols.rds` | `src/symbols.rds` unchanged. No `.o` or `.so` tracked (0). The 8 in `src/` are ignored by `src/.gitignore`. |

## R23 timing

The run used 1 worker, 1,000 trials, both scenarios, seed 900001, through `run_scenarios` with `return_variance = TRUE`.
The script is in the scratchpad (`timing.R`) and is not in the repo. Each library ran once, one after the other, on an
idle machine.

| Library | Wall time (2 scenarios) | Per scenario |
|---|---|---|
| lib/ref | 163.28 s | 81.6 s |
| lib/dev | 18.75 s | 9.4 s |

The speed-up is 8.71x, and the table and `.Random.seed` are `identical()`. lib/dev is above the 8 s target per 1,000
trials, by scenario: about 5 s for scenario 1 (null, both medians 6) and about 14 s for scenario 2 (Triplet median 9).
Scenario 2 trials stop for futility less often, so they run about three times as many PH looks.

**Where the time goes.** perf is unavailable (`perf_event_paranoid` = 4), so I counted draws instead. I traced the R
body's `sample_vs_ref_medians_ph` on a 200-trial R23 chunk, then timed R's vectorized `rgamma` and `rnorm` on the same
number of draws with the same shapes and scales. That is a lower bound on the C++ path's draw cost.

| Scenario | PH calls per trial | Gamma draws per 200 trials | C++ chunk | Draws alone | Share in draws |
|---|---|---|---|---|---|
| 1 | 10.3 | 24.7 M | 0.99 s | 0.92 s | 93% |
| 2 | 28.8 | 69.1 M | 2.77 s | 2.55 s | 92% |

About 92% of the C++ time is the 800 x K gamma draws themselves. Those draws are dead for the decision (P15), but they
must be made in order. The rest of the port (slicing, the Newton solver, medians, bookkeeping) is under 8%. The 8 s
target cannot be reached without removing or changing draws, which phase 1 forbids. This is the design's "draw count as
fidelity" item for phase 2. I changed no draw to gain speed.

## Files changed and why

| File | Change |
|---|---|
| `R/simulate_chunk_cpp.R` | Predicate restructured into `.chunk_cpp_common_checks`, `.chunk_cpp_hc_checks` and `.chunk_cpp_vsref_checks`. `.chunk_cpp_static_checks` picks the HC list or the vs-ref list by the frame's `compare_arms_option`. Packing split into a common part plus `.pack_hc` (step 3 fields, unchanged expressions) or `.pack_vsref`. New `cfg$mode` ("hc" or "vsref"). Test hook `.chunk_ph_sample_r`. |
| `src/simulate_chunk.cpp` | Cfg and unpack for mode "vsref". New `median_piecewise_r`, `ph_sample_core` and `ph_sample` (the PH sampler), `logical_mean`, `gates_pass_vsref`, `interim_vsref` and `trt_can_enroll`. The two vs-ref accrual breaks (SD:561-570, SD:640-646), look dispatch in the accrual and post-enrollment loops, and the vs-ref final (SD:736-833). Exported test hook `.chunk_ph_sample`. The HC code paths are unchanged apart from the HC final loop now being skipped in vsref mode. |
| `R/RcppExports.R`, `src/RcppExports.cpp` | Regenerated, additions only (gate 6). |
| `tests/testthat/test-simulate-chunk-cpp-vsref.R` (new) | 228 expectations, below. |
| `tests/testthat/test-simulate-chunk-cpp.R` | One step 3 test changed, below. |

`src/posterior_sampling.cpp`, `src/r_semantics.h`, `R/simulation_driver.R` and the harness are unchanged. `grep -n
'arma::randn\|sample_vs_ref_medians_ph_cpp' src/simulate_chunk.cpp R/simulate_chunk_cpp.R` returns nothing (exit 1).

## The PH posterior (P5, P15, P16)

- **The solver.** `ph_beta_mode_var_cpp` is called directly and unchanged, with every argument passed explicitly. It
  has exactly one definition, at `src/posterior_sampling.cpp:196`. `src/hybrid_trial_sim.cpp:88` is a forward declaration
  ending in `);` with no body, not a second function. `nm -C` on the installed `evolveTrial.so` shows one global `T
  ph_beta_mode_var_cpp(...)` symbol, and only `src/posterior_sampling.o` defines it. The R dispatcher reaches it through
  `R/RcppExports.R:103` -> `_evolveTrial_ph_beta_mode_var_cpp` (`src/RcppExports.cpp:98`) -> that definition.
- **Coercions, as in PD:102-113.** `E_C` and `E_T` become `arma::ivec` (the dispatcher's `as.integer`), and `PT_C`,
  `PT_T`, `alpha_prior` and `beta_prior` become `arma::vec`. `mu` is `coalesce_num(args$ph_loghr_prior_mean, 0)`.
  `sigma` is `max(1e-6, coalesce_num(args$ph_loghr_prior_sd, 1))`. `tol` and `max_iter` are
  `eval(formals(ph_beta_mode_var)$tol)` and `as.integer(eval(formals(ph_beta_mode_var)$max_iter))`. All of these are
  computed in R once per frame, so the C++ receives R's own parsed doubles.
- **Draw order.** First `num_samples` `R::rnorm(mean, sd)` draws, in index order. Then, for each sample, `exp(beta_i)`,
  the rate `(beta_prior[k] + PT_C[k]) + exp_beta * PT_T[k]`, `scale = 1 / rate`, and `R::rgamma(shape_k, scale)` over k
  in interval order, with `shape_k = (alpha[k] + E_C[k]) + E_T[k]`. `num_samples` is `args$num_posterior_draws` at the
  interim (IL:7) and the frame's `num_posterior_draws` at the final (SD:802). The predicate requires the two to be
  identical.
- **The medians are computed but not stored.** Nothing downstream reads `medCtrl` or `medTrt` on the supported path.
  `grep` over `R/` shows them read only at IL:34 (the non-PH branch of `calculate_current_probs_vs_ref`), SD:820 (the
  non-PH final) and `R/predictive_probabilities.R:178, 388`, which is not on the chunk path. `return_percentiles` only
  affects `final_n_raw`. I still evaluate both medians per sample with a port of `calculate_median_survival_piecewise`,
  because R evaluates them and their `if()` conditions error on NaN. The port stops where R would stop. The cost is
  inside the under-8% non-draw share.
- **Unit test before wiring.** The test traces the live `sample_vs_ref_medians_ph` through 6 R23 trials per scenario.
  The captured calls cover interim and final calls, with both original and rebalanced cutpoints. Each captured call is
  replayed from its entry `.Random.seed` through the R function and through `.chunk_ph_sample_r`. The comparison is
  `identical()` on the whole return list (`medCtrl`, `medTrt`, `logHR`) and on `.Random.seed` afterwards: 0
  differences. Empty control, empty treatment and both-empty slices (double zero vectors) are also identical.

## The supported-set predicate

A frame whose `compare_arms_option` is not TRUE is checked against the step 3 list. The test "the HC frame is checked
against exactly the step 3 items" asserts the 18 labels and their order on the R01 frame. A frame whose
`compare_arms_option` is TRUE is checked against the common items plus the items below. Every check stays wrapped in
`tryCatch` as in step 3. All 12 between-arm frames (R23 to R28, both scenarios, built by `chunk_frame`) satisfy every
item and pack as supported. Layer 1 ran the C++ path on all of them.

| Item | Check | Inventory rows |
|---|---|---|
| common 4, 5, 6, 7, 8, 12, 13b, 13g | unchanged step 3 functions. 13g was already by arm name over all of `arm_names`, so it covers two arms as written. | 2.5, 2.7, 2.11, 2.14, 2.18, 2.4, P13, 3.5, 4.1-4.7 |
| 11 | `EVOLVETRIAL_USE_CPP` truthy at dispatch (unchanged). It also switches `ph_beta_mode_var` to its long-double R fallback (P14). | 2.69, 4.30 |
| 13e (vs-ref form) | `num_posterior_draws` one finite number, at least 1, identical in args | 6.8 |
| V1 | `isTRUE(compare_arms_option)` in the frame and in args | 2.35 |
| V2 | `isTRUE(args$use_ph_model_vs_ref)` | 2.29, 2.42, 2.61, 2.39 |
| V3 | two distinct non-NA arm names, identical in args; `reference_arm_name` one of them, identical in args; exactly one experimental arm (refuses P6) | 2.36, 2.37, 2.57 |
| V4 | `identical(names(max_total_patients_per_arm), arm_names)`, finite, non-negative, identical in args | P19, 6.10, 2.58 |
| V5 | `!isTRUE(efficacy_stopping_rule_vs_ref)` in args and in the frame | 2.44 |
| V6 | `args$compare_arms_hr_margin` one finite number, `args$compare_arms_hr_eff_margin` NULL | 2.30, 2.31, 2.43 |
| V7 | both probability thresholds one finite number in args (interim) and in the frame (final) | 2.32, 2.45 |
| V8 | `compare_arms_target_arm` not in `names(args)`, and `args$compare_arms_target_arm` (with partial matching, as IL:105 reads it) NULL | 2.38 |
| V9 | `min_events_per_arm`, `min_median_followup_per_arm` and `min_person_time_frac_per_arm` each NULL or one unnamed finite number. The captures carry 11L to 22L, 3 and 0. | 2.60, 6.7 |
| V10 | `ph_loghr_prior_mean`, `ph_loghr_prior_sd` NULL or one unnamed finite number. The captures carry 0 and 1. | 2.63 |
| V11a | `num_posterior_draws` integral and at most `.Machine$integer.max` (it feeds `rnorm(n)` and `seq_len`) | PH:124, 127 |

For V11 I went through the inventory rows marked NOT EXERCISED that a vs-ref configuration can reach. The non-PH
branches 2.29, 2.39, 2.42 and 2.61 are refused by V2. 2.30 (NULL margin) and 2.31 (non-NULL efficacy margin) are refused
by V6. The 2.44 TRUE branch is refused by V5. The 2.60 named and positional gate branches are refused by V9, and the
2.38 override by V8. 2.36 and the stops in `gates_pass_for_both_arms` are refused by V3 and V4. The other rows are on the
generic code, which is ported: 2.37 (no recruiting experimental arm), 2.59 (`maxPT > 0` both ways), and both orders of
the reference arm. The test "reference arm listed second" runs that case and is identical. No row that R23 to R28
exercise is missing from the list.

The gate thresholds, the PT caps, `log_margin` (IL:18 and SD:808 expressions, computed separately) and `log_eff` (0) are
computed in R by the package's own `resolve_gate_vec` and the same expressions, once per frame. The interim futility
threshold comes from args and the final thresholds from the frame, kept apart as in step 3's C1.

### Tests added (`tests/testthat/test-simulate-chunk-cpp-vsref.R`, 228 expectations)

- PH sampler identity on traced R23 calls and on empty slices (above).
- R23 chunk identity (trials 3 to 12), both scenarios, under `return_variance` x `return_percentiles` (4
  combinations): `identical()` on the result and on `.Random.seed`, and `cfg$mode == "vsref"`.
- The sequential path (`seed = NULL`) on R23 is identical.
- The reference arm listed second (`ref_idx` 1) is supported and identical.
- Label sets: the HC list equals the 18 step 3 labels in order. The vs-ref list starts with the 8 common labels and
  contains V1 to V10 and V11a. The R23 frame has no violation.
- 27 modified copies of the R23 frame, covering V1 to V11a, 13e and the common items 4 to 8, 12, 13b and 13g, with two
  or three variants for V3, V4, V6, V7, V9 and 13g. For each copy: the targeted label is among the violations,
  `cfg$supported` is FALSE, the predicate is FALSE, `simulate_chunk_cpp` stops with "outside the supported set",
  `.Random.seed` is unchanged after the stop, and the dispatcher calls the R body and not the C++ body. The unmodified
  copy dispatches to C++, and to R with `evolveTrial.use_cpp_engine = FALSE`. With `EVOLVETRIAL_USE_CPP=FALSE` the
  dispatcher takes R and `simulate_chunk_cpp` stops.

### Step 3 test changed (item 3 of the brief)

In `test-simulate-chunk-cpp.R`, the test "an unsupported frame (between-arm) packs without error and falls back" is
renamed to "an unsupported between-arm frame (no HR margin) packs without error and falls back". Its one expectation
`expect_true(any(grepl("^1 ", env$cfg$violations)))` became `expect_true(any(grepl("^V6 ", ...)))`. The frame leaves
`compare_arms_hr_margin` at its NULL default, so it now fails V6 instead of item 1, and it still falls back with an
identical result and stream. No other step 3 expectation changed. The "1" modification in the step 3 predicate test (an
R01 frame with `compare_arms_option` TRUE) is now refused by the vs-ref items (V2, V3 and others) and still passes
unchanged. The R23 candidate is no longer refused anywhere: it is supported, and layer 1 and layer 0 run it through C++.

## Order of work and bisection episodes

I wrote the predicate, the packing, the PH sampler and the chunk path in one pass, as step 3 did, instead of gating
each piece separately as brief Section 4 orders. The first check after the build was `.Random.seed` and the result after
trial 1 of R23, then 5 trials, on both scenarios. All were identical. 60-trial chunks on all six between-arm designs
and on R01 were identical too. I wrote the tests after the code, including the PH unit test, which therefore did not
come first as the brief orders. It passed. The full layer 1 run passed on its first execution, so the bisection tool
never ran: there were no bisection episodes.

## Conflicts between the brief, plan, design or inventory and the code

- **D1 (brief Section 1).** "`src/hybrid_trial_sim.cpp:88` has a second function of the same name." It is a
  declaration of the one function defined at `posterior_sampling.cpp:196`, not a second definition (evidence above).
- **D2 (brief, frozen configuration).** `compare_arms_futility_margin` is 0 in the frame and in args for all 12
  between-arm frames, not NULL. It is read only on the non-PH path and in diagnostics messages, so the port does not use
  it and the predicate does not pin it.
- **D3 (brief, frozen configuration).** `num_posterior_draws_interim` is NULL in the frame, but `args` carries 800,
  because SD:370-374 fills it from `num_draws_interim`. It is not read on the vs-ref path (inventory 6.8).
- **D4 (brief "Read first" item 5).** `R/state_management.R:150-170` is not the reference and arm slicing. That is
  `slice_arm_data_at_time` at SM:115-149, ported in step 3 and reused. Lines 150-170 are the start of
  `gates_pass_for_both_arms`. The inventory row 2.20 points to SD:737. These are line references only.
- **D5 (brief Section 3, "13e for `num_posterior_draws`" as a common item).** Step 3's 13e checks both draw counts. To
  leave the HC set exactly as in step 3, I kept step 3's 13e on the HC path and gave the vs-ref path its own "13e
  num_posterior_draws" label. I also added checks the brief's list did not name, all of which narrow the set: V3 also
  checks the args copies of `arm_names` and `reference_arm_name` (the interim reads args, the body reads the frame), V7
  checks both copies of the thresholds, and V11a requires an integral draw count.
- **D6 (design Section 2 vs plan P6 and the brief).** The design says the control-cache draws spent when more than one
  experimental arm is evaluated "are kept". Plan P6 and the brief say the port refuses that configuration. I followed the
  brief: V2 and V3 refuse it.
- **D7 (design Section 1 target vs bit identity).** The target of under 8 s per 1,000 trials is not met: 9.4 s per
  scenario on average, about 14 s on scenario 2. About 92% of the C++ time is the gamma draws that bit identity requires.
  lib/ref measured 81.6 s per scenario against the design's 75.1 s baseline.
- **D8 (inventory 1.4 and 1.18).** "R23 K = 12 (capture)" is stale. The step 2a frames have K = 15 for R23 and R28, 9
  for R24, 12 for R25, 13 for R26 and 14 for R27. This was already known (I1, 6.20). The port does not hard-code K.
- **D9 (tooling).** The harness's reference-cache key hashes every `R/*.R` file except `simulation_driver.R` and
  `RcppExports.R`, including the port's own `R/simulate_chunk_cpp.R`. Every edit to that file therefore invalidates the
  cached R-body references. This is harmless (the gate rebuilt them in under 5 minutes). I did not change the harness.

## Open items and anything left unported

- **Step 3 O1 now covers R23 to R28 too.** On the vs-ref path the C++ engine never calls `interim_check`, so the certify
  capture path's probe-based engine check would see no probe lines for between-arm designs either. Stage A's instrument
  rewrite would hit the dispatcher. As in step 3, I changed nothing. The coordinator decides how layer 2 asserts its
  environment.
- **Unported (falls back to the R body through the dispatcher):** the vs-ref independent (non-PH) model; more than one
  experimental arm (P6); `compare_arms_target_arm`; interim vs-ref efficacy stopping; a non-NULL
  `compare_arms_hr_eff_margin`; a NULL `compare_arms_hr_margin`; named or vector vs-ref gates; non-integral draw counts.
  The step 3 list also still falls back: person-time milestones, `looks_stop_at_final = FALSE`, dropout on, infinite
  `max_trial_time`, diagnostics messages and the progress tick.
