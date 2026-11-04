# Default switch to the exact engine, and engine-equivalence tests (2026-11-04)

Worktree `~/research/evolveTrial-perf`, branch `perf/engine-cpp`, HEAD 3d5e776. Nothing committed, pushed or installed
into `~/R/library`. Test libraries are `dev/perf/lib/ds_before` (HEAD) and `dev/perf/lib/ds_after` (working tree),
built by `dev/perf/phase2/ds_install.sh` and run with `dev/perf/harness/run_pkg_tests.R`, the G4 method.

## Counts

| Item | Result |
|---|---|
| Engine-option read sites | 2 found (`R/simulate_chunk_exact.R`, in `.exact_engine_selected` and `.exact_requested`), 2 changed, plus 1 new read in a new helper `.exact_engine_explicit`. No other site reads `evolveTrial.engine` in R/. |
| Suite before (HEAD) | 1013 passed, 0 failed, 33 skipped, 0 errors, 62 s (matches G4) |
| Suite after the switch, before the new file | 1013 passed, 0 failed, 33 skipped (`ds_pkg_tests_mid.log`) |
| Suite after, with the new file | 1186 passed, 0 failed, 35 skipped, 0 errors, 106 s (`ds_pkg_tests_after.log`). The difference is the new file: 173 passed, 2 skipped (the two long tests) |
| New file runtime, routine | 42 s (largest pieces: the R body against exact on OCs 16 s, fallback checks 9.5 s, cpp against the R body 7.6 s, golden values 7.3 s) |
| New file, `EVOLVETRIAL_LONG_TESTS=1` | 210 passed, 0 failed, 0 skipped, 487 s (`ds_eq_long.log`). The long OC test took 350 s and the long cpp/R-body identity test (including the 2-worker parallel path) took 95 s |
| DESCRIPTION Version | 0.0.0.9001 to 0.0.0.9002 |

The man page was edited after the "after" suite ran. That edit changed text only (roxygen comment and Rd), not code.

## Task A, the code change

- The default of `getOption("evolveTrial.engine", ...)` in `.exact_engine_selected` is now `"exact"`.
- `.exact_requested`, which gates the L2 notice, now reads the option with no default. An unset option therefore never
  fires the notice. The hybrid (`run_hybrid_simulations_rcpp`) and binary (`run_simulation_binary`) paths call
  `.exact_note_unused`, so they stay silent under the default.
- In the driver (`R/simulation_driver.R`, near line 990), `.exact_engine_explicit()` (the option is non-NULL) now decides
  what happens when the engine is selected but the frame is unsupported (`.chunk_cpp_supported(env) && cfg$supported`
  is FALSE). With the option unset, `exact_engine` becomes FALSE and the phase 1 dispatcher runs: the C++ body if it is
  supported, otherwise the R body. That is exactly what `"cpp"` does, and nothing is drawn. With
  `options(evolveTrial.engine = "exact")` set explicitly, the existing stop is unchanged. The condition is the exact
  negation of `.chunk_cpp_assert_supported`'s failure condition, so under the default that stop can never fire.
- The distinction between "explicit" and "default" was possible: `getOption()` returns NULL when the option is unset.
  One consequence is that `options(evolveTrial.engine = NULL)` means "default".
- "cpp" is unchanged and stays bit-identical to the R body. The routine and long tests both check this end to end.
- Documentation: a new roxygen `@section Simulation engine` on `run_simulation_pure`, with the Rd edited by hand to
  match. roxygen2 8.0.0 is installed against RoxygenNote 7.3.3, and running `document()` would have rewritten unrelated
  man files. The header comments in `R/simulate_chunk_exact.R` and the driver are updated, and NEWS.md has a new
  0.0.0.9002 section.
- Existing tests: none failed after the switch, because every existing test that runs an engine either sets
  `evolveTrial.engine` explicitly or works at frame level through the `chunk_probe` seam, which returns before the
  engine choice. The only edit was a stale comment in the L2 test (`test-simulate-chunk-exact.R:352`).

## Where the brief conflicts with the code

1. **Under every package default the exact engine does not run.** The supported set requires
   `options(evolveTrial.looks_stop_at_final = TRUE)` (check item 4), and that option defaults to `FALSE`
   (`simulation_driver.R:251`). So a call with all options at their defaults still runs the R body, unchanged. The
   switch takes effect only for callers that set `looks_stop_at_final = TRUE`, as the gate harness and the tests do.
   I did not change that default because the brief did not cover it, and NEWS and the man page say this. If "exact by
   default" is meant to reach plain calls, flipping `looks_stop_at_final` is a separate ruling with its own semantic
   change.
2. **The brief's "seamless" path is not separate from hybrid.** Grep finds no separate seamless simulator that calls
   the notice. Seamless runs through `run_hybrid_simulations_rcpp`, which is covered.
3. **The brief described the existing posterior tests as comparing against the original R samplers.** They compare
   against the C++ samplers (`draw_posterior_hazard_samples_cpp`, `.chunk_ph_sample`). The new tests add the comparison
   against the R samplers.

## Task B, `tests/testthat/test-engine-equivalence.R` (11 tests)

Fixtures: the captures `args_R01.rds`, `args_R03.rds` and `args_R23.rds` are copied into `tests/testthat/fixtures/`
(about 95 KB), so the file does not depend on `dev/`. `fixtures/engine_golden.rds` (17 KB) holds the golden values.

1. **cpp against the R body, end to end through `run_scenarios`.** R01 and R23, both scenarios, `return_variance`
   TRUE and FALSE. The table and `.Random.seed` must be `identical()`. Tests that already cover this at chunk level, and
   are not repeated: in `test-simulate-chunk-cpp.R`, R01 with return_variance TRUE/FALSE and percentiles, plus the
   sequential path; in `test-simulate-chunk-cpp-vsref.R`, R23, plus the reference arm listed second. Nothing tested
   identity through the full driver and aggregation before. The long version adds R01 at 600 trials, R23 at 200 and
   the 2-worker parallel path.
2. **Exact posteriors against the original R samplers.** 40,000 draws from `draw_posterior_hazard_samples_r` with
   `calculate_median_survival_piecewise`, over 6 single-arm states and 16 thresholds. The states cover a single
   interval with a weak prior and with a strong prior, R01-like 8-interval states at an early and a late look, a
   zero-event state where most posterior mass lies beyond the grid, and a large prior with unequal interval lengths.
   For between-arm PH, `ph_beta_mode_var_r` with 100,000 `rnorm` draws over 4 states (sigma 1 and 0.2, a large hazard
   prior, a nonzero prior mean) and 5 margins, plus the mode and SD agreeing to 1e-6. One R23 state uses the full
   `sample_vs_ref_medians_ph_r` (20,000 draws). Each check requires |z| < 4 with z = (MC - p) / sqrt(max(p(1-p), 1/n)/n),
   at fixed seeds.
3. **Operating characteristics, exact against the R body.** Power or type I error, PET_Futility and Exp_N per arm and
   scenario, with z = (X_exact - X_R) / sqrt(Var_exact + Var_R). Var_* is the Monte Carlo variance of the estimate
   that `run_scenarios` returns, and the check is |z| < 4. Metrics with zero variance in both runs must be equal.
   One-sided checks: power_exact >= power_R - 4 SE, and PET_Futility_exact <= PET_Futility_R + 4 SE. Routine size:
   R01 at 4,000 exact trials against 200 R-body trials, R23 at 2,000 against 60. Long size: R01 and R03 at 40,000
   against 4,000, R23 at 20,000 against 1,000. All pass.
4. **Golden values.** Six entries (exact and cpp on R01, R03 and R23, fixed seeds, R03 with return_variance FALSE)
   store the full table and `.Random.seed`. The table is compared at tolerance 1e-10 (one changed trial moves a
   proportion by 1/n). The RNG state must be identical. To regenerate: `EVOLVETRIAL_UPDATE_GOLDEN=1`. The values were
   generated from engine code unchanged since 669ca3e.
5. **Default selection.** With the option unset, on a supported frame:
   - the output is identical to explicit "exact" and differs from "cpp";
   - exactly one `sample.int(.Machine$integer.max, 2)` draw is made per scenario;
   - no message is printed.

   On unsupported frames, the default output is identical to explicit "cpp" (table and RNG) and silent. The unsupported
   frames are `looks_stop_at_final` FALSE on R01 and on R23, R23 without its HR margin, `use_cpp_engine` FALSE, and
   `EVOLVETRIAL_USE_CPP` FALSE. An explicit "exact" still errors on an unsupported frame. The L2 notice is silent under
   the default, including through `run_simulation_binary`, and fires when the option is explicit.

## Files

Changed: `DESCRIPTION`, `NEWS.md`, `R/simulate_chunk_exact.R`, `R/simulation_driver.R`, `man/run_simulation_pure.Rd`,
`tests/testthat/test-simulate-chunk-exact.R`. New: `tests/testthat/test-engine-equivalence.R`,
`tests/testthat/fixtures/{args_R01,args_R03,args_R23,engine_golden}.rds`, `dev/perf/phase2/ds_install.sh`, and logs and
RDS files `dev/perf/phase2/ds_*`. Scratch: `dev/perf/build/ds_*`, `dev/perf/lib/ds_*` (deletable).
