# evolveTrial 0.0.0.9005

## Fast engines cover bidirectional designs and the symmetric rule variant (phase 3, step P3-2)

Both C++ engines (the draw-based `"cpp"` engine and the default exact engine)
now run interim efficacy stopping and the between-arm hazard-ratio margins,
which used to fall back to the R body.

* **`efficacy_stopping_rule_hc = TRUE`** (single-arm). At each interim look the
  efficacy probability `Pr(median > success value)` is checked before the
  futility probability, from the same posterior, and an arm that crosses
  `efficacy_threshold_hc_prob` stops for efficacy. The interim success value is
  `median_pfs_success_threshold_arms`, or `null_median_arms` when that is
  `NULL`, as in the R body. With the rule `TRUE` and
  `efficacy_threshold_hc_prob` `NULL`, `NA` or infinite, the R body never stops
  for efficacy, and neither do the engines (the result is the rule `FALSE`
  result). Frames where the R body itself errors (an interim success value that
  is `NULL` or `NA` with the rule on, a threshold of length above one) and a
  non-numeric threshold still run the R body.
* **`efficacy_stopping_rule_vs_ref = TRUE`** (between-arm, PH model).
  `Pr(log HR < log(1 - compare_arms_hr_eff_margin))` is checked before the
  futility probability at each interim look.
* **`compare_arms_hr_margin` and `compare_arms_hr_eff_margin`** may each be
  `NULL` (margin 0) or one finite number, so the `current`, `symmetric` and
  `no_margin` rule variants all run on the fast engines. A negative margin is
  clamped to 0 and one above 0.95 is capped at 0.95, as in the R body. `NA`, a
  vector of length above one and a non-numeric margin still run the R body.

An arm stopped for efficacy at an interim is treated downstream exactly as one
stopped for futility (no final analysis, its enrolment at the stop as its final
sample size). Power is early plus final efficacy, and `Var_PET_Efficacy` is
filled under `return_variance = TRUE`.

**Results change by design on these frames under the default engine.** The
bidirectional rows of the paper wrappers (interim efficacy stopping on) and the
symmetric rule variant used to run the R body and now run the exact engine.
`options(evolveTrial.engine = "cpp")` gives the draw-based engine, which is
bit-identical to the R body on every frame it admits, and
`options(evolveTrial.use_cpp_engine = FALSE)` runs the R body itself. Frames
admitted in 0.0.0.9004 give the same results as before, apart from the last
bit of `Pr_Reach_Max_N` (below).

How large the change is. The R body estimates the interim efficacy probability
from `num_posterior_draws` draws, and that noise at repeated looks adds false
efficacy stops. On the between-arm R23 design (efficacy threshold 0.95, 800
draws) the R body's type I error was about 0.004 above the exact engine's and
its expected sample size under the null about 0.2 patients (0.7 to 1%) lower.
The gap shrinks as the draw count grows (at 12,800 draws it is within Monte
Carlo error of 0), so the exact engine gives the value the R body tends to. The
single-arm bidirectional designs tested matched within the phase 2 margins.
Calibrations with a tight type I error constraint can move when rerun.

## Summary fix

`Pr_Reach_Max_N` is computed from the counts, `(n - early efficacy stops -
early futility stops) / n`, instead of `1 - PET_Efficacy - PET_Futility`, and
the `Var_` proportion variances clamp p to [0, 1]. The old form gave -1.1e-16
at 7 + 93 early stops in 100 trials, so `Var_Pr_Reach_Max_N` was negative.
Every engine shares this summary. Values move by at most about 1e-16.

# evolveTrial 0.0.0.9004

## Fast engines cover the default single-arm call (phase 3, step P3-1)

Both C++ engines (the draw-based `"cpp"` engine and the default exact engine)
now run six configurations that used to fall back to the R body. A call to
`run_simulation_pure()` with these options left at their defaults now runs on
the exact engine.

* **`looks_stop_at_final = FALSE` (the amd8 chronology).** Accrual closes at a
  finite `max_trial_time`, the look at the end of enrolment and the
  post-enrolment looks run only strictly before `T_final`, the last-slot and
  vs-reference phantom-arrival rules apply, and the chronology invariants stop
  the run with the R body's messages (for example with
  `min_follow_up_at_final < 0`).
* **`random_dropout = TRUE`** (the default). The uniform dropout time is kept as
  the censoring time. In the exact engine it is one extra draw per patient from
  the trial's own generator.
* **`max_trial_time = Inf`** (the default) under both chronologies. It needs
  every arm's `max_total_patients_per_arm` above 0 (otherwise the person-time
  gate denominator is `0 * Inf` and the R body errors, so such frames still run
  the R body).
* **`rebalance_after_events = NULL`** (the default): no interval rebalancing.
* **Median-threshold fallbacks.** `median_pfs_success_threshold_arms` and
  `median_pfs_futility_threshold_arms` may be `NULL` or lack the arm, giving the
  R body's `Inf` / `-Inf` fallbacks (final probability 0). At the interim the
  thresholds fall back to `null_median_arms` and `futility_median_arms` as
  before. An interim threshold vector that is not `NULL` must name the arm,
  because the R body errors otherwise.
* **`futility_threshold_hc_prob = NULL`** (the default): no interim futility
  stop.

Everything else in the supported set is unchanged (multi-arm vs historical
control, interim efficacy stopping, person-time milestones, unnamed caps, non-PH
between-arm models and the other excluded items still run the R body).

**Results change by design on these frames under the default engine.** Calls
that used to fall back now run the exact engine, whose posterior probabilities
are exact rather than Monte Carlo and whose patient-level draws come from
per-trial streams. This includes the paper wrappers on futility-only designs,
which leave `looks_stop_at_final` unset. `options(evolveTrial.engine = "cpp")`
gives the draw-based engine, which is bit-identical to the R body on every
frame it admits, and `options(evolveTrial.use_cpp_engine = FALSE)` runs the R
body itself. Frames the supported set already admitted in 0.0.0.9003 (for
example the atbo certify shims, `looks_stop_at_final = TRUE` with dropout off)
give the same results as before. A `looks_stop_at_final = TRUE` frame with
dropout on, `max_trial_time = Inf`, no rebalancing, a threshold fallback or no
interim futility rule used to run the R body and now runs the exact engine.

# evolveTrial 0.0.0.9003

## Merge of main (amd8 chronology) into the engine branch

Both lines of development called themselves 0.0.0.9002. This release merges
them. The two 0.0.0.9002 sections below are kept as written and labeled by
branch.

* **amd8 chronology on the `looks_stop_at_final = FALSE` path.** With
  `looks_stop_at_final = FALSE` (the default, also the default of
  `options(evolveTrial.looks_stop_at_final)`) the R body follows main's
  0.0.0.9002 chronology (BATON-C freeze Amendment 8): accrual closes at
  `max_trial_time`, interims run only strictly before
  `T_final = min(last enrolment + min_follow_up_at_final, max_trial_time)`,
  and the per-trial chronology invariants are asserted. This replaces the
  422eb40 schedule that `FALSE` reproduced on the engine branch, so results on
  the default path change exactly as main's 0.0.0.9002 section describes.
* **`looks_stop_at_final = TRUE` is unchanged.** The P6 chronology (ADAPT
  Amendment 2, the schedule both C++ engines implement) is textually the code of
  the engine branch at 37bb518. The amd8 rules and invariants do not apply on
  it, and results on these frames do not change.
* **Exact engine stays the default.** The supported set still admits only
  `looks_stop_at_final = TRUE` frames, so no amd8 frame reaches either C++
  engine. Default calls run the R body with the amd8 chronology.
* `random_dropout` keeps its default `TRUE`, the uniform dropout censoring that
  main's chronology tests assume.
* **Exact engine edge cases (Codex review, two P2 items).** The supported set
  now requires every prior shape and rate (`prior_alpha_params_model`,
  `prior_beta_params_model`) to be finite and `> 0`; other frames run the R body
  under the default and stop under an explicit `engine = "exact"`. A posterior
  scale `w / rate` that overflows (posterior rate near `1e-310`) is handled in
  log space and gives the correct tail, and an exact probability that is not
  finite now stops the run with an error instead of reading as inconclusive.
  `random_dropout` is now documented in `?run_simulation_pure`.

# evolveTrial 0.0.0.9002 (engine branch: exact engine as default)

## Changes to defaults

* **The exact chunk engine is now the default.** With `options(evolveTrial.engine)`
  unset, `run_simulation_pure()` and `run_scenarios()` compute posterior
  stopping probabilities in closed form (the phase 2 engine,
  `src/simulate_chunk_exact.cpp`) instead of from Monte Carlo posterior draws.
  This is the infinite-draw limit of the previous rule, so operating
  characteristics agree within Monte Carlo error; finite-draw runs stop slightly
  more often for futility. Results at a given seed differ from 0.0.0.9001, and
  they no longer depend on `num_workers` or the block size. PI ruling 2026-11-04.
  The supported set requires `options(evolveTrial.looks_stop_at_final = TRUE)`;
  that option still defaults to `FALSE`, so a call with every option at its
  default still runs the R body, unchanged.
* Configurations outside the supported set of the C++ engines (hybrid,
  seamless, binary and unsupported single-arm or between-arm frames) fall back
  silently to the draw-based engine under the default, as before. The one-time
  "exact is set but unused" notice now fires only when the option was set
  explicitly, and only an explicit `options(evolveTrial.engine = "exact")` turns
  an unsupported configuration into an error.
* `options(evolveTrial.engine = "cpp")` restores the draw-based C++ engine,
  bit-identical to the R body (`options(evolveTrial.use_cpp_engine = FALSE)`).

## Tests

* New `tests/testthat/test-engine-equivalence.R`: draw-based engine against the
  R body, exact tail probabilities against the original R posterior samplers,
  exact against R-body operating characteristics within Monte Carlo
  tolerances, golden-value snapshots of both engines, and default selection.
  Larger versions run when `EVOLVETRIAL_LONG_TESTS` is set.

# evolveTrial 0.0.0.9002 (main: amd8 chronology, BATON-C freeze Amendment 8)

## Chronology change (Amendment 8, item 1): calendar cutoff

`run_simulation_pure()` now enforces `max_trial_time` as a hard calendar
cutoff and fixes the final-analysis time. Results from 0.0.0.9001 and earlier
are not comparable for trials whose chronology reaches `max_trial_time` or
that ran interim looks at or after the final-analysis time.

* **Accrual closes at `max_trial_time`.** An arrival at or after the cutoff
  enrols nobody and triggers no interim; enrolment ends there. Previously the
  enrolment loop had no time condition, patients were enrolled after the
  cutoff, counted in N and dropped from the analysis, and interim looks ran
  after the cutoff on post-cutoff data.
* **Final-analysis time.** For a trial not stopped early the final analysis
  is at `T_final = min(last enrolment + min_follow_up_at_final,
  max_trial_time)` (unchanged formula; the last enrolment is now always
  before the cutoff).
* **Interims only before `T_final`.** Every interim look (calendar,
  person-time, backstop, the look at the moment enrolment ends, post-enrolment
  looks) runs only strictly before `T_final`. Previously post-enrolment looks
  ran until `max_trial_time` even when `T_final` was earlier, so a trial could
  stop after its own final analysis. With `min_follow_up_at_final = 0` (the
  default) the look at the moment enrolment ends coincides with `T_final` and
  no longer runs; this also changes results for `max_trial_time = Inf`
  callers that use the default.
* **Invariants asserted.** Each simulated trial stops with an error if any
  enrolment, stop or final-analysis time exceeds a finite `max_trial_time`,
  or any enrolment lies after the final-analysis time.
* RNG stream: an arrival at or after the cutoff consumes only the arrival
  draw. Trials whose chronology never reaches the cutoff and that ran no look
  at or after `T_final` reproduce 0.0.0.9001 draw for draw.
* Known limitation: with `min_follow_up_at_final = 0` and a finite cutoff, the
  look at the last arrival before the cutoff can coincide with `T_final`
  (whether that arrival is the last is only known when the next arrival lands
  at or after the cutoff); the same holds for some multi-arm partial-stop
  exits. With `min_follow_up_at_final > 0` in single-arm mode no look can
  reach `T_final`.
* The C++ hybrid path (`src/hybrid_trial_sim.cpp`) is unchanged.

# evolveTrial 0.0.0.9001

## Bug Fixes

Hazard-ratio (HR) margin correction for the between-arm (BA) comparison. The
benefit-side futility margin is now applied consistently across the R and C++
decision and forecasting paths.

* **F1 — Between-arm futility margin orientation.** The between-arm futility
  comparator now uses the benefit-side margin `Pr(HR >= 1 - delta_HR)`. The
  previous code used the harm-side `1 + delta_HR`, which buffered the futility
  decision in the wrong direction.

* **D2 — `rule_variant = "symmetric"` for between-arm/multi-arm.** The
  `"symmetric"` rule variant now buffers the between-arm efficacy comparator
  (efficacy requires `Pr(HR < 1 - delta_HR_eff)`). Previously this was a silent
  no-op for the BA/multi-arm path, so `"symmetric"` behaved identically to the
  unbuffered comparator there.

* **D3 — C++ hybrid/seamless simulator HR margin.** The C++ hybrid/seamless
  simulator now applies the benefit-side `delta_HR` margin in both the Stage-2
  between-arm decision and the SA-to-BA conversion predictive-probability
  forecast, scored through the proportional-hazards `logHR` model so the
  forecast uses the same rule as the realized decision.

* **D4 — Hybrid Stage-1 single-arm comparator.** The hybrid Stage-1 single-arm
  comparator is fixed to the control median `M_ctrl` (`hr_threshold_sa = 1.0`).
  It was previously `0.80`, which silently imposed an unintended HR margin on
  the single-arm phase.

## Code Quality

* Removed unused C++ helpers from `src/hybrid_trial_sim.cpp`:
  `compute_ba_posterior_internal` (and its forward declaration) and
  `compute_ba_posterior_internal_aggregated`. Both had zero call sites after
  the BA decision and conversion-PP paths were migrated to
  `compute_ba_ph_posterior_internal` (the PH-model benefit-side posterior).

# evolveTrial 0.0.0.9000

## Breaking Changes

* **Single-arm gate scaling**: Single-arm (`compare_arms_option = FALSE`) now applies proportional scaling to `min_events_per_arm` and `min_person_time_frac_per_arm` based on `randomization_probs`, matching vs-reference behavior. Arms with lower randomization probabilities automatically get proportionally lower gates when scalar gate values are provided.

* **Parameter naming harmonization**: Updated parameter names for consistency across comparison paths:
  - `efficacy_threshold_current_prob_hc` → `efficacy_threshold_hc_prob`
  - `posterior_futility_threshold_hc` → `futility_threshold_hc_prob`
  - `min_events_for_analysis` → `min_events_hc`
  - `min_median_followup` → `min_median_followup_hc`

## Deprecated Parameters

Old parameter names are still supported with deprecation warnings. Users should migrate to new names:

| Old Parameter Name | New Parameter Name | Context |
|-------------------|-------------------|---------|
| `efficacy_threshold_current_prob_hc` | `efficacy_threshold_hc_prob` | Single-arm efficacy threshold |
| `posterior_futility_threshold_hc` | `futility_threshold_hc_prob` | Single-arm futility threshold |
| `min_events_for_analysis` | `min_events_hc` | Single-arm minimum events gate |
| `min_median_followup` | `min_median_followup_hc` | Single-arm minimum median followup gate |

## New Features

* Added `resolve_gate_vec()` internal helper (in `R/gate_diagnostics.R`) for unified gate parameter resolution with optional proportional scaling across both comparison paths.

* Added comprehensive test suites for single-arm scenarios:
  - `test-single-arm-gates.R`: Basic edge cases with proportional scaling, rebalancing, and multiple gates
  - `test-single-arm-regression.R`: Null/alternative scenarios, early stopping, and final analysis validation
  - `test-path-parity.R`: Parity tests ensuring single-arm and vs-reference paths maintain consistent behavior

## Code Quality

* Removed unused functions: `run_single_arm_interim()`, `slice_arm_at_time()`, and `posterior_scalar_draws()` from `R/state_management.R` (no longer needed after refactoring).

* Consolidated gate resolution logic: Both vs-reference and single-arm paths now use the same `resolve_gate_vec()` helper, eliminating code duplication.

## Package Maintenance

* Fixed R CMD check warnings and notes for CRAN compliance:
  - Removed non-ASCII characters (emojis, Greek letters, smart quotes) from source code
  - Replaced `library()` calls with `requireNamespace()` + `::` for proper package dependencies
  - Reorganized DESCRIPTION: moved optional packages (gt, dplyr, openxlsx, etc.) to Suggests
  - Removed unused imports (ggpubr, survival, survminer)
  - Added VignetteBuilder field for proper vignette support

* Documentation improvements:
  - Added complete `@param` documentation for internal helper functions
  - Fixed documentation mismatches in `explore_early_stopping_from_cal()`
  - Added `@keywords internal` tags for non-exported functions

* Code quality enhancements:
  - Added `stats::rnorm` to namespace imports
  - Added missing global variable declarations
  - Updated .Rbuildignore to exclude development files
  - Improved test coverage with actual utility function tests

## Bug Fixes

* None in this development version

## New Features

* None in this development version (package in active development)
