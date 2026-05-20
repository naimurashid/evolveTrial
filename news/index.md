# Changelog

## evolveTrial 0.0.0.9001

### Bug Fixes

Hazard-ratio (HR) margin correction for the between-arm (BA) comparison.
The benefit-side futility margin is now applied consistently across the
R and C++ decision and forecasting paths.

- **F1 — Between-arm futility margin orientation.** The between-arm
  futility comparator now uses the benefit-side margin
  `Pr(HR >= 1 - delta_HR)`. The previous code used the harm-side
  `1 + delta_HR`, which buffered the futility decision in the wrong
  direction.

- **D2 — `rule_variant = "symmetric"` for between-arm/multi-arm.** The
  `"symmetric"` rule variant now buffers the between-arm efficacy
  comparator (efficacy requires `Pr(HR < 1 - delta_HR_eff)`). Previously
  this was a silent no-op for the BA/multi-arm path, so `"symmetric"`
  behaved identically to the unbuffered comparator there.

- **D3 — C++ hybrid/seamless simulator HR margin.** The C++
  hybrid/seamless simulator now applies the benefit-side `delta_HR`
  margin in both the Stage-2 between-arm decision and the SA-to-BA
  conversion predictive-probability forecast, scored through the
  proportional-hazards `logHR` model so the forecast uses the same rule
  as the realized decision.

- **D4 — Hybrid Stage-1 single-arm comparator.** The hybrid Stage-1
  single-arm comparator is fixed to the control median `M_ctrl`
  (`hr_threshold_sa = 1.0`). It was previously `0.80`, which silently
  imposed an unintended HR margin on the single-arm phase.

### Code Quality

- Removed unused C++ helpers from `src/hybrid_trial_sim.cpp`:
  `compute_ba_posterior_internal` (and its forward declaration) and
  `compute_ba_posterior_internal_aggregated`. Both had zero call sites
  after the BA decision and conversion-PP paths were migrated to
  `compute_ba_ph_posterior_internal` (the PH-model benefit-side
  posterior).

## evolveTrial 0.0.0.9000

### Breaking Changes

- **Single-arm gate scaling**: Single-arm
  (`compare_arms_option = FALSE`) now applies proportional scaling to
  `min_events_per_arm` and `min_person_time_frac_per_arm` based on
  `randomization_probs`, matching vs-reference behavior. Arms with lower
  randomization probabilities automatically get proportionally lower
  gates when scalar gate values are provided.

- **Parameter naming harmonization**: Updated parameter names for
  consistency across comparison paths:

  - `efficacy_threshold_current_prob_hc` → `efficacy_threshold_hc_prob`
  - `posterior_futility_threshold_hc` → `futility_threshold_hc_prob`
  - `min_events_for_analysis` → `min_events_hc`
  - `min_median_followup` → `min_median_followup_hc`

### Deprecated Parameters

Old parameter names are still supported with deprecation warnings. Users
should migrate to new names:

| Old Parameter Name | New Parameter Name | Context |
|----|----|----|
| `efficacy_threshold_current_prob_hc` | `efficacy_threshold_hc_prob` | Single-arm efficacy threshold |
| `posterior_futility_threshold_hc` | `futility_threshold_hc_prob` | Single-arm futility threshold |
| `min_events_for_analysis` | `min_events_hc` | Single-arm minimum events gate |
| `min_median_followup` | `min_median_followup_hc` | Single-arm minimum median followup gate |

### New Features

- Added [`resolve_gate_vec()`](../reference/resolve_gate_vec.md)
  internal helper (in `R/gate_diagnostics.R`) for unified gate parameter
  resolution with optional proportional scaling across both comparison
  paths.

- Added comprehensive test suites for single-arm scenarios:

  - `test-single-arm-gates.R`: Basic edge cases with proportional
    scaling, rebalancing, and multiple gates
  - `test-single-arm-regression.R`: Null/alternative scenarios, early
    stopping, and final analysis validation
  - `test-path-parity.R`: Parity tests ensuring single-arm and
    vs-reference paths maintain consistent behavior

### Code Quality

- Removed unused functions: `run_single_arm_interim()`,
  `slice_arm_at_time()`, and `posterior_scalar_draws()` from
  `R/state_management.R` (no longer needed after refactoring).

- Consolidated gate resolution logic: Both vs-reference and single-arm
  paths now use the same
  [`resolve_gate_vec()`](../reference/resolve_gate_vec.md) helper,
  eliminating code duplication.

### Package Maintenance

- Fixed R CMD check warnings and notes for CRAN compliance:
  - Removed non-ASCII characters (emojis, Greek letters, smart quotes)
    from source code
  - Replaced [`library()`](https://rdrr.io/r/base/library.html) calls
    with [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) +
    `::` for proper package dependencies
  - Reorganized DESCRIPTION: moved optional packages (gt, dplyr,
    openxlsx, etc.) to Suggests
  - Removed unused imports (ggpubr, survival, survminer)
  - Added VignetteBuilder field for proper vignette support
- Documentation improvements:
  - Added complete `@param` documentation for internal helper functions
  - Fixed documentation mismatches in
    [`explore_early_stopping_from_cal()`](../reference/explore_early_stopping_from_cal.md)
  - Added `@keywords internal` tags for non-exported functions
- Code quality enhancements:
  - Added [`stats::rnorm`](https://rdrr.io/r/stats/Normal.html) to
    namespace imports
  - Added missing global variable declarations
  - Updated .Rbuildignore to exclude development files
  - Improved test coverage with actual utility function tests

### Bug Fixes

- None in this development version

### New Features

- None in this development version (package in active development)
