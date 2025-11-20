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
