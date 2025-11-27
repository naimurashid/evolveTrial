# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Overview

**evolveTrial** is an R package for Bayesian adaptive platform/umbrella
clinical trial simulation and design analysis. It specializes in
progression-free survival (PFS) endpoints using piecewise exponential
models with both single-arm (vs historical control) and multi-arm
(vs-reference) decision logic.

## Package Structure

### Core Simulation Architecture

The simulation engine is built around three main components:

1.  **State Management** (`R/state_management.R`): The `make_state()`
    function creates trial state containers tracking arm status
    (recruiting/stopped), patient registries, and enrollment counts. The
    `slice_arm_data_at_time()` function provides time-based snapshots of
    patient data for interim analyses, computing interval metrics via
    [`calculate_interval_metrics_fast()`](reference/calculate_interval_metrics_fast.md).

2.  **Interim Logic** (`R/interim_logic.R`): The `interim_check()`
    function is the decision engine that evaluates efficacy and futility
    at scheduled looks. It branches on `compare_arms_option`:

    - **vs-reference path**: Compares experimental arms against a
      reference/control arm using `calculate_current_probs_vs_ref()`,
      evaluating per-arm gates (`min_events_per_arm`,
      `min_median_followup_per_arm`, `min_person_time_frac_per_arm`)
    - **single-arm path**: Evaluates arms independently against
      historical control targets using `calculate_current_probs_hc()`,
      with optional predictive probability fallbacks from
      `R/predictive_probabilities.R`

3.  **Simulation Driver** (`R/simulation_driver.R`): The
    [`run_simulation_pure()`](reference/run_simulation_pure.md) function
    is the workhorse that orchestrates Monte Carlo replicates, enrolling
    patients via piecewise exponential data generation
    (`R/data_generation.R`), applying interim checks at calendar beats
    or person-time milestones, and accumulating operating
    characteristics (type I error, power, PET, expected N).

### Statistical Helpers

- **Posterior Sampling** (`R/posterior_helpers.R`):
  [`draw_posterior_hazard_samples()`](reference/draw_posterior_hazard_samples.md)
  draws from Gamma posteriors for interval hazards;
  [`calculate_median_survival_piecewise()`](reference/calculate_median_survival_piecewise.md)
  converts hazard vectors to median survival times (handles open-ended
  tail extrapolation)
- **Predictive Probabilities** (`R/predictive_probabilities.R`):
  `calculate_predicted_success_prob_vs_hc()` performs forward simulation
  from current posterior to final analysis
- **PH Models**: When `use_ph_model_vs_ref = TRUE`, vs-reference
  comparisons use a joint proportional-hazards model with log hazard
  ratio priors

### Design Analysis Workflow

1.  **Parameter Grid Setup**: Define design grids (efficacy/futility
    thresholds, gates, sample sizes)
2.  **Simulation**: [`run_scenarios()`](reference/run_scenarios.md)
    wraps [`run_simulation_pure()`](reference/run_simulation_pure.md)
    for batch execution across scenarios (see `R/design_analysis.R`)
3.  **Calibration**: [`grid_calibrate()`](reference/grid_calibrate.md)
    searches efficacy/futility threshold combinations to meet target
    type I error constraints;
    [`calibrate_alpha()`](reference/calibrate_alpha.md) performs
    single-dimension searches
4.  **Reporting**:
    [`pretty_scenario_matrix()`](reference/pretty_scenario_matrix.md)
    pivots results to scenario × arm tables;
    [`export_scenario_table_to_excel()`](reference/export_scenario_table_to_excel.md)
    and
    [`export_scenario_table_to_png()`](reference/export_scenario_table_to_png.md)
    generate formatted outputs

### Gate Diagnostics

When early stopping fails to occur, it often means information gates
postpone the first informative look. Use
[`estimate_vsref_gate_timing()`](reference/estimate_vsref_gate_timing.md)
(from `R/gate_diagnostics.R`) to get lower-bound heuristics for when
vs-reference gates (`min_events_per_arm`, `min_median_followup_per_arm`,
`min_person_time_frac_per_arm`) can all be satisfied under deterministic
accrual. See README.Rmd for examples.

## Development Commands

### Package Development

``` r
# Iterative development (refreshes namespace without reinstall)
devtools::load_all()

# Run tests
devtools::test()

# Full R CMD check before publishing
devtools::check()

# Install locally to mirror user experience
devtools::install()
library(evolveTrial)

# Regenerate documentation site after API changes
pkgdown::build_site()

# Render README.Rmd to README.md
devtools::build_readme()
```

### Running Scripts

``` bash
# Isolated environment for scenario grids or one-off analyses
Rscript path/to/script.R
```

### Testing

``` r
# Run a single test file
testthat::test_file("tests/testthat/test-estimate-vsref-gates.R")

# Run tests matching a pattern
testthat::test_local(filter = "utils")
```

## Coding Conventions

- **Style**: Follow tidyverse conventions (two-space indentation, `{` on
  same line, `snake_case` identifiers like
  `min_person_time_frac_per_arm`)
- **Utilities**: Reuse package helpers (`%||%` for fallback,
  `coalesce_num()` for numeric coalescing,
  [`resolve_gate_vec()`](reference/resolve_gate_vec.md) for gate
  parameter resolution) instead of reimplementing
- **Diagnostics**: Keep
  [`message()`](https://rdrr.io/r/base/message.html) calls concise and
  prefixed (e.g., `[vsREF gate]`)
- **Argument Threading**: When adding parameters, add them to
  [`run_simulation_pure()`](reference/run_simulation_pure.md) formals so
  downstream helpers can access them via
  [`modifyList()`](https://rdrr.io/r/utils/modifyList.html)
- **Documentation**: Use roxygen2 markdown format; document all exported
  functions and user-facing helpers
- **Parameter Naming**: Use harmonized naming scheme across comparison
  paths (e.g., `efficacy_threshold_hc_prob`,
  `futility_threshold_hc_prob`)

## Simulation Configuration

Key design parameters threaded through
[`run_simulation_pure()`](reference/run_simulation_pure.md) and
`interim_check()`:

- **Arms & Comparison**: `arm_names`, `reference_arm_name`,
  `compare_arms_option` (TRUE = vs-reference, FALSE = single-arm)
- **Truth & Priors**: `weibull_shape_true_arms`,
  `weibull_median_true_arms`, `prior_alpha_params_model`,
  `prior_beta_params_model`
- **Thresholds** (harmonized naming):
  - vs-reference: `efficacy_threshold_vs_ref_prob`,
    `futility_threshold_vs_ref_prob`, `compare_arms_futility_margin`
  - single-arm: `efficacy_threshold_hc_prob`,
    `futility_threshold_hc_prob`, `null_median_arms`,
    `futility_median_arms`
  - **Deprecated (still supported with warnings)**:
    `efficacy_threshold_current_prob_hc`,
    `posterior_futility_threshold_hc`
- **Gates** (now use proportional scaling for both paths):
  - Per-arm gates: `min_events_per_arm`, `min_median_followup_per_arm`,
    `min_person_time_frac_per_arm` (used by both vs-reference and
    single-arm)
  - Legacy global gates: `min_events_hc`, `min_median_followup_hc`,
    `min_patients_for_analysis` (single-arm only)
  - **Deprecated (still supported with warnings)**:
    `min_events_for_analysis`, `min_median_followup`
  - **Proportional Scaling**: When `min_events_per_arm` or
    `min_person_time_frac_per_arm` are scalar, they are automatically
    scaled by `randomization_probs` to account for unbalanced
    randomization
- **Scheduling**: `interim_calendar_beat` (fixed spacing) or
  `person_time_milestones` (data-driven)
- **Sample Size**: `max_total_patients_per_arm`, `cohort_size_per_arm`,
  `overall_accrual_rate`, `randomization_probs`

## Gate Scaling Behavior

As of recent updates, **both vs-reference and single-arm paths apply
proportional gate scaling** when scalar gate values are provided:

- **What gets scaled**: `min_events_per_arm` and
  `min_person_time_frac_per_arm` (when provided as scalars)
- **What doesn’t get scaled**: `min_median_followup_per_arm` (followup
  time is independent of randomization ratio)
- **How scaling works**: If you specify `min_events_per_arm = 10` with
  `randomization_probs = c(Control = 0.3, Treatment = 0.7)`, the helper
  [`resolve_gate_vec()`](reference/resolve_gate_vec.md) (in
  `R/gate_diagnostics.R`) will scale gates proportionally: arms with
  lower randomization probabilities get proportionally lower gates
- **When to use named vectors**: To override proportional scaling,
  provide a named vector (e.g.,
  `min_events_per_arm = c(Control = 8, Treatment = 12)`)

This ensures fairness when comparing arms with unbalanced randomization,
and unifies behavior across comparison paths.

## Testing Approach

New decision logic or gating features require regression tests in
`tests/testthat/`. Prioritize: - Multi-arm comparisons with reference
arm switches - Proportional gate scaling across unbalanced arms (see
`test-single-arm-gates.R` and `test-path-parity.R`) - Predictive
probability fallback triggers - Use deterministic seeds
(`set.seed(4242)`) and low simulation counts (e.g.,
`num_simulations = 50`) for unit tests - Document required outputs (PET,
alpha, expected N) in `expect_*()` calls to catch behavior shifts -
Parity tests ensure single-arm and vs-reference paths maintain
consistent behavior

## Project Management

- **Dependencies**: Managed via `renv/` for reproducibility;
  `renv::restore()` to sync
- **Documentation**: `man/` holds function docs; `vignettes/` contains
  long-form guides; `_pkgdown.yml` configures the browsable site
- **Artifacts**: Keep ad-hoc scenario grids and analysis scripts in
  project root or `inst/`; exclude from package namespace
- **Commit Style**: Imperative, present-tense titles under 72 characters
  (e.g., “Add one-time interval rebalancing support”); group related
  code/docs/tests into single commits; link to upstream issues for
  context
