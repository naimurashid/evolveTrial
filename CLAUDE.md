# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**evolveTrial** is an R package for Bayesian adaptive platform/umbrella clinical trial simulation and design analysis. It specializes in progression-free survival (PFS) endpoints using piecewise exponential models with both single-arm (vs historical control) and multi-arm (vs-reference) decision logic.

## Package Structure

### Core Simulation Architecture

The simulation engine is built around three main components:

1. **State Management** (`R/state_management.R`): The `make_state()` function creates trial state containers tracking arm status (recruiting/stopped), patient registries, and enrollment counts. The `slice_arm_data_at_time()` function provides time-based snapshots of patient data for interim analyses, computing interval metrics via `calculate_interval_metrics_fast()`.

2. **Interim Logic** (`R/interim_logic.R`): The `interim_check()` function is the decision engine that evaluates efficacy and futility at scheduled looks. It branches on `compare_arms_option`:
   - **vs-reference path**: Compares experimental arms against a reference/control arm using `calculate_current_probs_vs_ref()`, evaluating per-arm gates (`min_events_per_arm`, `min_median_followup_per_arm`, `min_person_time_frac_per_arm`)
   - **single-arm path**: Evaluates arms independently against historical control targets using `calculate_current_probs_hc()`, with optional predictive probability fallbacks from `R/predictive_probabilities.R`

3. **Simulation Driver** (`R/simulation_driver.R`): The `run_simulation_pure()` function is the workhorse that orchestrates Monte Carlo replicates, enrolling patients via piecewise exponential data generation (`R/data_generation.R`), applying interim checks at calendar beats or person-time milestones, and accumulating operating characteristics (type I error, power, PET, expected N).

### Statistical Helpers

- **Posterior Sampling** (`R/posterior_helpers.R`): `draw_posterior_hazard_samples()` draws from Gamma posteriors for interval hazards; `calculate_median_survival_piecewise()` converts hazard vectors to median survival times (handles open-ended tail extrapolation)
- **Predictive Probabilities** (`R/predictive_probabilities.R`): `calculate_predicted_success_prob_vs_hc()` performs forward simulation from current posterior to final analysis
- **PH Models**: When `use_ph_model_vs_ref = TRUE`, vs-reference comparisons use a joint proportional-hazards model with log hazard ratio priors

### Design Analysis Workflow

1. **Parameter Grid Setup**: Define design grids (efficacy/futility thresholds, gates, sample sizes)
2. **Simulation**: `run_scenarios()` wraps `run_simulation_pure()` for batch execution across scenarios (see `R/design_analysis.R`)
3. **Calibration**: `grid_calibrate()` searches efficacy/futility threshold combinations to meet target type I error constraints; `calibrate_alpha()` performs single-dimension searches
4. **Reporting**: `pretty_scenario_matrix()` pivots results to scenario × arm tables; `export_scenario_table_to_excel()` and `export_scenario_table_to_png()` generate formatted outputs

### Gate Diagnostics

When early stopping fails to occur, it often means information gates postpone the first informative look. Use `estimate_vsref_gate_timing()` (from `R/gate_diagnostics.R`) to get lower-bound heuristics for when vs-reference gates (`min_events_per_arm`, `min_median_followup_per_arm`, `min_person_time_frac_per_arm`) can all be satisfied under deterministic accrual. See README.Rmd for examples.

## Development Commands

### Package Development
```r
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
```bash
# Isolated environment for scenario grids or one-off analyses
Rscript path/to/script.R
```

### Testing
```r
# Run a single test file
testthat::test_file("tests/testthat/test-estimate-vsref-gates.R")

# Run tests matching a pattern
testthat::test_local(filter = "utils")
```

## Coding Conventions

- **Style**: Follow tidyverse conventions (two-space indentation, `{` on same line, `snake_case` identifiers like `min_person_time_frac_per_arm`)
- **Utilities**: Reuse package helpers (`%||%` for fallback, `coalesce_num()` for numeric coalescing) instead of reimplementing
- **Diagnostics**: Keep `message()` calls concise and prefixed (e.g., `[vsREF gate]`)
- **Argument Threading**: When adding parameters, add them to `run_simulation_pure()` formals so downstream helpers can access them via `modifyList()`
- **Documentation**: Use roxygen2 markdown format; document all exported functions and user-facing helpers

## Simulation Configuration

Key design parameters threaded through `run_simulation_pure()` and `interim_check()`:

- **Arms & Comparison**: `arm_names`, `reference_arm_name`, `compare_arms_option` (TRUE = vs-reference, FALSE = single-arm)
- **Truth & Priors**: `weibull_shape_true_arms`, `weibull_median_true_arms`, `prior_alpha_params_model`, `prior_beta_params_model`
- **Thresholds**:
  - vs-reference: `efficacy_threshold_vs_ref_prob`, `futility_threshold_vs_ref_prob`, `compare_arms_futility_margin`
  - single-arm: `efficacy_threshold_current_prob_hc`, `posterior_futility_threshold_hc`, `null_median_arms`, `futility_median_arms`
- **Gates**:
  - vs-reference per-arm: `min_events_per_arm`, `min_median_followup_per_arm`, `min_person_time_frac_per_arm`
  - single-arm global: `min_events_for_analysis`, `min_median_followup`, `min_patients_for_analysis`
- **Scheduling**: `interim_calendar_beat` (fixed spacing) or `person_time_milestones` (data-driven)
- **Sample Size**: `max_total_patients_per_arm`, `cohort_size_per_arm`, `overall_accrual_rate`, `randomization_probs`

## Testing Approach

New decision logic or gating features require regression tests in `tests/testthat/`. Prioritize:
- Multi-arm comparisons with reference arm switches
- Proportional gate scaling across unbalanced arms
- Predictive probability fallback triggers
- Use deterministic seeds (`set.seed(4242)`) and low simulation counts (e.g., `num_simulations = 50`) for unit tests
- Document required outputs (PET, alpha, expected N) in `expect_*()` calls to catch behavior shifts

## Project Management

- **Dependencies**: Managed via `renv/` for reproducibility; `renv::restore()` to sync
- **Documentation**: `man/` holds function docs; `vignettes/` contains long-form guides; `_pkgdown.yml` configures the browsable site
- **Artifacts**: Keep ad-hoc scenario grids and analysis scripts in project root or `inst/`; exclude from package namespace
- **Commit Style**: Imperative, present-tense titles under 72 characters (e.g., "Add one-time interval rebalancing support"); group related code/docs/tests into single commits; link to upstream issues for context
