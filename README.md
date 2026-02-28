# evolveTrial

<!-- badges: start -->
[![R-CMD-check](https://github.com/naimurashid/evolveTrial/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/naimurashid/evolveTrial/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**evolveTrial** is an R package for designing and simulating Bayesian adaptive
clinical trials with time-to-event or binary endpoints. It supports single-arm,
multi-arm, and hybrid seamless (single-arm-to-between-arm) trial designs with
interim stopping rules based on posterior probabilities.

The package was developed for the ARPA-H ADAPT breast cancer platform trial and
accompanies the JASA paper:

> Rashid, N. (2026). "Constrained Bayesian Optimization for Calibration of
> Bayesian Adaptive Clinical Trials." *Journal of the American Statistical
> Association*.

## Features

- **Single-arm vs historical control**: evaluate efficacy against a historical benchmark with Bayesian posterior probability decision rules
- **Multi-arm adaptive trials**: compare experimental arms to a reference arm under proportional hazards
- **Hybrid seamless designs**: 4-state machine transitioning from single-arm monitoring to between-arm comparison with predictive-probability-based conversion decisions
- **Binary endpoints**: Simon two-stage designs with exact operating characteristics
- **Piecewise exponential hazards**: flexible time-to-event modeling with interval-specific rates
- **C++ acceleration**: posterior sampling, hazard computations, and hybrid trial simulation via Rcpp and RcppArmadillo
- **Design calibration**: grid search over decision thresholds with Pareto-optimal filtering
- **Operating characteristics**: type I error, power, expected sample size, probability of early termination, and conversion rates

## Installation

### C++ Compiler Prerequisites

evolveTrial contains compiled C++ code (via Rcpp/RcppArmadillo) that must be
built during installation. You need a working C++ toolchain.

#### Windows

Install **Rtools** before installing evolveTrial:

1. Download from <https://cran.r-project.org/bin/windows/Rtools/>
2. Choose the version matching your R version (e.g., Rtools44 for R 4.4.x)
3. Run the installer with default settings
4. Restart R/RStudio

Verify:

```r
Sys.which("make")
#> Should return a path like "C:/rtools44/usr/bin/make.exe"
```

#### macOS

Install **Xcode Command Line Tools**:

```bash
xcode-select --install
```

On Apple Silicon (M1/M2/M3/M4), you may also need gfortran for RcppArmadillo:

```bash
brew install gcc
```

Verify:

```r
system("clang++ --version")
```

#### Linux (Ubuntu/Debian)

```bash
sudo apt-get install r-base-dev build-essential liblapack-dev libblas-dev
```

#### Linux (Fedora/RHEL)

```bash
sudo dnf install R-devel gcc-c++ lapack-devel blas-devel
```

### Install from GitHub

R >= 4.1 is recommended due to Rcpp/RcppArmadillo usage.

```r
# Using pak (recommended)
install.packages("pak")
pak::pak("naimurashid/evolveTrial")

# Or using remotes
install.packages("remotes")
remotes::install_github("naimurashid/evolveTrial")

# Or using devtools
install.packages("devtools")
devtools::install_github("naimurashid/evolveTrial")
```

### Install from Local Source

```bash
git clone https://github.com/naimurashid/evolveTrial.git
```

```r
devtools::install("path/to/evolveTrial")
# Or for development:
devtools::load_all("path/to/evolveTrial")
```

### Troubleshooting

| Error | Fix |
|-------|-----|
| `compilation failed` or `make: not found` | Install Rtools (Windows) or Xcode CLI tools (macOS) |
| `'RcppArmadillo.h' file not found` | `install.packages("RcppArmadillo")` then retry |
| `library not found for -lgfortran` (macOS) | `brew install gcc` |
| Installation times out | `remotes::install_github("naimurashid/evolveTrial", build_vignettes = FALSE)` |

## Quick Start Examples

### Example 1: Single-Arm Survival Trial

Simulate a single-arm trial evaluating an experimental therapy against a
historical control median of 6 months, with an expected treatment median of
9 months.

```r
library(evolveTrial)

result <- run_simulation_pure(
  num_simulations           = 500,
  arm_names                 = "Experimental",
  reference_arm_name        = "Experimental",
  compare_arms_option       = "none",
  weibull_shape_true_arms   = c(Experimental = 1),
  weibull_median_true_arms  = c(Experimental = 9),
  interval_cutpoints_sim    = c(0, 3, 6, 9, 12, 15, 18, 21, 24),
  max_follow_up_sim         = 24,
  censor_max_time_sim       = 24,
  prior_alpha_params_model  = rep(0.01, 8),
  prior_beta_params_model   = rep(0.01, 8),
  num_posterior_draws       = 5000,
  cohort_size_per_arm       = c(Experimental = 1),
  max_total_patients_per_arm = c(Experimental = 60),
  efficacy_stopping_rule_hc = TRUE,
  efficacy_threshold_hc_prob = 0.95,
  futility_stopping_rule_hc  = TRUE,
  futility_threshold_hc_prob = 0.05,
  median_pfs_success_threshold_arms = c(Experimental = 6),
  overall_accrual_rate      = 2,
  randomization_probs       = c(Experimental = 1),
  min_events_hc             = 15,
  interim_calendar_beat     = 3
)

# Result is a data.frame -- inspect directly
print(result)
#>      Arm_Name True_Median Type_I_Error_or_Power PET_Efficacy PET_Futility ...
#> 1 Experimental           9                  0.82         0.45         0.02 ...
```

### Example 2: Binary Simon Two-Stage Design

Find an optimal Simon design and verify its operating characteristics via
simulation.

```r
library(evolveTrial)

# Find optimal Simon design: H0: p=0.20, H1: p=0.40
design <- find_simon_design(p0 = 0.20, p1 = 0.40, alpha = 0.05, beta = 0.20)
print(design)
#>   n1 r1  n  r n2 EN_null EN_alt PET_null PET_alt type1_error power criterion
#> 1 13  2 43 11 30    21.8   34.2     0.60    0.26       0.048  0.81   optimal

# Confirm via simulation using the Simon design parameters
result <- run_simulation_binary(
  num_simulations              = 2000,
  arm_names                    = "Arm_A",
  true_response_prob           = c(Arm_A = 0.40),
  n1_per_arm                   = c(Arm_A = design$n1),
  n_total_per_arm              = c(Arm_A = design$n),
  p0                           = 0.20,
  r1_per_arm                   = c(Arm_A = design$r1),
  r_per_arm                    = c(Arm_A = design$r),
  use_simon_rules              = TRUE
)

print(result)
#>   Arm_Name True_Response_Prob Type_I_Error_or_Power PET_Efficacy PET_Futility ...
```

### Example 3: Hybrid Seamless SA-to-BA Trial

Simulate a hybrid trial that begins with single-arm monitoring and can convert
to a between-arm comparison if sufficient evidence accumulates.

```r
library(evolveTrial)

# Define the hybrid design parameter vector
theta <- create_hybrid_theta(
  eff_sa            = 0.90,
  fut_sa            = 0.10,
  hr_threshold_sa   = 0.80,
  ev_sa             = 15,
  nmax_sa           = 40,
  conversion_trigger = "any_single_success",
  pp_go             = 0.70,
  pp_nogo           = 0.20,
  ss_method         = "predictive",
  max_additional_n  = 60,
  eff_ba            = 0.975,
  fut_ba            = 0.05,
  ev_ba             = 15,
  nmax_ba           = 80
)

# Base trial configuration
base_args <- list(
  interval_cutpoints_sim   = c(0, 3, 6, 9, 12, 15, 18, 21, 24),
  prior_alpha_params_model = rep(0.01, 8),
  prior_beta_params_model  = rep(0.01, 8),
  overall_accrual_rate     = 3,
  interim_calendar_beat    = 3,
  num_posterior_draws      = 5000
)

# Scenario: historical median 6 mo, reference 6 mo, experimental 9 mo
scenario <- list(historical_median = 6, ref_median = 6, exp_median = 9)

oc <- compute_hybrid_oc_rcpp(
  hybrid_theta    = theta,
  base_args       = base_args,
  scenario_params = scenario,
  num_simulations = 1000,
  seed            = 42,
  trial_mode      = "hybrid"
)

cat("Power:", round(oc$power, 3), "\n")
cat("Type I:", round(oc$type1, 3), "\n")
cat("E[N] under null:", round(oc$EN_null, 1), "\n")
cat("P(conversion):", round(oc$P_conversion, 3), "\n")
```

## Core API Reference

### Simulation Functions

| Function | Purpose |
|----------|---------|
| `run_simulation_pure()` | Simulate a single-arm or multi-arm trial with time-to-event endpoints |
| `run_scenarios()` | Run multiple scenario configurations over a shared baseline |
| `scenarios_from_grid()` | Generate factorial scenario grid from parameter choices |
| `run_simulation_binary()` | Simulate a binary-endpoint trial (one- or two-stage) |
| `compute_hybrid_oc_rcpp()` | Simulate hybrid seamless SA-to-BA trial via C++ engine |
| `create_hybrid_theta()` | Construct the hybrid design parameter vector |

### Design Calibration

| Function | Purpose |
|----------|---------|
| `grid_calibrate()` | Grid search over thresholds for optimal design calibration |
| `calibrate_alpha()` | Calibrate type I error across null scenarios |
| `explore_early_stopping_from_cal()` | Explore early stopping configurations from a calibration |
| `filter_early_grid()` | Filter designs meeting alpha/power constraints |
| `recommend_design_from_early()` | Select best design from feasible set |
| `adopt_calibration()` | Adopt a calibrated design into trial arguments |

### Binary Endpoint Helpers

| Function | Purpose |
|----------|---------|
| `find_simon_design()` | Find optimal or minimax Simon two-stage design |
| `simon_oc_exact()` | Exact operating characteristics for a Simon design |

### Visualization and Reporting

| Function | Purpose |
|----------|---------|
| `pretty_scenario_matrix()` | Pivot simulation results to wide-format summary |
| `plot_calibration()` | Power vs type I error with Pareto frontier |
| `plot_early_tradeoff()` | Early stopping trade-off visualization |
| `export_scenario_table_to_excel()` | Export results to formatted Excel workbook |
| `export_scenario_table_to_png()` | Export results to PNG via gt |

### Diagnostics

| Function | Purpose |
|----------|---------|
| `estimate_vsref_gate_timing()` | Estimate when information gates can first be satisfied |

### Decision Rule Functions

| Function | Purpose |
|----------|---------|
| `evaluate_sa_efficacy()` | Evaluate single-arm efficacy decision |
| `evaluate_sa_futility()` | Evaluate single-arm futility decision |
| `evaluate_ba_efficacy()` | Evaluate between-arm efficacy decision |
| `evaluate_ba_futility()` | Evaluate between-arm futility decision |
| `evaluate_conversion_trigger()` | Evaluate SA-to-BA conversion readiness |
| `make_conversion_decision()` | Execute conversion decision with SSR |
| `perform_ssr()` | Sample size re-estimation for between-arm phase |

## Package Structure

```
evolveTrial/
├── R/
│   ├── simulation_driver.R          # run_simulation_pure(), run_scenarios(), scenarios_from_grid()
│   ├── binary_endpoint.R            # run_simulation_binary(), find_simon_design(), simon_oc_exact()
│   ├── hybrid_trial.R               # create_hybrid_theta(), hybrid trial orchestration
│   ├── hybrid_sim_rcpp.R            # compute_hybrid_oc_rcpp()
│   ├── hybrid_decisions.R           # SA/BA efficacy/futility evaluators, conversion logic
│   ├── hybrid_ssr.R                 # perform_ssr()
│   ├── state_management.R           # Trial state containers
│   ├── design_analysis.R            # grid_calibrate(), plot_calibration(), pretty_scenario_matrix()
│   ├── gate_diagnostics.R           # estimate_vsref_gate_timing()
│   ├── interim_logic.R              # Stopping rule implementation
│   ├── posterior_helpers.R           # Bayesian posterior computation (R)
│   ├── posterior_cpp_dispatchers.R   # R-to-C++ dispatch for posterior draws
│   ├── predictive_probabilities.R   # Predictive probability computation
│   ├── data_generation.R            # Time-to-event data generation
│   ├── globals.R                    # Package-wide constants
│   ├── evolveTrial-package.R        # Package documentation
│   ├── zzz_imports.R                # Import declarations
│   └── RcppExports.R                # Auto-generated Rcpp bindings
├── src/
│   ├── hybrid_trial_sim.cpp         # C++ hybrid trial simulation engine
│   ├── posterior_sampling.cpp        # C++ posterior sampling with RcppArmadillo
│   └── predictive_probability.cpp   # C++ predictive probability calculations
├── tests/testthat/                  # 16 test files
├── vignettes/
│   └── design-overview.Rmd          # Methodology vignette
├── DESCRIPTION
├── NAMESPACE
└── LICENSE                          # MIT
```

## Companion Package: BATON

[**BATON**](https://github.com/naimurashid/BATON) (Bayesian Adaptive Trial
Optimization) provides constrained Bayesian optimization for calibrating
evolveTrial designs. While evolveTrial evaluates operating characteristics for a
given design configuration, BATON searches over the design space to find
configurations that satisfy frequentist constraints (e.g., type I error control)
while optimizing an objective (e.g., minimizing expected sample size or
maximizing power). Together, they form a complete pipeline for adaptive trial
design calibration.

```r
# Typical workflow:
# 1. Define trial structure in evolveTrial
# 2. Use BATON::bo_calibrate() to search over evolveTrial's parameter space
# 3. Extract optimal design with BATON::summarise_case_study()
```

## Running Tests

```r
devtools::test()
```

## Citation

If you use evolveTrial in your work, please cite:

```bibtex
@article{rashid2026constrained,
  title   = {Constrained {Bayesian} Optimization for Calibration of
             {Bayesian} Adaptive Clinical Trials},
  author  = {Rashid, Naim U.},
  journal = {Journal of the American Statistical Association},
  year    = {2026}
}
```

## License

MIT. See [LICENSE](LICENSE) for details.
