# evolveTrial

**evolveTrial** provides utilities for designing and simulating Bayesian
adaptive platform and umbrella clinical trials with time-to-event
endpoints. The package supports the ARPA-H ADAPT breast cancer platform
trial design.

## Key Features

- **Multi-arm adaptive trials**: Simulate trials comparing experimental
  arms to a reference arm with proportional hazards
- **Single-arm vs historical control**: Evaluate efficacy relative to
  historical benchmarks
- **Bayesian decision rules**: Interim stopping for efficacy or futility
  based on posterior probabilities
- **Piecewise exponential hazards**: Flexible time-to-event modeling
  with interval-specific hazards
- **Information gates**: Control interim timing via minimum events,
  median follow-up, and person-time requirements
- **Operating characteristics**: Type I error, power, expected sample
  size, and probability of early termination
- **C++ acceleration**: Posterior sampling and hazard computations via
  Rcpp for high performance

## Installation

### Prerequisites for Compiling C++ Code

This package contains C++ code (via Rcpp) that must be compiled during
installation. **Most users will need to install additional tools before
installing evolveTrial.**

#### Windows Users

You must install **Rtools** before installing packages with C++ code: 1.
Download Rtools from: <https://cran.r-project.org/bin/windows/Rtools/>
2. Choose the version matching your R version (e.g., Rtools44 for R
4.4.x) 3. Run the installer with default settings 4. Restart R/RStudio
after installation

To verify Rtools is installed correctly:

``` r
Sys.which("make")
# Should return a path like "C:/rtools44/usr/bin/make.exe"
```

#### macOS Users

You must install **Xcode Command Line Tools**:

``` bash
# Run this in Terminal (not R)
xcode-select --install
```

A dialog will appear - click “Install” and wait for completion (~5-10
minutes).

For **Apple Silicon Macs (M1/M2/M3)**, you may also need gfortran:

``` bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install gfortran
brew install gcc
```

To verify your toolchain is ready:

``` r
# In R, check for C++ compiler
system("clang++ --version")
```

#### Linux Users

Most Linux distributions include the necessary compilers. If not:

``` bash
# Ubuntu/Debian
sudo apt-get install r-base-dev

# Fedora/RHEL
sudo dnf install R-devel
```

### From GitHub (Recommended)

Install the development version from GitHub using one of these methods:

``` r
# Using pak (fastest, handles compilation automatically)
install.packages("pak")
pak::pak("naimurashid/evolveTrial")

# Or using remotes
install.packages("remotes")
remotes::install_github("naimurashid/evolveTrial")

# Or using devtools
install.packages("devtools")
devtools::install_github("naimurashid/evolveTrial")
```

### Dependencies

The package requires:

- **Core**: `data.table`, `ggplot2`, `magrittr`, `progress`,
  `tidyselect`
- **Suggested**: `dplyr`, `gt`, `openxlsx`, `ggrepel`, `knitr`,
  `rmarkdown`

All dependencies will be installed automatically when you install
evolveTrial.

### Development Installation

For development work with the latest features:

``` r
# Clone the repository
git clone https://github.com/naimurashid/evolveTrial.git

# Load in R session
devtools::load_all("path/to/evolveTrial")
```

### Troubleshooting Installation

**“Error: compilation failed”** or **“make: not found”** - Windows:
Rtools is not installed or not on PATH. Reinstall Rtools and restart
R. - macOS: Xcode Command Line Tools not installed. Run
`xcode-select --install` in Terminal.

**“fatal error: ‘RcppArmadillo.h’ file not found”**

``` r
# Install RcppArmadillo first
install.packages("RcppArmadillo")
# Then retry evolveTrial installation
```

**“ld: library not found for -lgfortran”** (macOS)

``` bash
# Install gfortran via Homebrew
brew install gcc
```

**Installation hangs or times out**

``` r
# Try installing without vignettes (faster)
remotes::install_github("naimurashid/evolveTrial", build_vignettes = FALSE)
```

**Still having issues?** - Check that your R version is 4.0 or higher:
`R.version.string` - Open an issue at:
<https://github.com/naimurashid/evolveTrial/issues>

## Diagnosing interim gating choices

When early stopping does not seem to occur during simulations, it is
often because the information gates (minimum events, median follow-up,
or required person-time) postpone the first informative interim look
until late in the trial. The helper
[`estimate_vsref_gate_timing()`](reference/estimate_vsref_gate_timing.md)
gives quick lower-bound heuristics for when those gates can all be
satisfied under deterministic accrual.

``` r
library(evolveTrial)

args <- list(
  arm_names = c("Doublet", "Triplet"),
  reference_arm_name = "Doublet",
  overall_accrual_rate = 3,
  randomization_probs = c(Doublet = 0.5, Triplet = 0.5),
  max_total_patients_per_arm = c(Doublet = 70, Triplet = 70),
  max_follow_up_sim = 24,
  min_events_per_arm = 8,
  min_median_followup_per_arm = 3,
  min_person_time_frac_per_arm = 0.15
)

estimate_vsref_gate_timing(args)
```

For the configuration above the joint lower bound is roughly 18 months,
meaning that with 3-month calendar beats the first eligible interim look
will occur near the sixth look. Adjusting the information gates (for
example by lowering the person-time fraction or the minimum median
follow-up) yields earlier opportunities for efficacy/futility stopping.

## Minimal Working Examples

### Example 1: Multi-arm Trial (Experimental vs Reference)

Simulate a two-arm trial comparing Triplet therapy against Doublet
(reference) with interim efficacy and futility stopping:

``` r
library(evolveTrial)

# Define trial arguments
args <- list(
  # Arm configuration
  arm_names = c("Doublet", "Triplet"),
  reference_arm_name = "Doublet",
  randomization_probs = c(Doublet = 0.5, Triplet = 0.5),
  overall_accrual_rate = 3,

  # Sample size limits
  max_total_patients_per_arm = c(Doublet = 70, Triplet = 70),
  max_follow_up_sim = 24,

  # Piecewise hazards (8 3-month intervals)
  pw_hazards = list(
    Doublet = rep(0.12, 8),  # Median ~6 months
    Triplet = rep(0.077, 8)  # Median ~9 months (HR = 0.67)
  ),
  interval_width = 3,

  # Bayesian decision thresholds
  efficacy_prob_threshold_vsref = 0.99,
  futility_prob_threshold_vsref = 0.10,

  # Information gates
  min_events_per_arm = 8,
  min_median_followup_per_arm = 3,
  min_person_time_frac_per_arm = 0.15,

  # Interim schedule
  beat_interval = 3,

  # MCMC settings
  n_posterior_draws = 5000
)

# Run Monte Carlo simulation (100 replicates for demo)
result <- run_single_scenario(args, n_sim = 100, seed = 2025)

# View key operating characteristics
summary(result)
#> Power: ~80%, Type I error: ~10%, Expected N: ~85 per arm
```

### Example 2: Single-arm vs Historical Control

Evaluate a single experimental arm against a historical median survival:

``` r
library(evolveTrial)

# Define trial arguments for single-arm design
args_sa <- list(
  # Arm configuration (vs historical control)
  arm_names = c("Control", "Experimental"),
  reference_arm_name = "Control",
  use_historical_control = TRUE,
  randomization_probs = c(Control = 0, Experimental = 1),
  overall_accrual_rate = 2,

  # Sample size
  max_total_patients_per_arm = c(Control = 0, Experimental = 60),
  max_follow_up_sim = 24,

  # Hazards: Historical median 6 months, experimental median 9 months
  pw_hazards = list(
    Control = rep(0.116, 8),      # Historical: median 6 mo
    Experimental = rep(0.077, 8)  # Expected: median 9 mo
  ),
  interval_width = 3,

  # Decision thresholds vs historical
  efficacy_prob_threshold_hc = 0.95,
  futility_prob_threshold_hc = 0.05,

  # Information gates
  min_events_hc = 15,

  # Interim schedule
  beat_interval = 3,
  n_posterior_draws = 5000
)

# Run simulation
result_sa <- run_single_scenario(args_sa, n_sim = 100, seed = 2025)

# View operating characteristics
summary(result_sa)
```

### Example 3: Grid Search for Optimal Thresholds

Find efficacy/futility thresholds that meet operating characteristic
targets:

``` r
library(evolveTrial)

# Base trial configuration
args <- list(
  arm_names = c("Control", "Experimental"),
  reference_arm_name = "Control",
  use_historical_control = TRUE,
  randomization_probs = c(Control = 0, Experimental = 1),
  overall_accrual_rate = 2,
  max_total_patients_per_arm = c(Control = 0, Experimental = 60),
  max_follow_up_sim = 24,
  pw_hazards = list(
    Control = rep(0.116, 8),
    Experimental = rep(0.077, 8)
  ),
  interval_width = 3,
  min_events_hc = 15,
  beat_interval = 3,
  n_posterior_draws = 5000
)

# Define threshold grid
grid <- expand.grid(
  efficacy_prob_threshold_hc = c(0.90, 0.95, 0.99),
  futility_prob_threshold_hc = c(0.05, 0.10, 0.15)
)

# Evaluate grid (use evaluate_hc_grid for historical control)
results <- evaluate_hc_grid(
  base_args = args,
  grid = grid,
  n_sim = 500,
  parallel = TRUE,
  seed = 2025
)

# Filter feasible designs (type I <= 0.10, power >= 0.80)
feasible <- filter_feasible_designs(
  results,
  alpha_cap = 0.10,
  power_floor = 0.80
)

# Recommend optimal design (minimizes expected N)
best <- recommend_design(feasible)
print(best)
```

## Key Functions

| Function                                                                  | Purpose                                          |
|---------------------------------------------------------------------------|--------------------------------------------------|
| `run_single_scenario()`                                                   | Simulate a single trial configuration            |
| [`run_scenarios()`](reference/run_scenarios.md)                           | Simulate multiple scenario configurations        |
| `evaluate_hc_grid()`                                                      | Grid search for single-arm vs historical control |
| [`evaluate_ph_grid()`](reference/evaluate_ph_grid.md)                     | Grid search for multi-arm proportional hazards   |
| `filter_feasible_designs()`                                               | Filter designs meeting alpha/power constraints   |
| `recommend_design()`                                                      | Select optimal design from feasible set          |
| [`estimate_vsref_gate_timing()`](reference/estimate_vsref_gate_timing.md) | Diagnose information gate timing                 |
| `summarise_arm_level_OCs()`                                               | Summarize operating characteristics              |

## Package Structure

    evolveTrial/
    ├── R/
    │   ├── simulation_driver.R      # Main simulation functions
    │   ├── design_analysis.R        # Grid search and optimization
    │   ├── gate_diagnostics.R       # Information gate utilities
    │   ├── interim_logic.R          # Stopping rule implementation
    │   ├── posterior_helpers.R      # Bayesian posterior computations
    │   └── data_generation.R        # Time-to-event data generation
    ├── src/
    │   └── *.cpp                    # C++ acceleration code
    └── vignettes/
        └── design-overview.Rmd      # Detailed methodology

## Related Packages

- **evolveBO**: Bayesian optimization for calibrating evolveTrial
  designs
- **adaptive-trial-bo-paper**: Research manuscript and case studies
