# Run hybrid simulations with Rcpp

Batch simulation interface that uses the Rcpp implementation.

## Usage

``` r
run_hybrid_simulations_rcpp(
  hybrid_theta,
  base_args,
  scenario_params,
  num_simulations = 1000,
  seed = NULL,
  return_raw = FALSE,
  trial_mode = "hybrid",
  efficacy_method = "posterior",
  futility_method = "posterior"
)
```

## Arguments

- hybrid_theta:

  Hybrid design parameters

- base_args:

  Base simulation arguments

- scenario_params:

  Scenario parameters

- num_simulations:

  Number of replications

- seed:

  Random seed

- return_raw:

  If TRUE, return raw DataFrame; if FALSE, aggregate

- trial_mode:

  Trial mode: "hybrid", "single_arm", "between_arm", or
  "dual_single_arm"

- efficacy_method:

  Decision method for efficacy: "posterior" or "predictive"

- futility_method:

  Decision method for futility: "posterior" or "predictive"

## Value

DataFrame with simulation results or aggregated OC
