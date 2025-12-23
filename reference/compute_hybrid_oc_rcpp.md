# Compute operating characteristics using Rcpp

High-level wrapper that converts R-style parameters to the format
expected by compute_hybrid_oc_cpp(). Supports multiple trial modes:

- "hybrid" (default): SA→conversion→BA seamless design

- "single_arm": SA only, vs historical control

- "between_arm": BA only, randomized comparison

- "dual_single_arm": Two independent SAs (both arms reported)

## Usage

``` r
compute_hybrid_oc_rcpp(
  hybrid_theta,
  base_args,
  scenario_params,
  num_simulations = 1000,
  seed = NULL,
  trial_mode = "hybrid",
  efficacy_method = "posterior",
  futility_method = "posterior",
  lambda_hist_per_arm = NULL
)
```

## Arguments

- hybrid_theta:

  Named list with hybrid design parameters

- base_args:

  Named list with base simulation arguments

- scenario_params:

  Named list with scenario parameters (medians)

- num_simulations:

  Number of MC replications

- seed:

  Optional random seed

- trial_mode:

  Trial mode: "hybrid", "single_arm", "between_arm", or
  "dual_single_arm"

- efficacy_method:

  Decision method for efficacy: "posterior" or "predictive"

- futility_method:

  Decision method for futility: "posterior" or "predictive"

- lambda_hist_per_arm:

  Optional list with per-arm historical lambdas (for dual_single_arm)

## Value

Named list with operating characteristics (including per-arm metrics for
dual_single_arm)
