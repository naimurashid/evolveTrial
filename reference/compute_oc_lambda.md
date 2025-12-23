# Compute OC with direct lambda (hazard rate) input

Lower-level interface for users who want to specify piecewise
exponential hazard rates directly rather than using median-based
conversion.

## Usage

``` r
compute_oc_lambda(
  theta,
  base_args,
  lambda_exp,
  lambda_ref,
  lambda_hist,
  lambda_exp_null = NULL,
  num_simulations = 1000,
  seed = NULL,
  trial_mode = "hybrid",
  efficacy_method = "posterior",
  futility_method = "posterior"
)
```

## Arguments

- theta:

  Named list with design parameters (eff_sa, fut_sa, etc.)

- base_args:

  Named list with interval_cutpoints_sim, overall_accrual_rate, etc.

- lambda_exp:

  Hazard rates for experimental arm (vector, length = n_intervals)

- lambda_ref:

  Hazard rates for reference arm (vector, length = n_intervals)

- lambda_hist:

  Hazard rates for historical control (vector, length = n_intervals)

- lambda_exp_null:

  Optional hazard rates for experimental under null (default =
  lambda_ref)

- num_simulations:

  Number of MC replications

- seed:

  Random seed

- trial_mode:

  Trial mode: "hybrid", "single_arm", "between_arm", "dual_single_arm"

- efficacy_method:

  "posterior" or "predictive"

- futility_method:

  "posterior" or "predictive"

## Value

Named list with operating characteristics

## Examples

``` r
if (FALSE) { # \dontrun{
# 3-interval PWE model
base_args <- list(
  interval_cutpoints_sim = c(0, 6, 12, 24),
  overall_accrual_rate = 2,
  max_follow_up_sim = 12,
  prior_alpha_params_model = rep(0.5, 3),
  prior_beta_params_model = rep(0.5, 3)
)

theta <- list(
  eff_sa = 0.90, fut_sa = 0.10,
  eff_ba = 0.90, fut_ba = 0.10,
  ev_sa = 10, ev_ba = 20,
  nmax_sa = 40, nmax_ba = 80,
  hr_threshold_sa = 0.7,
  pp_go = 0.5, pp_nogo = 0.3
)

# Single-arm trial vs historical
oc <- compute_oc_lambda(
  theta, base_args,
  lambda_exp = c(0.04, 0.05, 0.06),
  lambda_ref = c(0.08, 0.09, 0.10),
  lambda_hist = c(0.05, 0.06, 0.07),
  trial_mode = "single_arm",
  num_simulations = 1000
)
} # }
```
