# Evaluate a design across multiple scenarios

Merges each scenario override list with `base_args`, runs
[`run_simulation_pure()`](run_simulation_pure.md) for every scenario,
and binds the results into a single table.

## Usage

``` r
run_scenarios(
  base_args,
  scens,
  parallel = FALSE,
  seed = NULL,
  return_percentiles = FALSE,
  percentile_probs = c(0, 0.25, 0.5, 0.75, 0.9, 1)
)
```

## Arguments

- base_args:

  Named list of arguments accepted by
  [`run_simulation_pure()`](run_simulation_pure.md).

- scens:

  List of scenario override lists, typically from
  [`scenarios_from_grid()`](scenarios_from_grid.md).

- parallel:

  Logical; if `TRUE` uses
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) to
  distribute scenarios across cores.

- seed:

  Optional integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before simulations.

- return_percentiles:

  Logical; if `TRUE`, collect per-replicate sample sizes and return
  percentile summaries. Default `FALSE`.

- percentile_probs:

  Numeric vector of probabilities for percentile computation when
  `return_percentiles = TRUE`. Default
  `c(0, 0.25, 0.5, 0.75, 0.9, 1.0)`.

## Value

When `return_percentiles = FALSE` (default), a data.table/data.frame
containing the combined operating characteristic summaries with a
`scenario` column identifying the originating scenario index. When
`return_percentiles = TRUE`, a list with:

- summary:

  The combined summary data.table as when `return_percentiles = FALSE`

- percentiles:

  A list of percentile results, one per scenario

## Examples

``` r
if (FALSE) { # \dontrun{
base_args <- list(
  num_simulations = 200,
  arm_names = c("Doublet", "Triplet"),
  reference_arm_name = "Doublet",
  compare_arms_option = TRUE,
  weibull_shape_true_arms = c(Doublet = 1.2, Triplet = 1.2),
  weibull_median_true_arms = c(Doublet = 6, Triplet = 6),
  null_median_arms = c(Doublet = 6, Triplet = 6),
  futility_median_arms = c(Doublet = 6, Triplet = 6),
  interval_cutpoints_sim = seq(0, 24, by = 3),
  max_follow_up_sim = 24,
  censor_max_time_sim = 24,
  prior_alpha_params_model = rep(0.5, 8),
  prior_beta_params_model = rep(0.5, 8),
  num_posterior_draws = 400,
  cohort_size_per_arm = 1,
  max_total_patients_per_arm = c(Doublet = 60, Triplet = 60),
  min_patients_for_analysis = 10,
  efficacy_stopping_rule_hc = TRUE,
  efficacy_threshold_current_prob_hc = 0.95,
  posterior_futility_threshold_hc = 0.8,
  futility_stopping_rule_hc = TRUE,
  efficacy_threshold_vs_ref_prob = 0.98,
  futility_threshold_vs_ref_prob = 0.6,
  compare_arms_futility_margin = 0.4,
  overall_accrual_rate = 3,
  randomization_probs = c(Doublet = 1/3, Triplet = 2/3),
  min_follow_up_at_final = 0,
  min_events_for_analysis = 0,
  min_median_followup = 0,
  interim_calendar_beat = 3,
  pred_success_pp_threshold_hc = 1,
  pred_futility_pp_threshold_hc = 0,
  num_posterior_draws_pred = 200
)

scens <- scenarios_from_grid(list(
  weibull_median_true_arms = list(
    c(Doublet = 6, Triplet = 6),
    c(Doublet = 6, Triplet = 9)
  )
))

run_scenarios(base_args, scens, parallel = FALSE, seed = 123)
} # }
```
