# Run a full set of evolveTrial simulations for a design specification

`run_simulation_pure()` is the workhorse simulator for evolveTrial. It
enrols patients according to the supplied accrual plan, applies interim
gating/decision rules, optionally rebalances interval cut points, and
carries each replicate forward to the final analysis. Operating
characteristics are returned summarised per arm.

## Usage

``` r
run_simulation_pure(
  num_simulations,
  arm_names,
  reference_arm_name,
  compare_arms_option,
  weibull_shape_true_arms,
  weibull_median_true_arms,
  null_median_arms = NULL,
  futility_median_arms = NULL,
  interval_cutpoints_sim,
  max_follow_up_sim,
  censor_max_time_sim,
  prior_alpha_params_model,
  prior_beta_params_model,
  num_posterior_draws,
  num_posterior_draws_interim = NULL,
  cohort_size_per_arm,
  max_total_patients_per_arm,
  min_patients_for_analysis = NULL,
  efficacy_stopping_rule_hc = FALSE,
  efficacy_threshold_current_prob_hc = NULL,
  posterior_futility_threshold_hc = NULL,
  efficacy_threshold_hc_prob = NULL,
  futility_threshold_hc_prob = NULL,
  futility_stopping_rule_hc = FALSE,
  efficacy_stopping_rule_vs_ref = FALSE,
  futility_stopping_rule_vs_ref = FALSE,
  efficacy_threshold_vs_ref_prob = NULL,
  futility_threshold_vs_ref_prob = NULL,
  compare_arms_futility_margin = 0,
  compare_arms_hr_margin = NULL,
  use_ph_model_vs_ref = FALSE,
  ph_loghr_prior_mean = 0,
  ph_loghr_prior_sd = 1,
  median_pfs_success_threshold_arms = NULL,
  final_success_posterior_prob_threshold = 0.85,
  median_pfs_futility_threshold_arms = NULL,
  final_futility_posterior_prob_threshold = 0.85,
  overall_accrual_rate,
  randomization_probs,
  min_follow_up_at_final = 0,
  min_events_for_analysis = NULL,
  min_median_followup = NULL,
  min_events_hc = NULL,
  min_median_followup_hc = NULL,
  interim_calendar_beat = 2,
  diagnostics = FALSE,
  pred_success_pp_threshold_hc = 1,
  pred_futility_pp_threshold_hc = 0,
  num_posterior_draws_pred = 100,
  predictive_fast = FALSE,
  min_events_per_arm = NULL,
  min_median_followup_per_arm = NULL,
  min_person_time_frac_per_arm = 0,
  person_time_milestones = NULL,
  latest_calendar_look = Inf,
  rebalance_after_events = NULL,
  parallel_replicates = FALSE,
  num_workers = NULL,
  cluster_type = c("auto", "PSOCK", "FORK"),
  cluster = NULL,
  progress = interactive()
)
```

## Arguments

- num_simulations:

  Number of Monte Carlo replicates to run.

- arm_names:

  Character vector naming the trial arms.

- reference_arm_name:

  Character scalar naming the control/reference arm.

- compare_arms_option:

  Logical; `TRUE` evaluates vs-reference logic, `FALSE` evaluates arms
  independently against historical control targets.

- weibull_shape_true_arms:

  Named numeric vector of Weibull shape parameters under the truth.

- weibull_median_true_arms:

  Named numeric vector of true median PFS (months) for each arm.

- null_median_arms:

  Named numeric vector of null (historical control) medians used for
  single-arm evaluations.

- futility_median_arms:

  Named numeric vector of futility medians for single-arm logic.

- interval_cutpoints_sim:

  Numeric vector of interval boundaries (months) for piecewise
  exponential modelling.

- max_follow_up_sim:

  Maximum administrative follow-up time (months).

- censor_max_time_sim:

  Upper bound for random censoring draws (months).

- prior_alpha_params_model:

  Numeric vector of Gamma prior shape parameters for the piecewise
  exponential hazards.

- prior_beta_params_model:

  Numeric vector of Gamma prior rate parameters for the piecewise
  exponential hazards.

- num_posterior_draws:

  Number of posterior draws used for final analyses.

- num_posterior_draws_interim:

  Optional integer overriding the number of posterior draws used at
  interim looks.

- cohort_size_per_arm:

  Size of each enrolment batch per arm (typically 1).

- max_total_patients_per_arm:

  Named integer vector of per-arm sample size caps.

- min_patients_for_analysis:

  Minimum number of patients required to evaluate an arm in the
  single-arm path. If not specified, defaults to 0, allowing interim
  analyses to occur even with very few patients. Set this to a higher
  value to prevent interim analyses until a certain number of patients
  have been enrolled.

- efficacy_stopping_rule_hc:

  Logical; enable interim efficacy checks for the historical-control
  path.

- efficacy_threshold_current_prob_hc:

  **DEPRECATED**. Use `efficacy_threshold_hc_prob` instead. Interim
  success probability threshold for single-arm logic.

- posterior_futility_threshold_hc:

  **DEPRECATED**. Use `futility_threshold_hc_prob` instead. Interim
  futility probability threshold for single-arm logic.

- efficacy_threshold_hc_prob:

  Interim success probability threshold for single-arm logic (preferred
  harmonized name).

- futility_threshold_hc_prob:

  Interim futility probability threshold for single-arm logic (preferred
  harmonized name).

- futility_stopping_rule_hc:

  Logical; enable interim futility checks for the single-arm path.

- efficacy_stopping_rule_vs_ref:

  Logical; enable interim efficacy checks for the vs-reference path.

- futility_stopping_rule_vs_ref:

  Logical; enable interim futility checks for the vs-reference path.

- efficacy_threshold_vs_ref_prob:

  Posterior superiority threshold for vs-reference decisions.

- futility_threshold_vs_ref_prob:

  Posterior inferiority threshold for vs-reference decisions.

- compare_arms_futility_margin:

  Absolute median difference used when defining vs-reference futility.

- compare_arms_hr_margin:

  Optional hazard-ratio margin used when `use_ph_model_vs_ref = TRUE`.

- use_ph_model_vs_ref:

  Logical; use the proportional-hazards joint model for vs-reference
  comparisons.

- ph_loghr_prior_mean:

  Mean of the normal prior on the log hazard ratio (PH model).

- ph_loghr_prior_sd:

  Standard deviation of the normal prior on the log hazard ratio.

- median_pfs_success_threshold_arms:

  Named numeric vector of median PFS thresholds for declaring final
  success per arm.

- final_success_posterior_prob_threshold:

  Posterior probability threshold for final success declarations.

- median_pfs_futility_threshold_arms:

  Named numeric vector of median PFS futility thresholds for final
  analyses.

- final_futility_posterior_prob_threshold:

  Posterior probability threshold for final futility declarations.

- overall_accrual_rate:

  Expected accrual rate (patients per month).

- randomization_probs:

  Named numeric vector of randomisation probabilities.

- min_follow_up_at_final:

  Additional follow-up (months) required after last enrolment before the
  final analysis.

- min_events_for_analysis:

  **DEPRECATED**. Use `min_events_hc` instead. Minimum events required
  for interim review (global gate).

- min_median_followup:

  **DEPRECATED**. Use `min_median_followup_hc` instead. Minimum median
  follow-up required for interim review (global gate).

- min_events_hc:

  Minimum events required for single-arm interim review (preferred
  harmonized name).

- min_median_followup_hc:

  Minimum median follow-up required for single-arm interim review
  (preferred harmonized name).

- interim_calendar_beat:

  Calendar spacing (months) between scheduled interim looks when
  person-time milestones are not used.

- diagnostics:

  Logical; if `TRUE` prints interim diagnostic messages.

- pred_success_pp_threshold_hc:

  Predictive probability threshold for interim success in the single-arm
  predictive look (if enabled).

- pred_futility_pp_threshold_hc:

  Predictive probability threshold for interim futility in the
  single-arm predictive look (if enabled).

- num_posterior_draws_pred:

  Number of posterior draws used inside predictive probability
  calculations.

- predictive_fast:

  Logical; switch to the analytic predictive approximations.

- min_events_per_arm:

  Optional per-arm minimum event gate for vs-reference.

- min_median_followup_per_arm:

  Optional per-arm minimum median follow-up gate for vs-reference.

- min_person_time_frac_per_arm:

  Optional per-arm proportion of planned person-time required before
  evaluating vs-reference decisions.

- person_time_milestones:

  Optional numeric vector (fractions of total planned person-time)
  triggering interim looks.

- latest_calendar_look:

  Backstop calendar time for person-time schedules.

- rebalance_after_events:

  Optional integer; when non-`NULL` the piecewise cut points are
  re-estimated once that number of events has accrued.

- parallel_replicates:

  Logical; if `TRUE`, distribute Monte Carlo replicates across a
  parallel cluster.

- num_workers:

  Optional integer specifying the number of workers when
  `parallel_replicates = TRUE`. Defaults to
  `parallel::detectCores() - 1`.

- cluster_type:

  Type of parallel cluster to spawn when distributing replicates. One of
  `"auto"` (default, uses FORK on Unix, PSOCK on Windows), `"PSOCK"`, or
  `"FORK"`. FORK clusters are faster on Linux/macOS as they share memory
  and don't require package loading in workers.

- cluster:

  Optional pre-existing parallel cluster to reuse. When provided, the
  cluster is used for parallel execution and NOT stopped on exit. This
  enables cluster pooling for repeated calls (e.g., in Bayesian
  optimization). Create with
  [`evolveTrial::create_simulation_cluster()`](create_simulation_cluster.md)
  or
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).

- progress:

  Logical; show the simulation progress bar when running sequentially.
  Automatically disabled for parallel replicate execution.

## Value

A data frame with one row per arm and columns summarising operating
characteristics such as type I error / power, PETs, final decision
probabilities, and expected sample size.

## Details

When `parallel_replicates = TRUE`, results will vary based on
`num_workers` due to different random number stream partitioning. For
exact reproducibility across runs, use `parallel_replicates = FALSE`.

During package development with `devtools::load_all()`, parallel workers
will load the installed package version, not the development code. For
testing development changes, either use `parallel_replicates = FALSE` or
reinstall the package with `devtools::install()`.
