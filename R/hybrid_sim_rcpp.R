#' Rcpp-powered hybrid trial simulation wrapper
#'
#' Provides a high-level R interface to the C++ hybrid trial simulator,
#' converting median-based parameters to hazard rate parameters.
#'
#' @name hybrid_sim_rcpp
NULL

#' Compute operating characteristics using Rcpp
#'
#' High-level wrapper that converts R-style parameters to the format
#' expected by compute_hybrid_oc_cpp(). Supports multiple trial modes:
#' - "hybrid" (default): SA→conversion→BA seamless design
#' - "single_arm": SA only, vs historical control
#' - "between_arm": BA only, randomized comparison
#' - "dual_single_arm": Two independent SAs (both arms reported)
#'
#' @param hybrid_theta Named list with hybrid design parameters
#' @param base_args Named list with base simulation arguments
#' @param scenario_params Named list with scenario parameters (medians)
#' @param num_simulations Number of MC replications
#' @param seed Optional random seed
#' @param trial_mode Trial mode: "hybrid", "single_arm", "between_arm", or "dual_single_arm"
#' @param efficacy_method Decision method for efficacy: "posterior" or "predictive"
#' @param futility_method Decision method for futility: "posterior" or "predictive"
#' @param lambda_hist_per_arm Optional list with per-arm historical lambdas (for dual_single_arm)
#'
#' @return Named list with operating characteristics (including per-arm metrics for dual_single_arm)
#'
#' @export
compute_hybrid_oc_rcpp <- function(hybrid_theta, base_args, scenario_params,
                                    num_simulations = 1000, seed = NULL,
                                    trial_mode = "hybrid",
                                    efficacy_method = "posterior",
                                    futility_method = "posterior",
                                    lambda_hist_per_arm = NULL) {

  # Validate trial_mode
  valid_modes <- c("hybrid", "single_arm", "between_arm", "dual_single_arm")
  if (!trial_mode %in% valid_modes) {
    stop("trial_mode must be one of: ", paste(valid_modes, collapse = ", "))
  }

  # Validate decision methods
  valid_methods <- c("posterior", "predictive")
  if (!efficacy_method %in% valid_methods) {
    stop("efficacy_method must be 'posterior' or 'predictive'")
  }
  if (!futility_method %in% valid_methods) {
    stop("futility_method must be 'posterior' or 'predictive'")
  }

  # Extract intervals from base_args
  interval_cutpoints <- base_args$interval_cutpoints_sim

  # Validate interval_cutpoints
  if (is.null(interval_cutpoints) || length(interval_cutpoints) < 2) {
    stop("base_args$interval_cutpoints_sim must have at least 2 elements")
  }

  n_intervals <- length(interval_cutpoints) - 1

  # Validate required hybrid_theta parameters
  required_theta <- c("eff_sa", "fut_sa", "ev_sa", "nmax_sa", "pp_go", "pp_nogo",
                      "eff_ba", "fut_ba", "ev_ba", "nmax_ba")
  missing_theta <- setdiff(required_theta, names(hybrid_theta))
  if (length(missing_theta) > 0) {
    stop("hybrid_theta is missing required parameters: ", paste(missing_theta, collapse = ", "))
  }

  # Validate required scenario_params
  required_scenario <- c("historical_median", "ref_median", "exp_median")
  missing_scenario <- setdiff(required_scenario, names(scenario_params))
  if (length(missing_scenario) > 0) {
    stop("scenario_params is missing required parameters: ", paste(missing_scenario, collapse = ", "))
  }

  # Convert medians to piecewise exponential hazard rates
  # Using constant hazard per interval (exponential approximation)
  historical_median <- scenario_params$historical_median
  ref_median <- scenario_params$ref_median
  exp_median <- scenario_params$exp_median

  # Lambda = log(2) / median for exponential model (constant across intervals)
  lambda_hist <- rep(log(2) / historical_median, n_intervals)
  lambda_ref <- rep(log(2) / ref_median, n_intervals)
  lambda_exp_alt <- rep(log(2) / exp_median, n_intervals)
  lambda_exp_null <- lambda_ref  # Under null, exp has same hazard as ref

  # Build theta list for C++
  theta_cpp <- list(
    eff_sa = hybrid_theta$eff_sa,
    fut_sa = hybrid_theta$fut_sa,
    hr_threshold_sa = hybrid_theta$hr_threshold_sa %||% 0.8,
    ev_sa = as.integer(hybrid_theta$ev_sa),
    nmax_sa = as.integer(hybrid_theta$nmax_sa),
    conversion_trigger = hybrid_theta$conversion_trigger %||% "any_single_success",
    pp_go = hybrid_theta$pp_go,
    pp_nogo = hybrid_theta$pp_nogo,
    ss_method = hybrid_theta$ss_method %||% "posterior",
    max_additional_n = as.integer(hybrid_theta$max_additional_n %||% 60L),
    eff_ba = hybrid_theta$eff_ba,
    fut_ba = hybrid_theta$fut_ba,
    ev_ba = as.integer(hybrid_theta$ev_ba),
    nmax_ba = as.integer(hybrid_theta$nmax_ba),
    futility_action = hybrid_theta$futility_action %||% "drop_arm",
    prior_strength = hybrid_theta$prior_strength %||% 0.5,
    n_outer = as.integer(hybrid_theta$n_outer %||% 200),
    # New parameters for trial modes
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method
  )

  # Build base_args list for C++
  base_args_cpp <- list(
    n_intervals = n_intervals,
    interval_cutpoints_sim = interval_cutpoints,
    overall_accrual_rate = base_args$overall_accrual_rate %||% 2.0,
    max_follow_up_sim = base_args$max_follow_up_sim %||% 24,
    max_trial_time = base_args$max_trial_time %||% 72,
    prior_alpha_params_model = base_args$prior_alpha_params_model %||% rep(0.5, n_intervals),
    prior_beta_params_model = base_args$prior_beta_params_model %||% rep(0.5, n_intervals)
  )

  # Alternative scenario
  alt_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp_alt
  )

  # Null scenario (no treatment effect - null_global)
  null_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp_null
  )

  # Null between scenario (both arms better than historical, but equal to each other)
  # This tests type I error specifically for the between-arm comparison
  null_between_improvement <- scenario_params$null_between_improvement %||% 1.25
  null_between_median <- historical_median * null_between_improvement
  lambda_null_between <- rep(log(2) / null_between_median, n_intervals)
  null_between_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_null_between,
    lambda_exp = lambda_null_between
  )

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Call C++ OC computation for alt and null_global
  oc <- compute_hybrid_oc_cpp(
    n_sim = as.integer(num_simulations),
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = alt_scenario,
    null_scenario_list = null_scenario
  )

  # Run null_between scenario separately to get BA-specific type I error
  null_between_results <- run_hybrid_simulations_cpp(
    n_sim = as.integer(num_simulations),
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = null_between_scenario
  )

  # type1_between is the BA efficacy rate under null_between scenario
  # This is the false positive rate specifically for the between-arm comparison
  type1_between <- mean(null_between_results$ba_efficacy)

  # Return in expected format with per-arm metrics
  result <- list(
    power = oc$power,
    power_exp = oc$power_exp,
    power_ref = oc$power_ref,
    type1 = oc$type1,
    type1_exp = oc$type1_exp,
    type1_ref = oc$type1_ref,
    type1_between = type1_between,
    EN_null = oc$EN_null,
    EN_alt = oc$EN_alt,
    EN_alt_exp = oc$EN_alt_exp,
    EN_alt_ref = oc$EN_alt_ref,
    EN_null_exp = oc$EN_null_exp,
    EN_null_ref = oc$EN_null_ref,
    EN_total = 0.5 * oc$EN_null + 0.5 * oc$EN_alt,
    P_conversion = oc$conversion_rate,
    P_conversion_null = oc$conversion_rate_null,
    P_conversion_alt = oc$conversion_rate,
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method,
    var_power = oc$power * (1 - oc$power) / num_simulations,
    var_type1 = oc$type1 * (1 - oc$type1) / num_simulations,
    var_type1_between = type1_between * (1 - type1_between) / num_simulations,
    results_by_scenario = list(
      alternative = list(
        success_rate = oc$power,
        success_rate_exp = oc$power_exp,
        success_rate_ref = oc$power_ref,
        EN_total = oc$EN_alt,
        EN_exp = oc$EN_alt_exp,
        EN_ref = oc$EN_alt_ref,
        conversion_rate = oc$conversion_rate
      ),
      null_global = list(
        success_rate = oc$type1,
        success_rate_exp = oc$type1_exp,
        success_rate_ref = oc$type1_ref,
        EN_total = oc$EN_null,
        EN_exp = oc$EN_null_exp,
        EN_ref = oc$EN_null_ref,
        conversion_rate = oc$conversion_rate_null
      ),
      null_between = list(
        ba_success_rate = type1_between,
        EN_total = mean(null_between_results$total_n),
        conversion_rate = mean(null_between_results$converted)
      )
    )
  )

  # For dual_single_arm mode, add clarifying documentation
  if (trial_mode == "dual_single_arm") {
    result$note <- "power_exp/type1_exp = experimental arm, power_ref/type1_ref = reference arm"
  }

  result
}

#' Run hybrid simulations with Rcpp
#'
#' Batch simulation interface that uses the Rcpp implementation.
#'
#' @param hybrid_theta Hybrid design parameters
#' @param base_args Base simulation arguments
#' @param scenario_params Scenario parameters
#' @param num_simulations Number of replications
#' @param seed Random seed
#' @param return_raw If TRUE, return raw DataFrame; if FALSE, aggregate
#' @param trial_mode Trial mode: "hybrid", "single_arm", "between_arm", or "dual_single_arm"
#' @param efficacy_method Decision method for efficacy: "posterior" or "predictive"
#' @param futility_method Decision method for futility: "posterior" or "predictive"
#'
#' @return DataFrame with simulation results or aggregated OC
#'
#' @export
run_hybrid_simulations_rcpp <- function(hybrid_theta, base_args, scenario_params,
                                         num_simulations = 1000, seed = NULL,
                                         return_raw = FALSE,
                                         trial_mode = "hybrid",
                                         efficacy_method = "posterior",
                                         futility_method = "posterior") {

  # Validate trial_mode
  valid_modes <- c("hybrid", "single_arm", "between_arm", "dual_single_arm")
  if (!trial_mode %in% valid_modes) {
    stop("trial_mode must be one of: ", paste(valid_modes, collapse = ", "))
  }

  # Extract intervals
  interval_cutpoints <- base_args$interval_cutpoints_sim

  # Validate interval_cutpoints
  if (is.null(interval_cutpoints) || length(interval_cutpoints) < 2) {
    stop("base_args$interval_cutpoints_sim must have at least 2 elements")
  }

  n_intervals <- length(interval_cutpoints) - 1

  # Validate required hybrid_theta parameters
  required_theta <- c("eff_sa", "fut_sa", "ev_sa", "nmax_sa", "pp_go", "pp_nogo",
                      "eff_ba", "fut_ba", "ev_ba", "nmax_ba")
  missing_theta <- setdiff(required_theta, names(hybrid_theta))
  if (length(missing_theta) > 0) {
    stop("hybrid_theta is missing required parameters: ", paste(missing_theta, collapse = ", "))
  }

  # Validate required scenario_params
  required_scenario <- c("historical_median", "ref_median", "exp_median")
  missing_scenario <- setdiff(required_scenario, names(scenario_params))
  if (length(missing_scenario) > 0) {
    stop("scenario_params is missing required parameters: ", paste(missing_scenario, collapse = ", "))
  }

  # Convert medians to lambdas
  historical_median <- scenario_params$historical_median
  ref_median <- scenario_params$ref_median
  exp_median <- scenario_params$exp_median

  lambda_hist <- rep(log(2) / historical_median, n_intervals)
  lambda_ref <- rep(log(2) / ref_median, n_intervals)
  lambda_exp <- rep(log(2) / exp_median, n_intervals)

  # Build theta for C++
  theta_cpp <- list(
    eff_sa = hybrid_theta$eff_sa,
    fut_sa = hybrid_theta$fut_sa,
    hr_threshold_sa = hybrid_theta$hr_threshold_sa %||% 0.8,
    ev_sa = as.integer(hybrid_theta$ev_sa),
    nmax_sa = as.integer(hybrid_theta$nmax_sa),
    conversion_trigger = hybrid_theta$conversion_trigger %||% "any_single_success",
    pp_go = hybrid_theta$pp_go,
    pp_nogo = hybrid_theta$pp_nogo,
    ss_method = hybrid_theta$ss_method %||% "posterior",
    max_additional_n = as.integer(hybrid_theta$max_additional_n %||% 60L),
    eff_ba = hybrid_theta$eff_ba,
    fut_ba = hybrid_theta$fut_ba,
    ev_ba = as.integer(hybrid_theta$ev_ba),
    nmax_ba = as.integer(hybrid_theta$nmax_ba),
    futility_action = hybrid_theta$futility_action %||% "drop_arm",
    prior_strength = hybrid_theta$prior_strength %||% 0.5,
    n_outer = as.integer(hybrid_theta$n_outer %||% 200),
    # Trial mode parameters
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method
  )

  # Build base_args for C++
  base_args_cpp <- list(
    n_intervals = n_intervals,
    interval_cutpoints_sim = interval_cutpoints,
    overall_accrual_rate = base_args$overall_accrual_rate %||% 2.0,
    max_follow_up_sim = base_args$max_follow_up_sim %||% 24,
    max_trial_time = base_args$max_trial_time %||% 72,
    prior_alpha_params_model = base_args$prior_alpha_params_model %||% rep(0.5, n_intervals),
    prior_beta_params_model = base_args$prior_beta_params_model %||% rep(0.5, n_intervals)
  )

  # Scenario for C++
  scenario_cpp <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Call C++ batch simulation
  results <- run_hybrid_simulations_cpp(
    n_sim = as.integer(num_simulations),
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = scenario_cpp
  )

  if (return_raw) {
    return(results)
  }

  # Aggregate results
  # Guard against NaN when no conversions occurred
  n_converted <- sum(results$converted)
  ba_eff_rate <- if (n_converted > 0) {
    mean(results$ba_efficacy[results$converted])
  } else {
    NA_real_
  }

  list(
    power = mean(results$is_success),
    EN = mean(results$total_n),
    conversion_rate = mean(results$converted),
    sa_efficacy_rate = mean(results$sa_efficacy),
    ba_efficacy_rate = ba_eff_rate,
    outcome_distribution = table(results$outcome),
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method
  )
}

#' Compute OC with direct lambda (hazard rate) input
#'
#' Lower-level interface for users who want to specify piecewise exponential
#' hazard rates directly rather than using median-based conversion.
#'
#' @param theta Named list with design parameters (eff_sa, fut_sa, etc.)
#' @param base_args Named list with interval_cutpoints_sim, overall_accrual_rate, etc.
#' @param lambda_exp Hazard rates for experimental arm (vector, length = n_intervals)
#' @param lambda_ref Hazard rates for reference arm (vector, length = n_intervals)
#' @param lambda_hist Hazard rates for historical control (vector, length = n_intervals)
#' @param lambda_exp_null Optional hazard rates for experimental under null (default = lambda_ref)
#' @param num_simulations Number of MC replications
#' @param seed Random seed
#' @param trial_mode Trial mode: "hybrid", "single_arm", "between_arm", "dual_single_arm"
#' @param efficacy_method "posterior" or "predictive"
#' @param futility_method "posterior" or "predictive"
#'
#' @return Named list with operating characteristics
#'
#' @examples
#' \dontrun{
#' # 3-interval PWE model
#' base_args <- list(
#'   interval_cutpoints_sim = c(0, 6, 12, 24),
#'   overall_accrual_rate = 2,
#'   max_follow_up_sim = 12,
#'   prior_alpha_params_model = rep(0.5, 3),
#'   prior_beta_params_model = rep(0.5, 3)
#' )
#'
#' theta <- list(
#'   eff_sa = 0.90, fut_sa = 0.10,
#'   eff_ba = 0.90, fut_ba = 0.10,
#'   ev_sa = 10, ev_ba = 20,
#'   nmax_sa = 40, nmax_ba = 80,
#'   hr_threshold_sa = 0.7,
#'   pp_go = 0.5, pp_nogo = 0.3
#' )
#'
#' # Single-arm trial vs historical
#' oc <- compute_oc_lambda(
#'   theta, base_args,
#'   lambda_exp = c(0.04, 0.05, 0.06),
#'   lambda_ref = c(0.08, 0.09, 0.10),
#'   lambda_hist = c(0.05, 0.06, 0.07),
#'   trial_mode = "single_arm",
#'   num_simulations = 1000
#' )
#' }
#'
#' @export
compute_oc_lambda <- function(theta, base_args,
                               lambda_exp, lambda_ref, lambda_hist,
                               lambda_exp_null = NULL,
                               num_simulations = 1000, seed = NULL,
                               trial_mode = "hybrid",
                               efficacy_method = "posterior",
                               futility_method = "posterior") {

  # Validate trial_mode
  valid_modes <- c("hybrid", "single_arm", "between_arm", "dual_single_arm")
  if (!trial_mode %in% valid_modes) {
    stop("trial_mode must be one of: ", paste(valid_modes, collapse = ", "))
  }

  # Validate decision methods
  valid_methods <- c("posterior", "predictive")
  if (!efficacy_method %in% valid_methods) {
    stop("efficacy_method must be 'posterior' or 'predictive'")
  }
  if (!futility_method %in% valid_methods) {
    stop("futility_method must be 'posterior' or 'predictive'")
  }

  # Get interval info
  interval_cutpoints <- base_args$interval_cutpoints_sim

  # Validate interval_cutpoints
  if (is.null(interval_cutpoints) || length(interval_cutpoints) < 2) {
    stop("base_args$interval_cutpoints_sim must have at least 2 elements")
  }

  n_intervals <- length(interval_cutpoints) - 1

  # Validate lambda lengths
  if (length(lambda_exp) != n_intervals) {
    stop("lambda_exp must have length ", n_intervals)
  }
  if (length(lambda_ref) != n_intervals) {
    stop("lambda_ref must have length ", n_intervals)
  }
  if (length(lambda_hist) != n_intervals) {
    stop("lambda_hist must have length ", n_intervals)
  }

  # Default null = reference

  if (is.null(lambda_exp_null)) {
    lambda_exp_null <- lambda_ref
  }

  # Build theta list for C++
  theta_cpp <- list(
    eff_sa = theta$eff_sa,
    fut_sa = theta$fut_sa,
    hr_threshold_sa = theta$hr_threshold_sa %||% 0.7,
    ev_sa = as.integer(theta$ev_sa %||% 10),
    nmax_sa = as.integer(theta$nmax_sa %||% 40),
    conversion_trigger = theta$conversion_trigger %||% "any_single_success",
    k_required = as.integer(theta$k_required %||% 1),
    pp_go = theta$pp_go %||% 0.5,
    pp_nogo = theta$pp_nogo %||% 0.3,
    eff_ba = theta$eff_ba,
    fut_ba = theta$fut_ba,
    ev_ba = as.integer(theta$ev_ba %||% 20),
    nmax_ba = as.integer(theta$nmax_ba %||% 80),
    futility_action = theta$futility_action %||% "drop_arm",
    n_outer = as.integer(theta$n_outer %||% 200),
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method
  )

  # Build base_args for C++
  base_args_cpp <- list(
    n_intervals = n_intervals,
    interval_cutpoints_sim = interval_cutpoints,
    overall_accrual_rate = base_args$overall_accrual_rate %||% 2.0,
    max_follow_up_sim = base_args$max_follow_up_sim %||% 24,
    look_interval = base_args$look_interval %||% 3.0,
    max_trial_time = base_args$max_trial_time %||% 72,
    prior_alpha_params_model = base_args$prior_alpha_params_model %||% rep(0.5, n_intervals),
    prior_beta_params_model = base_args$prior_beta_params_model %||% rep(0.5, n_intervals)
  )

  # Alternative scenario
  alt_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp
  )

  # Null scenario
  null_scenario <- list(
    lambda_hist = lambda_hist,
    lambda_ref = lambda_ref,
    lambda_exp = lambda_exp_null
  )

  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)

  # Call C++ OC computation
  oc <- compute_hybrid_oc_cpp(
    n_sim = as.integer(num_simulations),
    theta_list = theta_cpp,
    base_args_list = base_args_cpp,
    scenario_params_list = alt_scenario,
    null_scenario_list = null_scenario
  )

  # Return with mode info
  c(oc, list(
    trial_mode = trial_mode,
    efficacy_method = efficacy_method,
    futility_method = futility_method
  ))
}
