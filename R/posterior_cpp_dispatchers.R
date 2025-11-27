# posterior_cpp_dispatchers.R
# Dispatcher functions that route to C++ or R implementations

#' Check if C++ posterior sampling should be used
#'
#' Checks environment variable EVOLVETRIAL_USE_CPP to determine
#' whether to use C++ implementation (default: TRUE)
#'
#' @return Logical indicating whether to use C++ version
#' @keywords internal
.use_cpp_posterior <- function() {
  use_cpp <- Sys.getenv("EVOLVETRIAL_USE_CPP", "TRUE")
  return(toupper(use_cpp) %in% c("TRUE", "1", "YES"))
}

#' Draw posterior hazard samples (DISPATCHER)
#'
#' Routes to C++ or R implementation based on EVOLVETRIAL_USE_CPP.
#' Default: C++ for better performance.
#'
#' @param num_intervals Number of intervals
#' @param events_per_interval Events per interval
#' @param person_time_per_interval Person-time per interval
#' @param prior_alpha_params Gamma prior shapes
#' @param prior_beta_params Gamma prior rates
#' @param num_samples Number of samples to draw
#' @return Matrix of posterior hazard samples
#' @export
draw_posterior_hazard_samples <- function(
    num_intervals,
    events_per_interval,
    person_time_per_interval,
    prior_alpha_params,
    prior_beta_params,
    num_samples = 1000
) {
  if (.use_cpp_posterior()) {
    # Use C++ version (faster)
    draw_posterior_hazard_samples_cpp(
      num_intervals = as.integer(num_intervals),
      events_per_interval = as.integer(events_per_interval),
      person_time_per_interval = as.numeric(person_time_per_interval),
      prior_alpha_params = as.numeric(prior_alpha_params),
      prior_beta_params = as.numeric(prior_beta_params),
      num_samples = as.integer(num_samples)
    )
  } else {
    # Use original R version (for debugging/fallback)
    draw_posterior_hazard_samples_r(
      num_intervals, events_per_interval, person_time_per_interval,
      prior_alpha_params, prior_beta_params, num_samples
    )
  }
}

#' Draw posterior hazard samples (ORIGINAL R VERSION)
#'
#' @keywords internal
draw_posterior_hazard_samples_r <- function(
    num_intervals,
    events_per_interval,
    person_time_per_interval,
    prior_alpha_params,
    prior_beta_params,
    num_samples = 1000
) {
  posterior_hazard_samples <- matrix(NA, nrow = num_samples, ncol = num_intervals)

  for (j in 1:num_intervals) {
    posterior_alpha_j <- prior_alpha_params[j] + events_per_interval[j]
    posterior_beta_j <- prior_beta_params[j] + person_time_per_interval[j]

    if (posterior_beta_j <= 0) {
      warning(paste("Posterior beta parameter for interval", j, "is non-positive."))
      posterior_beta_j <- 1e-6
    }

    posterior_hazard_samples[, j] <- rgamma(num_samples, shape = posterior_alpha_j, rate = posterior_beta_j)
  }

  return(posterior_hazard_samples)
}

#' Compute mode and variance of log-HR posterior for PH model (DISPATCHER)
#'
#' @export
ph_beta_mode_var <- function(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, mu, sigma,
                              tol = 1e-6, max_iter = 50) {
  if (.use_cpp_posterior()) {
    # Use C++ version
    ph_beta_mode_var_cpp(
      E_C = as.integer(E_C),
      PT_C = as.numeric(PT_C),
      E_T = as.integer(E_T),
      PT_T = as.numeric(PT_T),
      alpha_prior = as.numeric(alpha_prior),
      beta_prior = as.numeric(beta_prior),
      mu = mu,
      sigma = sigma,
      tol = tol,
      max_iter = as.integer(max_iter)
    )
  } else {
    # Use original R version
    ph_beta_mode_var_r(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, mu, sigma, tol, max_iter)
  }
}

#' Compute mode and variance of log-HR posterior (ORIGINAL R VERSION)
#'
#' @keywords internal
ph_beta_mode_var_r <- function(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, mu, sigma,
                                tol = 1e-6, max_iter = 50) {
  total_events_ctrl <- sum(E_C)
  total_events_trt  <- sum(E_T)
  total_pt_ctrl <- sum(PT_C)
  total_pt_trt  <- sum(PT_T)
  beta <- log((total_events_trt + 0.5) / (total_pt_trt + 0.5)) -
    log((total_events_ctrl + 0.5) / (total_pt_ctrl + 0.5))

  sigma2 <- sigma^2

  for (iter in seq_len(max_iter)) {
    exp_beta <- exp(beta)
    denom <- beta_prior + PT_C + exp_beta * PT_T
    g <- exp_beta * PT_T / denom
    grad <- sum(E_T) - sum((alpha_prior + E_C + E_T) * g) - (beta - mu) / sigma2
    hess <- -sum((alpha_prior + E_C + E_T) * (g - g^2)) - 1 / sigma2
    step <- grad / hess
    beta_new <- beta - step

    if (is.nan(beta_new) || is.infinite(beta_new)) break
    beta <- beta_new
    if (abs(step) < tol) break
  }

  exp_beta <- exp(beta)
  denom <- beta_prior + PT_C + exp_beta * PT_T
  g <- exp_beta * PT_T / denom
  hess <- -sum((alpha_prior + E_C + E_T) * (g - g^2)) - 1 / sigma2
  var_beta <- if (hess < 0) min(1e6, max(1e-6, -1 / hess)) else 1

  list(mean = beta, sd = sqrt(var_beta))
}

#' Sample medians for vs-ref independent model (DISPATCHER)
#'
#' @export
sample_vs_ref_medians_independent <- function(slCtrl, slTrt, args, num_samples) {
  if (.use_cpp_posterior()) {
    # Use C++ version
    interval_lengths <- diff(args$interval_cutpoints_sim)
    E_C <- slCtrl$metrics$events_per_interval
    PT_C <- slCtrl$metrics$person_time_per_interval
    E_T <- slTrt$metrics$events_per_interval
    PT_T <- slTrt$metrics$person_time_per_interval
    alpha_prior <- args$prior_alpha_params_model
    beta_prior  <- args$prior_beta_params_model
    num_samples <- num_samples %||% args$num_posterior_draws

    sample_vs_ref_medians_independent_cpp(
      E_C = as.integer(E_C),
      PT_C = as.numeric(PT_C),
      E_T = as.integer(E_T),
      PT_T = as.numeric(PT_T),
      alpha_prior = as.numeric(alpha_prior),
      beta_prior = as.numeric(beta_prior),
      interval_lengths = as.numeric(interval_lengths),
      num_samples = as.integer(num_samples)
    )
  } else {
    # Use original R version
    sample_vs_ref_medians_independent_r(slCtrl, slTrt, args, num_samples)
  }
}

#' Sample medians for vs-ref independent model (ORIGINAL R VERSION)
#'
#' @keywords internal
sample_vs_ref_medians_independent_r <- function(slCtrl, slTrt, args, num_samples) {
  interval_lengths <- diff(args$interval_cutpoints_sim)
  num_samples <- num_samples %||% args$num_posterior_draws

  lamC <- draw_posterior_hazard_samples_r(
    num_intervals = length(interval_lengths),
    events_per_interval = slCtrl$metrics$events_per_interval,
    person_time_per_interval = slCtrl$metrics$person_time_per_interval,
    prior_alpha_params = args$prior_alpha_params_model,
    prior_beta_params  = args$prior_beta_params_model,
    num_samples = num_samples
  )
  lamT <- draw_posterior_hazard_samples_r(
    num_intervals = length(interval_lengths),
    events_per_interval = slTrt$metrics$events_per_interval,
    person_time_per_interval = slTrt$metrics$person_time_per_interval,
    prior_alpha_params = args$prior_alpha_params_model,
    prior_beta_params  = args$prior_beta_params_model,
    num_samples = num_samples
  )

  # PERFORMANCE: Use C++ matrix version instead of apply() for 20-30x speedup
  med_ctrl <- calculate_median_survival_matrix_cpp(lamC, interval_lengths)
  med_trt  <- calculate_median_survival_matrix_cpp(lamT, interval_lengths)

  list(medCtrl = med_ctrl, medTrt = med_trt, logHR = NULL)
}

#' Sample medians for vs-ref PH model (DISPATCHER)
#'
#' @export
sample_vs_ref_medians_ph <- function(slCtrl, slTrt, args, num_samples) {
  if (.use_cpp_posterior()) {
    # Use C++ version
    interval_lengths <- diff(args$interval_cutpoints_sim)
    E_C <- slCtrl$metrics$events_per_interval
    PT_C <- slCtrl$metrics$person_time_per_interval
    E_T <- slTrt$metrics$events_per_interval
    PT_T <- slTrt$metrics$person_time_per_interval
    alpha_prior <- args$prior_alpha_params_model
    beta_prior  <- args$prior_beta_params_model
    mu <- args$ph_loghr_prior_mean %||% 0
    sigma <- max(1e-6, args$ph_loghr_prior_sd %||% 1)
    num_samples <- num_samples %||% args$num_posterior_draws

    sample_vs_ref_medians_ph_cpp(
      E_C = as.integer(E_C),
      PT_C = as.numeric(PT_C),
      E_T = as.integer(E_T),
      PT_T = as.numeric(PT_T),
      alpha_prior = as.numeric(alpha_prior),
      beta_prior = as.numeric(beta_prior),
      interval_lengths = as.numeric(interval_lengths),
      mu = mu,
      sigma = sigma,
      num_samples = as.integer(num_samples)
    )
  } else {
    # Use original R version
    sample_vs_ref_medians_ph_r(slCtrl, slTrt, args, num_samples)
  }
}

#' Sample medians for vs-ref PH model (ORIGINAL R VERSION)
#'
#' @keywords internal
sample_vs_ref_medians_ph_r <- function(slCtrl, slTrt, args, num_samples) {
  interval_lengths <- diff(args$interval_cutpoints_sim)
  E_C <- slCtrl$metrics$events_per_interval
  PT_C <- slCtrl$metrics$person_time_per_interval
  E_T <- slTrt$metrics$events_per_interval
  PT_T <- slTrt$metrics$person_time_per_interval
  alpha_prior <- args$prior_alpha_params_model
  beta_prior  <- args$prior_beta_params_model
  beta_prior_mean <- args$ph_loghr_prior_mean %||% 0
  beta_prior_sd   <- max(1e-6, args$ph_loghr_prior_sd %||% 1)
  num_samples <- num_samples %||% args$num_posterior_draws

  beta_params <- ph_beta_mode_var_r(
    E_C = E_C, PT_C = PT_C,
    E_T = E_T, PT_T = PT_T,
    alpha_prior = alpha_prior,
    beta_prior  = beta_prior,
    mu = beta_prior_mean,
    sigma = beta_prior_sd
  )
  beta_draws <- rnorm(num_samples, mean = beta_params$mean, sd = beta_params$sd)
  shapes <- alpha_prior + E_C + E_T
  med_ctrl <- med_trt <- numeric(num_samples)

  for (i in seq_len(num_samples)) {
    exp_beta <- exp(beta_draws[i])
    rates <- beta_prior + PT_C + exp_beta * PT_T
    lambda <- rgamma(length(shapes), shape = shapes, rate = rates)
    med_ctrl[i] <- calculate_median_survival_piecewise(lambda, interval_lengths)
    med_trt[i]  <- calculate_median_survival_piecewise(lambda * exp_beta, interval_lengths)
  }

  list(medCtrl = med_ctrl, medTrt = med_trt, logHR = beta_draws)
}

