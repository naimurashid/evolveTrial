
# --- HELPER FUNCTION: Calculate Median Survival for a Piecewise Exponential Model ---
# --- UPDATED: Calculate Median Survival for a Piecewise Exponential Model ---
calculate_median_survival_piecewise <- function(hazard_rates, interval_lengths) {
  #' Computes the median survival time for a piecewise exponential model.
  #' If the 0.5 survival is not reached by the end of the last interval,
  #' we *continue past the last cutpoint* with the last interval's hazard (open-ended tail).
  #' Only return Inf if the last hazard is exactly zero.
  
  if (length(hazard_rates) != length(interval_lengths)) {
    stop("hazard_rates and interval_lengths must have the same length.")
  }
  
  cumulative_hazard <- 0.0
  current_time <- 0.0
  
  for (i in seq_along(hazard_rates)) {
    lambda_j <- hazard_rates[i]
    delta_t_j <- interval_lengths[i]
    cumulative_hazard_end_interval <- cumulative_hazard + (lambda_j * delta_t_j)
    # Survival at end of interval:
    # S(t_end) = exp(- cumulative_hazard_end_interval)
    if (cumulative_hazard_end_interval >= log(2)) {
      # Median occurs *within* this interval
      if (lambda_j == 0) {
        # If hazard is zero but we've already crossed log(2) (numerically), clamp to start
        # (this path is practically avoided by the check above unless due to rounding)
        return(current_time)
      }
      time_in_interval_to_median <- (log(2) - cumulative_hazard) / lambda_j
      return(current_time + time_in_interval_to_median)
    }
    # Move to next interval
    cumulative_hazard <- cumulative_hazard_end_interval
    current_time <- current_time + delta_t_j
  }
  
  # If we get here, the end-of-grid survival is still > 0.5.
  # Extend beyond last cutpoint with last hazard.
  lambda_last <- tail(hazard_rates, 1)
  if (lambda_last > 0) {
    time_beyond_grid <- (log(2) - cumulative_hazard) / lambda_last
    return(current_time + time_beyond_grid)
  } else {
    # Truly zero tail hazard => median is infinite
    return(Inf)
  }
}


# --- HELPER FUNCTION: Draw Posterior Hazard Samples ---
draw_posterior_hazard_samples <- function(
    num_intervals,
    events_per_interval,
    person_time_per_interval,
    prior_alpha_params,
    prior_beta_params,
    num_samples = 1000
) {
  #' Draws samples from the posterior distribution of hazard rates for each interval
  #' in a Bayesian piecewise exponential model with Gamma priors.
  
  if (!all(sapply(list(events_per_interval, person_time_per_interval, prior_alpha_params, prior_beta_params),
                  length) == num_intervals)) {
    stop("All input vectors must have length equal to num_intervals.")
  }
  
  posterior_hazard_samples <- matrix(NA, nrow = num_samples, ncol = num_intervals)
  
  for (j in 1:num_intervals) {
    posterior_alpha_j <- prior_alpha_params[j] + events_per_interval[j]
    posterior_beta_j <- prior_beta_params[j] + person_time_per_interval[j]
    
    if (posterior_beta_j <= 0) {
      warning(paste("Posterior beta parameter for interval", j, "is non-positive. Setting to a small value (1e-6)."))
      posterior_beta_j <- 1e-6
    }
    
    posterior_hazard_samples[, j] <- rgamma(num_samples, shape = posterior_alpha_j, rate = posterior_beta_j)
  }
  
  return(posterior_hazard_samples)
}

final_vsref_probs_abs <- function(med_arm, med_ref, margin_abs) {
  list(
    p_eff_ref = mean(med_arm >  med_ref + margin_abs),
    p_fut_ref = mean(med_arm <= med_ref - margin_abs)
  )
}

# --- VS-REF MEDIAN SAMPLERS ---------------------------------------------------

sample_vs_ref_medians <- function(slCtrl, slTrt, args, num_samples) {
  num_samples <- num_samples %||% args$num_posterior_draws
  if (isTRUE(args$use_ph_model_vs_ref)) {
    sample_vs_ref_medians_ph(slCtrl, slTrt, args, num_samples)
  } else {
    sample_vs_ref_medians_independent(slCtrl, slTrt, args, num_samples)
  }
}

sample_vs_ref_medians_independent <- function(slCtrl, slTrt, args, num_samples) {
  interval_lengths <- diff(args$interval_cutpoints_sim)
  lamC <- draw_posterior_hazard_samples(
    num_intervals = length(interval_lengths),
    events_per_interval = slCtrl$metrics$events_per_interval,
    person_time_per_interval = slCtrl$metrics$person_time_per_interval,
    prior_alpha_params = args$prior_alpha_params_model,
    prior_beta_params  = args$prior_beta_params_model,
    num_samples = num_samples
  )
  lamT <- draw_posterior_hazard_samples(
    num_intervals = length(interval_lengths),
    events_per_interval = slTrt$metrics$events_per_interval,
    person_time_per_interval = slTrt$metrics$person_time_per_interval,
    prior_alpha_params = args$prior_alpha_params_model,
    prior_beta_params  = args$prior_beta_params_model,
    num_samples = num_samples
  )
  med_ctrl <- apply(lamC, 1, calculate_median_survival_piecewise, interval_lengths = interval_lengths)
  med_trt  <- apply(lamT,  1, calculate_median_survival_piecewise, interval_lengths = interval_lengths)
  list(medCtrl = med_ctrl, medTrt = med_trt)
}

sample_vs_ref_medians_ph <- function(slCtrl, slTrt, args, num_samples) {
  interval_lengths <- diff(args$interval_cutpoints_sim)
  E_C <- slCtrl$metrics$events_per_interval
  PT_C <- slCtrl$metrics$person_time_per_interval
  E_T <- slTrt$metrics$events_per_interval
  PT_T <- slTrt$metrics$person_time_per_interval
  alpha_prior <- args$prior_alpha_params_model
  beta_prior  <- args$prior_beta_params_model
  beta_prior_mean <- coalesce_num(args$ph_loghr_prior_mean, 0)
  beta_prior_sd   <- max(1e-6, coalesce_num(args$ph_loghr_prior_sd, 1))
  
  beta_params <- ph_beta_mode_var(
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

ph_beta_mode_var <- function(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, mu, sigma,
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
