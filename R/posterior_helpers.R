
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

