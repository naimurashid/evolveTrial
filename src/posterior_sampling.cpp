// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declarations
double calculate_median_survival_piecewise_impl(const vec& hazard_rates, const vec& interval_lengths);

// =============================================================================
// FUNCTION 1: Draw Posterior Hazard Samples
// =============================================================================

//' Draw posterior hazard samples (C++ implementation)
//'
//' Draws samples from the posterior distribution of hazard rates for each interval
//' in a Bayesian piecewise exponential model with Gamma priors.
//'
//' @param num_intervals Number of intervals in piecewise model
//' @param events_per_interval Integer vector of observed events per interval
//' @param person_time_per_interval Numeric vector of person-time at risk per interval
//' @param prior_alpha_params Numeric vector of Gamma prior shape parameters
//' @param prior_beta_params Numeric vector of Gamma prior rate parameters
//' @param num_samples Number of posterior samples to draw
//' @return Matrix with num_samples rows and num_intervals columns of hazard samples
//' @export
// [[Rcpp::export]]
arma::mat draw_posterior_hazard_samples_cpp(
    int num_intervals,
    const arma::ivec& events_per_interval,
    const arma::vec& person_time_per_interval,
    const arma::vec& prior_alpha_params,
    const arma::vec& prior_beta_params,
    int num_samples
) {
  // Validate inputs
  if (events_per_interval.n_elem != num_intervals ||
      person_time_per_interval.n_elem != num_intervals ||
      prior_alpha_params.n_elem != num_intervals ||
      prior_beta_params.n_elem != num_intervals) {
    stop("All input vectors must have length equal to num_intervals");
  }

  // Pre-allocate output matrix (num_samples × num_intervals)
  arma::mat posterior_hazard_samples(num_samples, num_intervals);

  // Loop over intervals and samples
  // Use R's RNG for reproducibility with same seeds
  for (int j = 0; j < num_intervals; j++) {
    double posterior_alpha_j = prior_alpha_params(j) + events_per_interval(j);
    double posterior_beta_j = prior_beta_params(j) + person_time_per_interval(j);

    // Guard against non-positive beta
    if (posterior_beta_j <= 0) {
      Rcpp::warning("Posterior beta parameter for interval " + std::to_string(j) +
                    " is non-positive. Setting to 1e-6.");
      posterior_beta_j = 1e-6;
    }

    // Draw num_samples Gamma samples for this interval using R's RNG
    // R::rgamma uses shape/scale parameterization
    double scale_j = 1.0 / posterior_beta_j;
    for (int i = 0; i < num_samples; i++) {
      posterior_hazard_samples(i, j) = R::rgamma(posterior_alpha_j, scale_j);
    }
  }

  return posterior_hazard_samples;
}

// =============================================================================
// FUNCTION 2: Calculate Median Survival (Single Vector)
// =============================================================================

//' Calculate median survival for piecewise exponential model (C++ implementation)
//'
//' Computes the median survival time for a piecewise exponential model.
//' If the 0.5 survival is not reached by the end of the last interval,
//' continues past the last cutpoint with the last interval's hazard (open-ended tail).
//' Only returns Inf if the last hazard is exactly zero.
//'
//' @param hazard_rates Numeric vector of hazard rates for each interval
//' @param interval_lengths Numeric vector of interval lengths (durations)
//' @return Median survival time (numeric scalar, possibly Inf)
//' @export
// [[Rcpp::export]]
double calculate_median_survival_piecewise_cpp(
    const arma::vec& hazard_rates,
    const arma::vec& interval_lengths
) {
  if (hazard_rates.n_elem != interval_lengths.n_elem) {
    stop("hazard_rates and interval_lengths must have the same length");
  }

  return calculate_median_survival_piecewise_impl(hazard_rates, interval_lengths);
}

// Internal implementation (for reuse in other functions)
double calculate_median_survival_piecewise_impl(
    const arma::vec& hazard_rates,
    const arma::vec& interval_lengths
) {
  const double LOG_2 = 0.69314718055994530942;  // log(2) precomputed
  double cumulative_hazard = 0.0;
  double current_time = 0.0;
  int num_intervals = hazard_rates.n_elem;

  // Search for interval containing median
  for (int i = 0; i < num_intervals; i++) {
    double lambda_j = hazard_rates(i);
    double delta_t_j = interval_lengths(i);
    double cumulative_hazard_end = cumulative_hazard + lambda_j * delta_t_j;

    if (cumulative_hazard_end >= LOG_2) {
      // Median occurs within this interval
      if (lambda_j == 0.0) {
        return current_time;  // Edge case: zero hazard
      }
      double time_in_interval = (LOG_2 - cumulative_hazard) / lambda_j;
      return current_time + time_in_interval;
    }

    // Move to next interval
    cumulative_hazard = cumulative_hazard_end;
    current_time += delta_t_j;
  }

  // Median not reached by end of grid; extend with last hazard
  double lambda_last = hazard_rates(num_intervals - 1);
  if (lambda_last > 0.0) {
    double time_beyond_grid = (LOG_2 - cumulative_hazard) / lambda_last;
    return current_time + time_beyond_grid;
  } else {
    return arma::datum::inf;  // Infinite survival if last hazard is zero
  }
}

// =============================================================================
// FUNCTION 3: Calculate Median Survival (Entire Matrix - Vectorized)
// =============================================================================

//' Calculate median survival for multiple hazard samples (C++ implementation)
//'
//' Vectorized version that processes an entire matrix of posterior hazard samples.
//' Each row is a posterior sample; each column is an interval.
//'
//' @param hazard_samples Matrix with num_samples rows and num_intervals columns
//' @param interval_lengths Numeric vector of interval lengths (durations)
//' @return Numeric vector of median survival times (length num_samples)
//' @export
// [[Rcpp::export]]
arma::vec calculate_median_survival_matrix_cpp(
    const arma::mat& hazard_samples,
    const arma::vec& interval_lengths
) {
  int num_samples = hazard_samples.n_rows;
  int num_intervals = hazard_samples.n_cols;

  if (interval_lengths.n_elem != num_intervals) {
    stop("interval_lengths must have length equal to ncol(hazard_samples)");
  }

  arma::vec medians(num_samples);

  // Process each sample (row)
  for (int i = 0; i < num_samples; i++) {
    arma::vec hazards = hazard_samples.row(i).t();  // Extract row as column vector
    medians(i) = calculate_median_survival_piecewise_impl(hazards, interval_lengths);
  }

  return medians;
}

// =============================================================================
// FUNCTION 4: PH Model - Log-HR Posterior Mode and Variance
// =============================================================================

//' Compute mode and variance of log-HR posterior for PH model (C++ implementation)
//'
//' Uses Newton-Raphson iteration to find the mode of the log-HR posterior
//' under a proportional hazards model with Gamma prior on baseline hazards
//' and normal prior on log-HR.
//'
//' @param E_C Integer vector of control events per interval
//' @param PT_C Numeric vector of control person-time per interval
//' @param E_T Integer vector of treatment events per interval
//' @param PT_T Numeric vector of treatment person-time per interval
//' @param alpha_prior Numeric vector of Gamma prior shape parameters
//' @param beta_prior Numeric vector of Gamma prior rate parameters
//' @param mu Log-HR prior mean
//' @param sigma Log-HR prior SD
//' @param tol Convergence tolerance (default 1e-6)
//' @param max_iter Maximum Newton-Raphson iterations (default 50)
//' @return List with "mean" and "sd" of log-HR posterior
//' @export
// [[Rcpp::export]]
List ph_beta_mode_var_cpp(
    const arma::ivec& E_C,
    const arma::vec& PT_C,
    const arma::ivec& E_T,
    const arma::vec& PT_T,
    const arma::vec& alpha_prior,
    const arma::vec& beta_prior,
    double mu,
    double sigma,
    double tol = 1e-6,
    int max_iter = 50
) {
  // Initialize beta using marginal counts
  double total_events_ctrl = sum(E_C);
  double total_events_trt = sum(E_T);
  double total_pt_ctrl = sum(PT_C);
  double total_pt_trt = sum(PT_T);

  double beta = std::log((total_events_trt + 0.5) / (total_pt_trt + 0.5)) -
                std::log((total_events_ctrl + 0.5) / (total_pt_ctrl + 0.5));

  double sigma2 = sigma * sigma;

  // Newton-Raphson iteration
  for (int iter = 0; iter < max_iter; iter++) {
    double exp_beta = std::exp(beta);
    arma::vec denom = beta_prior + PT_C + exp_beta * PT_T;
    arma::vec g = exp_beta * PT_T / denom;

    double grad = sum(E_T) - sum((alpha_prior + conv_to<vec>::from(E_C + E_T)) % g) - (beta - mu) / sigma2;
    double hess = -sum((alpha_prior + conv_to<vec>::from(E_C + E_T)) % (g - g % g)) - 1.0 / sigma2;

    double step = grad / hess;
    double beta_new = beta - step;

    // Check for numerical issues
    if (!std::isfinite(beta_new)) {
      break;
    }

    beta = beta_new;

    // Check convergence
    if (std::abs(step) < tol) {
      break;
    }
  }

  // Calculate variance
  double exp_beta = std::exp(beta);
  arma::vec denom = beta_prior + PT_C + exp_beta * PT_T;
  arma::vec g = exp_beta * PT_T / denom;
  double hess = -sum((alpha_prior + conv_to<vec>::from(E_C + E_T)) % (g - g % g)) - 1.0 / sigma2;

  double var_beta;
  if (hess < 0) {
    var_beta = std::min(1e6, std::max(1e-6, -1.0 / hess));
  } else {
    var_beta = 1.0;
  }

  return List::create(
    Named("mean") = beta,
    Named("sd") = std::sqrt(var_beta)
  );
}

// =============================================================================
// FUNCTION 5: Sample Medians for PH Model (End-to-End)
// =============================================================================

//' Sample medians for vs-ref PH model (C++ implementation)
//'
//' End-to-end posterior sampler for proportional hazards model.
//' Samples log-HR from posterior, then baseline hazards, then computes medians.
//'
//' @param E_C Integer vector of control events per interval
//' @param PT_C Numeric vector of control person-time per interval
//' @param E_T Integer vector of treatment events per interval
//' @param PT_T Numeric vector of treatment person-time per interval
//' @param alpha_prior Numeric vector of Gamma prior shape parameters
//' @param beta_prior Numeric vector of Gamma prior rate parameters
//' @param interval_lengths Numeric vector of interval durations
//' @param mu Log-HR prior mean
//' @param sigma Log-HR prior SD
//' @param num_samples Number of posterior samples
//' @return List with "medCtrl", "medTrt", and "logHR" vectors
//' @export
// [[Rcpp::export]]
List sample_vs_ref_medians_ph_cpp(
    const arma::ivec& E_C,
    const arma::vec& PT_C,
    const arma::ivec& E_T,
    const arma::vec& PT_T,
    const arma::vec& alpha_prior,
    const arma::vec& beta_prior,
    const arma::vec& interval_lengths,
    double mu,
    double sigma,
    int num_samples
) {
  // Get log-HR posterior parameters
  List beta_params = ph_beta_mode_var_cpp(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, mu, sigma);
  double beta_mean = beta_params["mean"];
  double beta_sd = beta_params["sd"];

  // Draw log-HR samples
  arma::vec beta_draws = randn<vec>(num_samples) * beta_sd + beta_mean;

  // Pre-compute shapes (constant across samples)
  arma::vec shapes = alpha_prior + conv_to<vec>::from(E_C + E_T);

  // Pre-allocate output vectors
  arma::vec med_ctrl(num_samples);
  arma::vec med_trt(num_samples);

  // Sample medians
  for (int i = 0; i < num_samples; i++) {
    double exp_beta = std::exp(beta_draws(i));
    arma::vec rates = beta_prior + PT_C + exp_beta * PT_T;

    // Draw baseline hazards
    arma::vec lambda(shapes.n_elem);
    for (int j = 0; j < shapes.n_elem; j++) {
      lambda(j) = R::rgamma(shapes(j), 1.0 / rates(j));  // shape, scale parameterization
    }

    // Calculate medians
    med_ctrl(i) = calculate_median_survival_piecewise_impl(lambda, interval_lengths);
    med_trt(i) = calculate_median_survival_piecewise_impl(lambda * exp_beta, interval_lengths);
  }

  return List::create(
    Named("medCtrl") = med_ctrl,
    Named("medTrt") = med_trt,
    Named("logHR") = beta_draws
  );
}

// =============================================================================
// FUNCTION 6: Sample Medians for Independent Model (End-to-End)
// =============================================================================

//' Sample medians for vs-ref independent model (C++ implementation)
//'
//' End-to-end posterior sampler for independent hazards model.
//' Used for single-arm trials and multi-arm trials with independent hazards.
//'
//' @param E_C Integer vector of control events per interval
//' @param PT_C Numeric vector of control person-time per interval
//' @param E_T Integer vector of treatment events per interval
//' @param PT_T Numeric vector of treatment person-time per interval
//' @param alpha_prior Numeric vector of Gamma prior shape parameters
//' @param beta_prior Numeric vector of Gamma prior rate parameters
//' @param interval_lengths Numeric vector of interval durations
//' @param num_samples Number of posterior samples
//' @return List with "medCtrl" and "medTrt" vectors (no logHR for independent model)
//' @export
// [[Rcpp::export]]
List sample_vs_ref_medians_independent_cpp(
    const arma::ivec& E_C,
    const arma::vec& PT_C,
    const arma::ivec& E_T,
    const arma::vec& PT_T,
    const arma::vec& alpha_prior,
    const arma::vec& beta_prior,
    const arma::vec& interval_lengths,
    int num_samples
) {
  int num_intervals = interval_lengths.n_elem;

  // Sample control arm hazards
  arma::mat lamC = draw_posterior_hazard_samples_cpp(
    num_intervals, E_C, PT_C, alpha_prior, beta_prior, num_samples
  );

  // Sample treatment arm hazards (independently)
  arma::mat lamT = draw_posterior_hazard_samples_cpp(
    num_intervals, E_T, PT_T, alpha_prior, beta_prior, num_samples
  );

  // Calculate medians for both arms
  arma::vec med_ctrl = calculate_median_survival_matrix_cpp(lamC, interval_lengths);
  arma::vec med_trt = calculate_median_survival_matrix_cpp(lamT, interval_lengths);

  return List::create(
    Named("medCtrl") = med_ctrl,
    Named("medTrt") = med_trt,
    Named("logHR") = R_NilValue  // NULL for independent model
  );
}
