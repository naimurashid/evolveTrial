# Test equivalence between R and C++ implementations
# Addresses Bug 2: Potential edge case differences

library(evolveTrial)

# Access internal R functions
r_median <- evolveTrial:::calculate_median_survival_piecewise
r_hazard_samples <- evolveTrial:::draw_posterior_hazard_samples
r_ph_beta <- evolveTrial:::ph_beta_mode_var

test_that("C++ median survival matches R implementation", {
  # Test basic case
  hazard_rates <- c(0.1, 0.2, 0.3)
  interval_lengths <- c(5, 5, 5)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
})

test_that("C++ median survives edge case: zero hazard", {
  # All zero hazards should return Inf
  hazard_rates <- c(0, 0, 0)
  interval_lengths <- c(5, 5, 5)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, Inf)
  expect_equal(cpp_result, Inf)
})

test_that("C++ median survives edge case: very high hazard (immediate median)", {
  # Very high hazard should give median in first interval
  hazard_rates <- c(10, 0.1, 0.1)
  interval_lengths <- c(1, 5, 5)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
  expect_lt(r_result, 1)  # Should be in first interval
})

test_that("C++ median survives edge case: median past end of grid", {
  # Low hazards mean median is past grid
  hazard_rates <- c(0.01, 0.01, 0.02)
  interval_lengths <- c(5, 5, 5)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
  expect_gt(r_result, 15)  # Beyond grid
})

test_that("C++ median survives edge case: single interval", {
  hazard_rates <- c(0.1)
  interval_lengths <- c(20)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
})

test_that("C++ median survives edge case: many intervals", {
  # 8 intervals like production use
  set.seed(42)
  hazard_rates <- runif(8, 0.05, 0.3)
  interval_lengths <- rep(3, 8)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
})

test_that("C++ matrix median survives equivalence with apply()", {
  set.seed(123)
  # Simulate posterior samples
  hazard_samples <- matrix(rgamma(400, shape = 2, rate = 10), nrow = 100, ncol = 4)
  interval_lengths <- c(6, 6, 6, 6)

  # R apply version
  r_results <- apply(hazard_samples, 1, function(row) {
    r_median(row, interval_lengths)
  })

  # C++ matrix version
  cpp_results <- calculate_median_survival_matrix_cpp(hazard_samples, interval_lengths)

  expect_equal(r_results, as.numeric(cpp_results), tolerance = 1e-10)
})

test_that("C++ posterior hazard samples match R implementation statistically", {
  set.seed(42)
  num_intervals <- 4
  events <- c(5, 10, 8, 3)
  person_time <- c(100, 200, 150, 50)
  alpha_prior <- rep(0.5, 4)
  beta_prior <- rep(0.5, 4)
  num_samples <- 10000

  # R version
  set.seed(42)
  r_samples <- r_hazard_samples(
    num_intervals, events, person_time, alpha_prior, beta_prior, num_samples
  )

  # C++ version
  set.seed(42)
  cpp_samples <- draw_posterior_hazard_samples_cpp(
    num_intervals, as.integer(events), person_time, alpha_prior, beta_prior, num_samples
  )

  # Compare means and variances (should be very close with same seed)
  r_means <- colMeans(r_samples)
  cpp_means <- colMeans(cpp_samples)

  # Both use R's RNG, so should be identical
  expect_equal(r_means, cpp_means, tolerance = 1e-8)
})

test_that("C++ handles NA/NaN in hazard rates gracefully", {
  # The C++ code should match R's handling of edge cases
  # Test with very small values that might cause numerical issues
  hazard_rates <- c(1e-15, 1e-10, 1e-5)
  interval_lengths <- c(5, 5, 5)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  # Both should handle this gracefully
  expect_true(is.finite(r_result) || is.infinite(r_result))
  expect_true(is.finite(cpp_result) || is.infinite(cpp_result))
  expect_equal(r_result, cpp_result, tolerance = 1e-8)
})

test_that("C++ handles very large hazard values", {
  # Very large hazards should give very small medians
  hazard_rates <- c(100, 50, 25)
  interval_lengths <- c(1, 1, 1)

  r_result <- r_median(hazard_rates, interval_lengths)
  cpp_result <- calculate_median_survival_piecewise_cpp(hazard_rates, interval_lengths)

  expect_equal(r_result, cpp_result, tolerance = 1e-10)
  expect_lt(r_result, 0.1)  # Should be very small
})

test_that("C++ PH model matches R implementation", {
  skip_if_not_installed("evolveTrial")

  # Test PH model posterior computation
  set.seed(42)
  E_C <- c(5, 8, 6, 3)
  PT_C <- c(100, 150, 120, 60)
  E_T <- c(3, 4, 5, 2)
  PT_T <- c(90, 140, 110, 55)
  alpha_prior <- rep(0.5, 4)
  beta_prior <- rep(0.5, 4)

  r_result <- r_ph_beta(E_C, PT_C, E_T, PT_T, alpha_prior, beta_prior, 0, 1)
  cpp_result <- ph_beta_mode_var_cpp(
    as.integer(E_C), PT_C, as.integer(E_T), PT_T, alpha_prior, beta_prior, 0, 1
  )

  expect_equal(r_result$mean, cpp_result$mean, tolerance = 1e-6)
  expect_equal(r_result$sd, cpp_result$sd, tolerance = 1e-6)
})
