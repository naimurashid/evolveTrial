# Package index

## All functions

- [`CONVERSION_TRIGGERS`](CONVERSION_TRIGGERS.md) : Conversion triggers
- [`FUTILITY_ACTIONS`](FUTILITY_ACTIONS.md) : Futility actions
- [`HYBRID_STATES`](HYBRID_STATES.md) : hybrid_trial.R Core hybrid
  single-arm to between-arm Bayesian adaptive trial simulator
- [`adopt_calibration()`](adopt_calibration.md) : Adopt a calibrated
  design configuration
- [`apply_futility_action()`](apply_futility_action.md) : Apply futility
  action to trial state
- [`apply_recommended_to_args()`](apply_recommended_to_args.md) : Apply
  a recommended early-stopping configuration to the argument list
- [`calculate_median_survival_matrix_cpp()`](calculate_median_survival_matrix_cpp.md)
  : Calculate median survival for multiple hazard samples (C++
  implementation)
- [`calculate_median_survival_piecewise_cpp()`](calculate_median_survival_piecewise_cpp.md)
  : Calculate median survival for piecewise exponential model (C++
  implementation)
- [`calculate_weibull_scale()`](calculate_weibull_scale.md) : Calculate
  Weibull scale from median survival time
- [`calibrate_alpha()`](calibrate_alpha.md) : Calibrate interim and
  final thresholds for single-arm designs
- [`check_conversion_trigger()`](check_conversion_trigger.md) : Check if
  conversion trigger is met
- [`classify_trial_outcome()`](classify_trial_outcome.md) : Classify
  final trial outcome
- [`compare_sa_vs_hybrid()`](compare_sa_vs_hybrid.md) : Compare SA-only
  vs Hybrid decision
- [`compare_simon_to_bo()`](compare_simon_to_bo.md) : Validate BO
  calibration against Simon enumeration
- [`compile_hybrid_results()`](compile_hybrid_results.md) : Compile
  hybrid trial results
- [`compute_ba_posterior()`](compute_ba_posterior.md) : Compute
  between-arm posterior probability
- [`compute_ba_posterior_exponential()`](compute_ba_posterior_exponential.md)
  : Compute exact PP for exponential model (for validation)
- [`compute_ba_posterior_pwe_mc()`](compute_ba_posterior_pwe_mc.md) :
  Monte Carlo comparison for PWE posteriors
- [`compute_interval_metrics()`](compute_interval_metrics.md) : Compute
  interval-specific metrics from registry
- [`compute_p_between_arm()`](compute_p_between_arm.md) : Compute
  between-arm posterior probability
- [`compute_p_single_arm()`](compute_p_single_arm.md) : Compute
  single-arm posterior probability
- [`compute_pp_curve()`](compute_pp_curve.md) : Compute predictive
  probability curve
- [`compute_pp_curve_posterior()`](compute_pp_curve_posterior.md) :
  Compute PP curve using posterior method
- [`compute_pp_curve_predictive()`](compute_pp_curve_predictive.md) :
  Compute PP curve using predictive probability method
- [`compute_pp_posterior()`](compute_pp_posterior.md) : Compute PP using
  posterior method (fast approximation)
- [`compute_pp_posterior_single()`](compute_pp_posterior_single.md) :
  Compute posterior-based PP for single N_add
- [`compute_pp_predictive()`](compute_pp_predictive.md) : Compute PP
  using predictive method (Monte Carlo)
- [`compute_pp_predictive_full()`](compute_pp_predictive_full.md) :
  Compute predictive probability for a single N_add value
- [`compute_pp_predictive_pwe()`](compute_pp_predictive_pwe.md) :
  Compute predictive probability with full PWE model
- [`compute_pwe_median()`](compute_pwe_median.md) : Compute median
  survival from PWE hazards
- [`compute_pwe_median_survival()`](compute_pwe_median_survival.md) :
  Compute median survival from PWE hazards
- [`compute_trial_oc()`](compute_trial_oc.md) : Compute operating
  characteristics for a single trial
- [`count_arm_events()`](count_arm_events.md) : Count events for an arm
  up to analysis time
- [`create_decision_report()`](create_decision_report.md) : Create
  structured decision report
- [`create_hybrid_state()`](create_hybrid_state.md) : Create initial
  hybrid trial state
- [`create_hybrid_theta()`](create_hybrid_theta.md) : Create hybrid
  design parameter structure
- [`create_simulation_cluster()`](create_simulation_cluster.md) : Create
  a reusable evolveTrial simulation cluster
- [`draw_posterior_hazard_samples_cpp()`](draw_posterior_hazard_samples_cpp.md)
  : Draw posterior hazard samples (C++ implementation)
- [`draw_posterior_response_rate()`](draw_posterior_response_rate.md) :
  Draw posterior samples for binary response rate
- [`estimate_vsref_gate_timing()`](estimate_vsref_gate_timing.md) :
  Estimate when the vs-reference interim gates can be satisfied
- [`evaluate_ba_efficacy()`](evaluate_ba_efficacy.md) : Evaluate
  between-arm efficacy
- [`evaluate_ba_futility()`](evaluate_ba_futility.md) : Evaluate
  between-arm futility
- [`evaluate_conversion_trigger()`](evaluate_conversion_trigger.md) :
  Evaluate conversion trigger
- [`evaluate_ph_grid()`](evaluate_ph_grid.md) : Evaluate PH-based grid
  of designs
- [`evaluate_sa_efficacy()`](evaluate_sa_efficacy.md) :
  hybrid_decisions.R Decision rules for hybrid single-arm to between-arm
  trials
- [`evaluate_sa_futility()`](evaluate_sa_futility.md) : Evaluate
  single-arm futility for an arm
- [`exp_arms_from_args()`](exp_arms_from_args.md) : Extract experimental
  arm names from an argument list
- [`expected_information_gain()`](expected_information_gain.md) :
  Expected information gain from additional patients
- [`explore_early_stopping_from_cal()`](explore_early_stopping_from_cal.md)
  : Explore early stopping knobs around a calibrated design
- [`export_scenario_table_to_excel()`](export_scenario_table_to_excel.md)
  : Export a scenario summary table to Excel
- [`export_scenario_table_to_png()`](export_scenario_table_to_png.md) :
  Render the scenario summary table to a PNG image
- [`filter_early_grid()`](filter_early_grid.md) : Filter early-stopping
  designs by operating targets
- [`find_optimal_n_add()`](find_optimal_n_add.md) : Optimal N_add finder
- [`find_simon_design()`](find_simon_design.md) : Find optimal Simon
  design
- [`generate_decision_summary()`](generate_decision_summary.md) :
  Generate human-readable decision summary
- [`get_historical_hazard()`](get_historical_hazard.md) : Get historical
  hazard for an arm
- [`grid_calibrate()`](grid_calibrate.md) : Calibrate historical-control
  thresholds over a grid
- [`handle_sa_futility()`](handle_sa_futility.md) : Handle single-arm
  futility for an arm
- [`handle_state_between()`](handle_state_between.md) : Handle
  STATE_BETWEEN phase
- [`handle_state_consider_conversion()`](handle_state_consider_conversion.md)
  : Handle STATE_CONSIDER_CONVERSION phase
- [`handle_state_single()`](handle_state_single.md) : Handle
  STATE_SINGLE phase
- [`is_trial_success()`](is_trial_success.md) : Determine if trial
  concluded with overall success
- [`make_conversion_decision()`](make_conversion_decision.md) : Make
  conversion decision based on PP curve
- [`perform_ssr()`](perform_ssr.md) : hybrid_ssr.R Sample Size
  Re-Estimation methods for hybrid trials
- [`ph_beta_mode_var()`](ph_beta_mode_var.md) : Compute mode and
  variance of log-HR posterior for PH model (DISPATCHER)
- [`ph_beta_mode_var_cpp()`](ph_beta_mode_var_cpp.md) : Compute mode and
  variance of log-HR posterior for PH model (C++ implementation)
- [`plot_calibration()`](plot_calibration.md) : Plot power versus type I
  error for calibration grids
- [`plot_early_tradeoff()`](plot_early_tradeoff.md) : Plot
  early-stopping trade-offs
- [`pretty_scenario_matrix()`](pretty_scenario_matrix.md) : Summarise
  simulation output by scenario and arm
- [`prob_response_below()`](prob_response_below.md) : Calculate
  posterior probability that response rate is below threshold
- [`prob_response_exceeds()`](prob_response_exceeds.md) : Calculate
  posterior probability that response rate exceeds threshold
- [`recommend_design_from_early()`](recommend_design_from_early.md) :
  Recommend a single early-stopping design
- [`release_cluster()`](release_cluster.md) : Release an evolveTrial
  simulation cluster
- [`run_scenarios()`](run_scenarios.md) : Evaluate a design across
  multiple scenarios
- [`run_simulation_binary()`](run_simulation_binary.md) : Run binary
  endpoint trial simulation
- [`run_simulation_pure()`](run_simulation_pure.md) : Run a full set of
  evolveTrial simulations for a design specification
- [`sample_vs_ref_medians_independent()`](sample_vs_ref_medians_independent.md)
  : Sample medians for vs-ref independent model (DISPATCHER)
- [`sample_vs_ref_medians_independent_cpp()`](sample_vs_ref_medians_independent_cpp.md)
  : Sample medians for vs-ref independent model (C++ implementation)
- [`sample_vs_ref_medians_ph()`](sample_vs_ref_medians_ph.md) : Sample
  medians for vs-ref PH model (DISPATCHER)
- [`sample_vs_ref_medians_ph_cpp()`](sample_vs_ref_medians_ph_cpp.md) :
  Sample medians for vs-ref PH model (C++ implementation)
- [`scenarios_from_grid()`](scenarios_from_grid.md) : Build scenario
  overrides from a grid of design choices
- [`select_arms_for_ba()`](select_arms_for_ba.md) : Select arms for
  between-arm comparison
- [`set_futility_medians()`](set_futility_medians.md) : Adjust futility
  medians for experimental arms
- [`simon_oc_exact()`](simon_oc_exact.md) : Calculate Simon design
  operating characteristics analytically
- [`simulate_future_arm_data()`](simulate_future_arm_data.md) : Simulate
  future arm data under PWE model
- [`simulate_future_arm_pwe()`](simulate_future_arm_pwe.md) : Simulate
  future arm data under PWE model
- [`simulate_pwe_survival()`](simulate_pwe_survival.md) : Simulate
  survival time from PWE model
- [`simulate_pwe_time()`](simulate_pwe_time.md) : Simulate survival time
  from PWE model
- [`test_armadillo_random()`](test_armadillo_random.md) : Test Armadillo
  matrix operations
- [`test_rcpp_sum()`](test_rcpp_sum.md) : Test Rcpp setup
- [`update_hybrid_state()`](update_hybrid_state.md) : Update hybrid
  trial state based on current conditions
- [`update_posteriors()`](update_posteriors.md) : Update posterior
  parameters from trial data
- [`validate_exponential_ba()`](validate_exponential_ba.md) : Validate
  MC against closed-form for exponential
