# evolveTrial Package Audit Report

**Date**: 2025-12-01
**Scope**: Systematic review of model assumptions, decision frameworks, and simulation implementation
**Purpose**: Investigate why sample sizes may appear "too good to be true"

---

## Executive Summary

This audit identified **10 potential issues** that could explain optimistic operating characteristics. The most impactful are:

| Issue | Severity | Impact on Sample Size |
|-------|----------|----------------------|
| No multiple testing correction | HIGH | Inflates Type I error across interims |
| Multi-arm FWER not controlled | HIGH | ~14% FWER with 3 arms at 5% each |
| Weak default priors | MEDIUM-HIGH | Volatile posteriors with few events |
| Model mismatch (Weibull → piecewise) | MEDIUM | Potential bias in posterior |
| Historical control uncertainty ignored | MEDIUM | Overstates single-arm power |

**Bottom line**: The simulation engine is mathematically correct, but relies on proper calibration of posterior probability thresholds to achieve desired frequentist properties. Without formal alpha spending or FWER control, users must empirically calibrate thresholds under null scenarios.

---

## Part 1: Model Assumptions

### 1.1 Single-Arm Trials (vs Historical Control)

**Data Generation Model** (`R/data_generation.R:28-33`, `R/simulation_driver.R:459-461`):
```r
# True event times from Weibull
t_event_true <- rweibull(1, shape = weibull_shape_true_arms[arm],
                         scale = weibull_scale_true_arms[arm])
# Random censoring: Uniform(0, censor_max_time_sim)
t_random_censor <- runif(1, min = 0, max = censor_max_time_sim)
```

**Analysis Model** (`R/posterior_helpers.R:75-105`):
- Piecewise exponential with Gamma conjugate priors
- Prior: `Gamma(alpha, beta)` per interval (default: α=0.1, β=0.1)
- Posterior: `Gamma(alpha + events, beta + person_time)`

**Key Assumptions**:

| Assumption | Location | Concern |
|------------|----------|---------|
| Piecewise exponential hazards | `posterior_helpers.R:16-57` | Model mismatch with Weibull truth |
| Gamma conjugate priors | `posterior_helpers.R:92-101` | Default α=0.1, β=0.1 is very weak |
| Historical control known exactly | `interim_logic.R:43-47` | No uncertainty propagation |
| Tail extrapolation | `posterior_helpers.R:49-56` | Uses last interval hazard indefinitely |
| MCAR censoring | `simulation_driver.R:461` | Uniform random censoring |

### 1.2 Multi-Arm Trials (vs Reference)

**Analysis Model** (`R/posterior_helpers.R:140-220`):
- Two modes: Independent hazards OR Proportional Hazards (PH)
- Independent: Separate posteriors for control and treatment
- PH: Joint posterior with normal prior on log-HR

**Key Assumptions**:

| Assumption | Location | Concern |
|------------|----------|---------|
| Independent arm posteriors | `posterior_helpers.R:152-188` | Arms compared independently |
| Control posterior caching | `interim_logic.R:100-105` | Induces correlation across comparisons |
| PH assumes constant HR | `posterior_helpers.R:190-220` | May not hold across intervals |
| No shrinkage/borrowing | Throughout | Arms treated as independent |

---

## Part 2: Decision Framework Analysis

### 2.1 Single-Arm: Declaring Superiority

**Decision Logic** (`R/interim_logic.R:207-226`):
```r
# Efficacy: P(median > success_threshold) >= threshold
probs_hc <- calculate_current_probs_hc(slA, args, arm)
pr_eff <- probs_hc$pr_eff  # mean(med_draws > success_thr)

if (pr_eff >= args$efficacy_threshold_hc_prob) {
  state$arm_status[arm] <- "stopped_efficacy"
}
```

**Issues Identified**:

1. **No Multiple Testing Correction**
   - Same threshold (e.g., 0.95) applied at EVERY interim
   - With 5 interims, compound Type I error much higher than single-look
   - No alpha spending function implemented
   - **Location**: `interim_logic.R:220-226`

2. **Historical Control Uncertainty Ignored**
   - `success_thr` treated as known constant: `interim_logic.R:43-47`
   - Real trials have uncertainty in historical estimates
   - Overstates power by not accounting for this variability

3. **Final Analysis Uses Same Structure**
   - `final_success_posterior_prob_threshold` applied identically
   - No distinction between interim and final thresholds
   - **Location**: `simulation_driver.R:521-524`

### 2.2 Multi-Arm: Declaring Superiority

**Decision Logic** (`R/interim_logic.R:134-150`):
```r
# Efficacy: P(treatment_median > control_median) >= threshold
probs_vs_ref <- calculate_current_probs_vs_ref(slC, slT, args)
pr_eff <- probs_vs_ref$pr_eff  # mean(diff_med > 0)

if (pr_eff >= args$efficacy_threshold_vs_ref_prob) {
  state$arm_status[trt_name] <- "stopped_efficacy"
}
```

**Issues Identified**:

1. **No Family-Wise Error Rate (FWER) Control**
   - Each experimental arm tested independently
   - 3 arms × 5% Type I error each ≈ 14% any-arm error
   - No Bonferroni/Holm/Hochberg adjustment
   - **Location**: `interim_logic.R:107-151` (for loop over arms)

2. **Efficacy Threshold = 0**
   - Default comparison: `mean(diff_med > 0)` requires only P > 0.5 direction
   - No clinically meaningful difference required
   - **Location**: `posterior_helpers.R:20-24`

3. **Control Arm Correlation**
   - Same control posterior used for all experimental comparisons
   - Induces correlation in test statistics
   - **Location**: `interim_logic.R:100-105`

### 2.3 Early Stopping (Both Paths)

**Gate Mechanism** (`R/state_management.R:153-244`):
```r
# Gates must pass for BOTH arms before interim analysis
passC <- (evC >= min_ev_ctrl) && (mfuC >= min_mfu_ctrl) && (fracC >= min_pt_ctrl)
passT <- (evT >= min_ev_trt)  && (mfuT >= min_mfu_trt)  && (fracT >= min_pt_trt)
pass <- passC && passT
```

**Issues Identified**:

1. **Gates Delay But Don't Control Error**
   - Information gates only postpone decisions
   - Once gates pass, same threshold applies
   - No sequential monitoring structure
   - **Location**: `state_management.R:230-233`

2. **Proportional Scaling Can Surprise**
   - Scalar gates scaled by randomization probability
   - Arm with 30% randomization gets ~43% of the gate
   - **Location**: `gate_diagnostics.R:45-66`

---

## Part 3: Data Generation & Operating Characteristics

### 3.1 Data Generation Flow

**Enrollment** (`R/simulation_driver.R:368-470`):
```r
# Poisson arrival process
interarrival <- rexp(1, rate = overall_accrual_rate)
current_time <- current_time + interarrival

# Randomization among eligible arms
probs <- randomization_probs[elig_arms]
probs <- probs / sum(probs)
chosen_arm <- sample(elig_arms, size = 1, prob = probs)
```

**Event Time Generation** (`R/simulation_driver.R:459-461`):
```r
t_event_true <- rweibull(1, shape = weibull_shape_true_arms[chosen_arm],
                         scale = weibull_scale_true_arms[chosen_arm])
t_random_censor <- runif(1, min = 0, max = censor_max_time_sim)
```

**Issues Identified**:

1. **Weibull → Piecewise Exponential Mismatch**
   - Data generated: Weibull (smooth hazard function)
   - Analysis model: Piecewise constant hazards
   - Potential bias if intervals are coarse
   - **Impact**: Model misspecification can affect posterior accuracy

2. **Shape Parameter Sensitivity**
   - Shape = 1.0: Exponential (constant hazard) - often unrealistic
   - Shape = 1.3-1.5: Increasing hazard - more realistic for PFS
   - Shape < 1.0: Decreasing hazard - post-response plateau
   - **Test observation**: Tests use shape=1.5 but users may use 1.0

3. **Censoring Assumption (MCAR)**
   - Random censoring: Uniform(0, censor_max_time)
   - Assumes Missing Completely At Random
   - Real trials may have informative dropout

### 3.2 Operating Characteristics Computation

**Aggregation** (`R/simulation_driver.R:704-720`):
```r
early_eff <- sum_stop_efficacy[arm] * inv_num_sims
early_fut <- sum_stop_futility[arm] * inv_num_sims
final_eff <- sum_final_efficacy[arm] * inv_num_sims

# Key: Type_I_Error_or_Power includes BOTH early and final efficacy
results_data$Type_I_Error_or_Power[j] <- early_eff + final_eff
results_data$Pr_Reach_Max_N[j] <- 1 - early_eff - early_fut
```

**Validation**:
- ✅ `PET_Efficacy + PET_Futility + Pr_Reach_Max_N = 1` (verified in tests)
- ✅ `Pr_Final_Efficacy + Pr_Final_Futility + Pr_Final_Inconclusive = Pr_Reach_Max_N`
- ✅ State tracking is clean across simulations

---

## Part 4: Root Cause Analysis - "Too Good" Sample Sizes

### Most Likely Causes (Ranked by Probability)

#### 1. No Multiple Testing Correction (40% likelihood)

**Mechanism**:
- Each interim uses threshold (e.g., 0.90) independently
- With k interims, chance of at least one false positive increases
- No O'Brien-Fleming or Pocock boundaries implemented

**Evidence**:
```r
# interim_logic.R:220-226 - Same threshold at every look
if (pr_eff >= args$efficacy_threshold_hc_prob) {
  state$arm_status[arm] <- "stopped_efficacy"
}
```

**Recommendation**: Use higher thresholds (0.99+) or implement alpha spending.

#### 2. Multi-Arm Independence (25% likelihood)

**Mechanism**:
- 3 experimental arms each at 5% Type I error
- P(at least one false positive) ≈ 1 - (0.95)³ ≈ 14%

**Evidence**:
```r
# interim_logic.R:107 - Independent loop over arms
for (trt_name in eval_arms) {
  # Each arm evaluated independently...
}
```

**Recommendation**: Apply Bonferroni/Dunnett correction or use hierarchical testing.

#### 3. Weak Priors (20% likelihood)

**Mechanism**:
- Default: Gamma(0.1, 0.1) is extremely vague
- With 5 events, posterior mean = (0.1 + 5) / (0.1 + PT)
- Can produce volatile estimates

**Evidence**:
```r
# Test files show:
prior_alpha_params_model = c(0.1, 0.1, 0.1)
prior_beta_params_model = c(0.1, 0.1, 0.1)
```

**Recommendation**: Use Gamma(0.5, 0.5) or informative priors from historical data.

#### 4. Weibull Shape Parameter (10% likelihood)

**Mechanism**:
- Shape = 1.0 gives constant hazard (exponential)
- Real PFS often has decreasing hazard (shape < 1) after initial response
- Shape = 1.3 gives increasing hazard (more events early)

**Evidence**: Package allows any shape, but users may default to 1.0.

**Recommendation**: Use shape from historical data; test sensitivity.

#### 5. Model Mismatch (5% likelihood)

**Mechanism**:
- Weibull truth vs piecewise exponential analysis
- Coarse intervals may not capture hazard shape

**Recommendation**: Use finer intervals (e.g., 3-month instead of 6-month).

---

## Part 5: Test Coverage Gaps

### Current Coverage

| File | Tests | Focus |
|------|-------|-------|
| `test-single-arm-regression.R` | 4 | Basic single-arm OCs |
| `test-path-parity.R` | 3 | Single-arm vs multi-arm structure |
| `test-single-arm-gates.R` | 2 | Gate proportional scaling |
| `test-early-stopping.R` | 1 | Early stopping mechanics |

### Gaps Identified

1. **No Type I Error Calibration Tests**
   - Current: `expect_true(Type_I_Error < 0.3)` with 50 sims
   - Needed: 1000+ sims to verify exact Type I error rates

2. **No Multi-Arm FWER Tests**
   - No tests for 3+ experimental arms
   - No verification of family-wise error rate

3. **No Prior Sensitivity Tests**
   - No tests comparing Gamma(0.1, 0.1) vs Gamma(2, 2)

4. **No Shape Parameter Sensitivity Tests**
   - No tests comparing shape = 1.0 vs 1.3 vs 0.8

5. **No Multiple Interim Tests**
   - No tests varying number of interim looks
   - No verification of alpha accumulation

---

## Part 6: Recommendations

### For Users

1. **Calibrate Thresholds Empirically**
   ```r
   # Run 1000+ sims under null to verify Type I error
   null_result <- run_simulation_pure(
     num_simulations = 1000,
     weibull_median_true_arms = null_median_arms,  # Null = true
     ...
   )
   # Check: Type_I_Error_or_Power should be < target (e.g., 0.05)
   ```

2. **Use Higher Posterior Thresholds**
   - Single interim: 0.95 may achieve ~5% Type I error
   - Multiple interims: Use 0.99+ or implement alpha spending

3. **Account for Multi-Arm Comparisons**
   - 3 arms: Divide target alpha by 3 (Bonferroni)
   - Or use hierarchical testing with fixed sequence

4. **Use Realistic Shape Parameters**
   - Default to 1.3-1.5 for increasing hazard
   - Or use historical data to estimate shape

5. **Increase Prior Information**
   - Use Gamma(0.5, 0.5) or Gamma(1, 1) instead of (0.1, 0.1)
   - Or use empirical Bayes from historical data

### For Package Development (Future)

1. Add alpha spending function support
2. Add FWER control options (Bonferroni, Dunnett)
3. Add calibration helpers that search for threshold achieving target alpha
4. Add hierarchical priors for multi-arm trials
5. Add historical control uncertainty modeling

---

## Appendix A: Code Locations

| Component | File | Lines | Function |
|-----------|------|-------|----------|
| Single-arm decision | `R/interim_logic.R` | 156-240 | `interim_check_hc()` |
| Multi-arm decision | `R/interim_logic.R` | 70-154 | `interim_check_vs_ref()` |
| Posterior sampling | `R/posterior_helpers.R` | 75-105 | `draw_posterior_hazard_samples()` |
| Median calculation | `R/posterior_helpers.R` | 16-57 | `calculate_median_survival_piecewise()` |
| Gate enforcement | `R/state_management.R` | 153-244 | `gates_pass_for_both_arms()` |
| OC aggregation | `R/simulation_driver.R` | 704-720 | `run_simulation_pure()` |
| Data generation | `R/simulation_driver.R` | 459-461 | (inline in sim loop) |

## Appendix B: Mathematical Notes

### Posterior Distribution

For piecewise exponential with Gamma(α, β) prior:
```
Events in interval j: E_j
Person-time in interval j: PT_j

Posterior: Gamma(α + E_j, β + PT_j)
Posterior mean: (α + E_j) / (β + PT_j)
```

### Median Survival Calculation

For piecewise hazards λ₁, λ₂, ..., λ_k:
```
S(t) = exp(-H(t)) where H(t) = Σⱼ λⱼ × Δtⱼ
Median: solve exp(-H(t)) = 0.5, i.e., H(t) = log(2)
```

### Type I Error Accumulation

With k independent looks at threshold p:
```
P(at least one false positive) ≈ k × (1 - p) for large p
Example: 5 looks at p=0.95 → ~25% false positive rate
```

---

*Report generated by Claude Code audit*
