# Hybrid Single-Arm to Multi-Arm Bayesian Adaptive Trial Design Specification

**Version**: 1.1 **Date**: 2025-12-20 **Status**: IMPLEMENTED - All Core
Components Complete

------------------------------------------------------------------------

## 1. Design Overview

### 1.1 Purpose

This specification defines a hybrid trial design for evolveTrial
that: 1. Begins with **single-arm (SA) monitoring** of each experimental
arm against pre-specified historical controls 2. Upon meeting SA
efficacy criteria, evaluates **conversion** to a between-arm (BA)
comparison 3. Uses **predictive probability** to determine if additional
enrollment for BA comparison is worthwhile 4. If converted, continues
with **seamless** BA monitoring using all accumulated data

### 1.2 Key Design Decisions

| Decision                  | Choice                                             | Rationale                                                 |
|---------------------------|----------------------------------------------------|-----------------------------------------------------------|
| SA efficacy criteria      | Configurable: “any”, “all”, or “k_of_K”            | Flexibility for different trial contexts                  |
| Sample size method        | Both predictive AND posterior                      | Predictive for production, posterior for fast prototyping |
| SA futility handling      | Configurable: “stop_trial”, “drop_arm”, “continue” | Per-arm flexibility                                       |
| Data use after conversion | Seamless (use all data)                            | FDA-accepted, maximizes information                       |
| Survival model            | Piecewise exponential (PWE)                        | Matches evolveTrial, handles non-constant hazards         |

### 1.3 Literature Foundation

- **BASIC Design** (Thall et al. 2023) - Bayesian adaptive
  synthetic-control
- **MAMS Literature** - Multi-arm multi-stage seamless designs
- **SSR Literature** - Predictive probability for sample size
  re-estimation (Chuang-Stein 2007)
- **FDA Guidance** (2019) - Adaptive designs with prospective
  specification

------------------------------------------------------------------------

## 2. Mathematical Framework

### 2.1 Survival Model

**Piecewise Exponential (PWE)** with K intervals: - Hazard rates: λ_k
for arm k, interval j - Gamma conjugate priors: λ\_{k,j} ~ Gamma(α_0,
β_0)

**Posterior Update**:

    α_post = α_prior + n_events
    β_post = β_prior + person_time

### 2.2 Historical Control

Historical control hazard is a **fixed constant** (not modeled as
random):

    λ_hist = log(2) / median_hist

### 2.3 Within-Arm Comparison (vs Historical)

Probability that arm k beats historical by threshold c_k:

    p_k^single = P(λ_k < c_k × λ_hist | data)

For exponential (K=1): **Closed-form Gamma CDF**

``` r
p_single <- pgamma(c_k * lambda_hist, shape = a_post, rate = b_post)
```

For PWE (K\>1): **Posterior sampling** using existing evolveTrial
infrastructure

### 2.4 Between-Arm Comparison

Probability that arm A is superior to arm B:

    p_AB^between = P(λ_A / λ_B < 1 | data)

For exponential: **Closed-form F-distribution**

``` r
p_between <- pf(
  q = (b_A / b_B) * (a_B / a_A),
  df1 = 2 * a_A,
  df2 = 2 * a_B
)
```

For PWE: **Posterior sampling** using
[`sample_vs_ref_medians_ph()`](reference/sample_vs_ref_medians_ph.md)

------------------------------------------------------------------------

## 3. State Machine Architecture

### 3.1 States

    STATE_SINGLE              # Primary SA monitoring
    STATE_CONSIDER_CONVERSION # Evaluating conversion via PP
    STATE_BETWEEN             # Full BA comparison phase
    STATE_STOP                # Terminal state

### 3.2 Transitions

                        ┌─────────────────┐
                        │  STATE_SINGLE   │
                        └────────┬────────┘
                                 │
              ┌──────────────────┼──────────────────┐
              │                  │                  │
              ▼                  ▼                  ▼
      ┌────────────────┐ ┌────────────────┐ ┌─────────────────────┐
      │ All arms       │ │ Max N reached  │ │ SA trigger met      │
      │ futile         │ │ (no trigger)   │ │                     │
      └───────┬────────┘ └───────┬────────┘ └──────────┬──────────┘
              │                  │                     │
              │                  │                     ▼
              │                  │          ┌─────────────────────────┐
              │                  │          │ STATE_CONSIDER_CONVERSION│
              │                  │          └──────────┬──────────────┘
              │                  │                     │
              │                  │         ┌───────────┴───────────┐
              │                  │         │                       │
              │                  │         ▼                       ▼
              │                  │  ┌─────────────┐         ┌──────────┐
              │                  │  │ PP >= pp_go │         │PP < nogo │
              │                  │  └──────┬──────┘         └────┬─────┘
              │                  │         │                     │
              │                  │         ▼                     │
              │                  │  ┌─────────────────┐          │
              │                  │  │  STATE_BETWEEN  │          │
              │                  │  └────────┬────────┘          │
              │                  │           │                   │
              ▼                  ▼           ▼                   ▼
           ┌────────────────────────────────────────────────────────┐
           │                        STATE_STOP                      │
           └────────────────────────────────────────────────────────┘

### 3.3 State Logic

**STATE_SINGLE**: - Enroll patients, compute SA posteriors at each
interim - Check SA efficacy: `p_k^single > eff_sa` - Check SA futility:
`p_k^single < fut_sa` - Apply futility action per arm configuration -
Check transition trigger

**STATE_CONSIDER_CONVERSION**: - Compute PP curve for candidate N_add
values - If `PP(N_add) >= pp_go` for some viable N → proceed - If
`max(PP) < pp_nogo` → stop (not worth continuing) - Else: ambiguous
region (default to no-go)

**STATE_BETWEEN**: - Continue enrollment to selected N_max - Check BA
efficacy: `p_AB^between > eff_ba` - Check BA futility:
`p_AB^between < fut_ba`

**STATE_STOP**: - Report single-arm and/or between-arm conclusions

------------------------------------------------------------------------

## 4. Decision Rules

### 4.1 Single-Arm Rules

| Rule        | Condition                     | Action                    |
|-------------|-------------------------------|---------------------------|
| SA Efficacy | `P(HR < c_k | data) > eff_sa` | Mark arm as SA-successful |
| SA Futility | `P(HR < c_k | data) < fut_sa` | Apply futility_action     |

**Default thresholds**: - `eff_sa = 0.90` - `fut_sa = 0.10` -
`c_k = 0.80` (20% improvement target)

### 4.2 Transition Triggers

| Trigger              | Condition                             |
|----------------------|---------------------------------------|
| `any_single_success` | At least 1 arm has SA efficacy        |
| `all_single_success` | All active arms have SA efficacy      |
| `k_of_K`             | At least k of K arms have SA efficacy |

### 4.3 Conversion Rules

| Decision  | Condition                         | Action                          |
|-----------|-----------------------------------|---------------------------------|
| GO        | `PP(N_add) >= pp_go` for viable N | Proceed to STATE_BETWEEN        |
| NO-GO     | `max(PP) < pp_nogo`               | Stop with SA conclusions        |
| AMBIGUOUS | `pp_nogo <= max(PP) < pp_go`      | Default to NO-GO (configurable) |

**Default thresholds**: - `pp_go = 0.70` - `pp_nogo = 0.20`

### 4.4 Between-Arm Rules

| Rule        | Condition                      | Action              |
|-------------|--------------------------------|---------------------|
| BA Efficacy | `P(HR_AB < 1 | data) > eff_ba` | Declare BA efficacy |
| BA Futility | `P(HR_AB < 1 | data) < fut_ba` | Declare BA futility |

**Default thresholds**: - `eff_ba = 0.975` - `fut_ba = 0.05`

------------------------------------------------------------------------

## 5. Sample Size Re-Estimation

### 5.1 Predictive Probability Method (Recommended)

**Algorithm**:

    For each candidate N_add in {10, 20, ..., max_additional_n}:
      success_count = 0
      For i = 1 to n_outer (default 1000):
        1. Draw λ_true from current posterior
        2. Simulate future patients under λ_true
        3. Update posterior with future data (seamless)
        4. Compute P(BA success | final data)
        5. If P(BA success) >= eff_ba: success_count++
      PP(N_add) = success_count / n_outer

    Find smallest N_add where PP(N_add) >= pp_go

**Properties**: - Accounts for posterior uncertainty -
Bayesian-consistent - More conservative than conditional power

### 5.2 Posterior Method (Fast Approximation)

**Algorithm**:

    1. Extract point estimate θ̂ from current posterior mode
    2. For each candidate N_add:
       Calculate conditional power: CP(N_add | θ̂)
    3. Find smallest N_add where CP >= target_power

**Properties**: - Analytical (milliseconds) - Can be overconfident with
sparse data - Use for rapid prototyping only

### 5.3 Recommendation

Use **predictive probability** for production calibration and final
design selection. Use **posterior method** for quick feasibility checks
during development.

------------------------------------------------------------------------

## 6. Key Parameters

### 6.1 Phase 1 (Single-Arm)

| Parameter    | Symbol  | Range          | Default | Description                              |
|--------------|---------|----------------|---------|------------------------------------------|
| SA efficacy  | eff_sa  | \[0.80, 0.99\] | 0.90    | Posterior prob threshold for SA efficacy |
| SA futility  | fut_sa  | \[0.01, 0.20\] | 0.10    | Posterior prob threshold for SA futility |
| HR threshold | c_k     | \[0.60, 0.90\] | 0.80    | Target HR vs historical                  |
| Min events   | ev_sa   | \[5, 25\]      | 15      | Minimum events before SA interim         |
| Max N        | nmax_sa | \[30, 80\]     | 40      | Maximum N in SA phase per arm            |

### 6.2 Conversion

| Parameter        | Symbol    | Range          | Default | Description                    |
|------------------|-----------|----------------|---------|--------------------------------|
| PP go            | pp_go     | \[0.50, 0.90\] | 0.70    | PP threshold to proceed to BA  |
| PP no-go         | pp_nogo   | \[0.10, 0.40\] | 0.20    | PP threshold to stop           |
| Trigger          | trigger   | categorical    | “any”   | “any”, “all”, or “k_of_K”      |
| Max additional N | max_add_n | \[30, 100\]    | 60      | Max additional patients for BA |

### 6.3 Phase 2 (Between-Arm)

| Parameter   | Symbol  | Range           | Default | Description                              |
|-------------|---------|-----------------|---------|------------------------------------------|
| BA efficacy | eff_ba  | \[0.95, 0.999\] | 0.975   | Posterior prob threshold for BA efficacy |
| BA futility | fut_ba  | \[0.01, 0.10\]  | 0.05    | Posterior prob threshold for BA futility |
| Min events  | ev_ba   | \[10, 35\]      | 15      | Minimum events before BA interim         |
| Max N       | nmax_ba | \[40, 100\]     | 80      | Maximum N in BA phase per arm            |

### 6.4 Structural

| Parameter      | Symbol         | Range          | Default | Description                   |
|----------------|----------------|----------------|---------|-------------------------------|
| PWE intervals  | n_intervals    | \[4, 12\]      | 8       | Number of piecewise intervals |
| Prior strength | prior_strength | \[0.15, 0.60\] | 0.45    | Gamma prior concentration     |

------------------------------------------------------------------------

## 7. Operating Characteristics

### 7.1 Primary Metrics

| Metric        | Definition                        | Constraint |
|---------------|-----------------------------------|------------|
| power         | P(reject H0 in SA or BA) under H1 | \>= 0.80   |
| type1         | P(reject H0 in SA or BA) under H0 | \<= 0.10   |
| type1_between | P(BA efficacy) under null_between | \<= 0.05   |
| EN_total      | Expected total N across phases    | minimize   |

### 7.2 Secondary Metrics

| Metric       | Definition                              | Purpose                         |
|--------------|-----------------------------------------|---------------------------------|
| P_conversion | P(proceed to BA phase)                  | \>= 0.10 (ensure BA evaluation) |
| P_eff_sa     | P(SA efficacy) under H1                 | Monitor SA stopping             |
| P_fut_sa     | P(SA futility) under H0                 | Monitor early stopping          |
| EN_sa        | Expected N in SA phase                  | Cost analysis                   |
| EN_ba_cond   | Expected N in BA phase given conversion | Cost analysis                   |
| ET_total     | Expected trial duration                 | Timeline planning               |

------------------------------------------------------------------------

## 8. Calibration Scenarios

### 8.1 Null Scenarios (Type I Error Control)

**null_global**: Both arms equal to historical

    true_hr_vs_hist = c(A = 1.0, B = 1.0)
    true_hr_between = 1.0

**null_between**: Both beat historical, no BA difference

    true_hr_vs_hist = c(A = 0.8, B = 0.8)
    true_hr_between = 1.0

### 8.2 Alternative Scenarios (Power)

**alt_both_different**: A better than B, both beat historical

    true_hr_vs_hist = c(A = 0.67, B = 0.89)
    true_hr_between = 0.75

**alt_strong_difference**: Strong BA difference

    true_hr_vs_hist = c(A = 0.6, B = 1.0)
    true_hr_between = 0.6

### 8.3 Edge Cases

**one_arm_futile**: One arm futile

    true_hr_vs_hist = c(A = 0.67, B = 1.2)
    true_hr_between = 0.55

------------------------------------------------------------------------

## 9. BATON Integration

### 9.1 Objective Function

    objective = EN_total + conversion_penalty_weight × P_conversion_under_null_between

**Rationale**: Penalize unnecessary conversions when SA succeeds but BA
shouldn’t.

### 9.2 Constraints

``` r
constraints <- list(
  power = c("ge", 0.80),
  type1 = c("le", 0.10),
  type1_between = c("le", 0.05),
  P_conversion = c("ge", 0.10),
  P_eff_sa = c("le", 0.40)
)
```

### 9.3 Warmstart Budget

| Stage   | Budget | Purpose                                  |
|---------|--------|------------------------------------------|
| Stage 1 | 250    | Exploration (14D)                        |
| Stage 2 | 200    | Refinement (10D after fixing structural) |
| Stage 3 | 150    | Validation (7-8D)                        |

------------------------------------------------------------------------

## 10. Implementation Files

### 10.1 evolveTrial Package

| File                   | Purpose                             |
|------------------------|-------------------------------------|
| `R/hybrid_trial.R`     | Core simulator with 4-state machine |
| `R/hybrid_ssr.R`       | Predictive/posterior SSR methods    |
| `R/hybrid_decisions.R` | SA/BA decision rules, trigger logic |

### 10.2 BATON Integration

| File                               | Purpose                          |
|------------------------------------|----------------------------------|
| `R/warmstart_wrappers.R`           | `create_hybrid_wrapper_simple()` |
| `R/scenario_configs.R`             | Hybrid bounds, constraints       |
| `R/aggregation.R`                  | Hybrid metrics aggregation       |
| `R/hybrid_calibration_scenarios.R` | 5 calibration scenarios          |

### 10.3 Testing

| File                                       | Purpose                                 |
|--------------------------------------------|-----------------------------------------|
| `tests/testthat/test-hybrid-simulator.R`   | Parameter validation, state transitions |
| `tests/testthat/test-hybrid-ssr.R`         | PP algorithms                           |
| `tests/testthat/test-hybrid-integration.R` | Full BO calibration                     |

------------------------------------------------------------------------

## 11. Validation Strategy

### 11.1 Unit Tests

- Closed-form F-distribution matches MC for exponential case
- State transitions occur correctly
- PP increases with sample size
- PP approaches 1 when true effect is large
- PP approaches 0 when true effect is null

### 11.2 Regression Tests

- Hybrid with `pp_go=1.0` matches pure SA design
- Hybrid with instant conversion matches pure BA design
- Type I error controlled across scenarios

### 11.3 Integration Tests

- Full BO calibration completes without errors
- CSV loading and wrapper creation work correctly
- Parallel execution produces consistent results

------------------------------------------------------------------------

## 12. Recommended Defaults for ADAPT

For the ARPA-H ADAPT breast cancer platform trial:

1.  **SSR method**: Predictive probability
2.  **Conversion trigger**: `any_single_success`
3.  **Thresholds**: `pp_go=0.70`, `pp_nogo=0.20`
4.  **Futility action**: `drop_arm`
5.  **P_conversion constraint**: \>= 0.10
6.  **HR threshold**: 0.80 (20% improvement)

------------------------------------------------------------------------

## 13. Implementation Status

### 13.1 Completed Components

| Component             | File                                      | Status   | Tests                  |
|-----------------------|-------------------------------------------|----------|------------------------|
| Core simulator        | `R/hybrid_trial.R`                        | COMPLETE | 60 passing             |
| SSR methods           | `R/hybrid_ssr.R`                          | COMPLETE | Included in unit tests |
| Decision rules        | `R/hybrid_decisions.R`                    | COMPLETE | Included in unit tests |
| BATON wrapper         | `R/warmstart_wrappers.R`                  | COMPLETE | Integration test       |
| Bounds/constraints    | `R/scenario_configs.R`                    | COMPLETE | \-                     |
| Aggregation           | `R/aggregation.R`                         | COMPLETE | \-                     |
| Calibration scenarios | `R/hybrid_calibration_scenarios.R`        | COMPLETE | \-                     |
| Unit tests            | `tests/testthat/test-hybrid-simulator.R`  | COMPLETE | 60 tests               |
| Integration test      | `scripts/test_hybrid_integration.R`       | COMPLETE | 9 test sections        |
| Production CSV        | `scenarios/cohortB_hybrid_production.csv` | COMPLETE | 11 scenarios           |

### 13.2 Key Functions Implemented

**evolveTrial Package:** - `create_hybrid_theta()` - Parameter structure
with validation - `create_hybrid_state()` - State initialization with
registries - `update_hybrid_state()` - Main state machine driver -
`handle_state_single()` - SA phase logic -
`handle_state_consider_conversion()` - PP evaluation -
`handle_state_between()` - BA phase logic -
`compute_pp_curve_predictive()` - Full Monte Carlo PP -
`compute_pp_curve_posterior()` - Fast analytical PP -
`evaluate_sa_efficacy()` - SA efficacy decision -
`evaluate_sa_futility()` - SA futility decision -
`evaluate_conversion_trigger()` - Trigger evaluation -
`evaluate_ba_efficacy()` - BA efficacy decision -
`apply_futility_action()` - Futility handling -
`validate_exponential_ba()` - Closed-form validation -
`compile_hybrid_results()` - Result compilation

**BATON Integration:** - `create_hybrid_wrapper_simple()` - CSV-based
wrapper factory - `get_hybrid_bounds()` - Parameter bounds by stage -
`get_hybrid_constraints()` - OC constraints -
`get_hybrid_calibration_scenarios()` - 5 calibration scenarios -
`aggregate_hybrid_results()` - Multi-scenario aggregation -
`compute_hybrid_objective()` - Objective with conversion penalty

### 13.3 Usage

**Run unit tests:**

``` bash
cd evolveTrial
Rscript -e "testthat::test_file('tests/testthat/test-hybrid-simulator.R')"
```

**Run integration test:**

``` bash
cd adaptive-trial-bo-paper
Rscript scripts/test_hybrid_integration.R
```

**Run production calibration:**

``` bash
cd adaptive-trial-bo-paper
Rscript scripts/run_warmstart_scenarios.R \
  scenarios/cohortB_hybrid_production.csv \
  --filter-ids cohortB_hybrid_optimal
```

**Dry run (validate scenarios):**

``` bash
Rscript scripts/run_warmstart_scenarios.R \
  scenarios/cohortB_hybrid_production.csv \
  --dry-run
```

### 13.4 Production CSV Scenarios

| Scenario ID                 | Description                   | Conversion Trigger |
|-----------------------------|-------------------------------|--------------------|
| cohortB_hybrid_optimal      | Minimize EN under alternative | any_single_success |
| cohortB_hybrid_minimax      | Balance null/alt EN           | any_single_success |
| cohortB_hybrid_fleming      | Balanced EN optimization      | any_single_success |
| cohortB_hybrid_all_success  | Require all arms succeed      | all_single_success |
| cohortB_hybrid_k_of_K       | k of K arms succeed           | k_of_K             |
| cohortB_hybrid_posterior    | Fast posterior SSR            | any_single_success |
| cohortB_hybrid_stop_on_fut  | Stop on any futility          | any_single_success |
| cohortB_hybrid_continue_fut | Continue despite futility     | any_single_success |
| cohortB_hybrid_strong       | Strong effect (HR=0.60)       | any_single_success |
| cohortB_hybrid_moderate     | Moderate effect (HR=0.82)     | any_single_success |
| cohortB_hybrid_weibull13    | Weibull shape=1.3             | any_single_success |

------------------------------------------------------------------------

**Document maintained by**: evolveTrial development team **Last
updated**: 2025-12-20
