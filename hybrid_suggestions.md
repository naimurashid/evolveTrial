# Hybrid Single-to-Between Arm Bayesian Adaptive Trial Design: Complete Specification

## Table of Contents

1. [Design Overview](#1-design-overview)
2. [Mathematical Framework](#2-mathematical-framework)
3. [Decision Rules](#3-decision-rules)
4. [State Machine Architecture](#4-state-machine-architecture)
5. [Predictive Probability Engine](#5-predictive-probability-engine)
6. [Software Architecture](#6-software-architecture)
7. [BATON Calibration](#7-baton-calibration)
8. [Implementation Phases](#8-implementation-phases)
9. [Testing Strategy](#9-testing-strategy)

---

## 1. Design Overview

### 1.1 Motivation

In resource-constrained settings, investigators often prefer within-arm comparisons to historical controls because they require fewer samples and complete more quickly than between-arm comparisons. However, between-arm comparisons provide stronger evidence. This hybrid design allows:

1. **Initial phase**: Monitor each arm against its historical benchmark (single-arm view) while simultaneously tracking between-arm contrasts
2. **Transition decision**: Once single-arm efficacy is established, use accumulated data to determine whether extending to a powered between-arm comparison is worthwhile
3. **Extension phase**: If warranted, continue enrollment to complete the between-arm comparison

### 1.2 Key Design Features

| Feature | Description |
|---------|-------------|
| **Randomization** | Always randomize among experimental arms from trial start |
| **Dual monitoring** | Compute both within-arm (vs historical) and between-arm posteriors at every interim |
| **Historical controls** | Fixed, pre-specified benchmarks (not modeled as random) |
| **Transition trigger** | Configurable: any arm succeeds, all arms succeed, or k-of-K |
| **Conversion decision** | Based on predictive probability of between-arm success |
| **Sample size adaptation** | Determine additional N needed for between-arm comparison at transition |

### 1.3 Trial Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                              TRIAL FLOW                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌─────────────┐                                                             │
│  │ Trial Start │                                                             │
│  └──────┬──────┘                                                             │
│         │                                                                    │
│         ▼                                                                    │
│  ┌─────────────────────────────────────────────────────────────────────┐    │
│  │                     STATE: SINGLE-ARM PHASE                          │    │
│  │                                                                       │    │
│  │  • Randomize patients to arms A, B, ... (equal or RAR)               │    │
│  │  • At each interim, compute:                                          │    │
│  │      - P(HR_k^hist < c_k | data) for each arm k                      │    │
│  │      - P(HR_AB < 1 | data) [passive monitoring]                      │    │
│  │  • Apply single-arm efficacy/futility rules                          │    │
│  │                                                                       │    │
│  └─────────────────────────┬───────────────────────────────────────────┘    │
│                            │                                                 │
│         ┌──────────────────┼──────────────────┐                             │
│         │                  │                  │                             │
│         ▼                  ▼                  ▼                             │
│  ┌─────────────┐   ┌─────────────┐   ┌─────────────────────────────────┐   │
│  │All arms fail│   │ Max N with  │   │  Single-arm efficacy trigger    │   │
│  │  futility   │   │  no trigger │   │  (e.g., any arm succeeds)       │   │
│  └──────┬──────┘   └──────┬──────┘   └───────────────┬─────────────────┘   │
│         │                 │                          │                      │
│         ▼                 ▼                          ▼                      │
│  ┌─────────────────────────────┐   ┌─────────────────────────────────────┐ │
│  │    STATE: STOP (futility)   │   │   STATE: CONSIDER CONVERSION        │ │
│  │    No efficacy claims       │   │                                      │ │
│  └─────────────────────────────┘   │   • Compute π_pred(N_add) for       │ │
│                                     │     candidate additional sample     │ │
│                                     │     sizes                           │ │
│                                     │   • Determine go/no-go              │ │
│                                     └───────────────┬─────────────────────┘ │
│                                                     │                       │
│                              ┌──────────────────────┴────────────────┐      │
│                              │                                       │      │
│                              ▼                                       ▼      │
│               ┌───────────────────────────┐      ┌────────────────────────┐│
│               │  π_pred < pp_nogo         │      │  π_pred ≥ pp_go        ││
│               │  (not worth continuing)   │      │  (worth continuing)    ││
│               └─────────────┬─────────────┘      └───────────┬────────────┘│
│                             │                                │             │
│                             ▼                                ▼             │
│               ┌───────────────────────────┐    ┌─────────────────────────┐ │
│               │  STATE: STOP              │    │  STATE: BETWEEN-ARM     │ │
│               │  Report single-arm        │    │  PHASE                  │ │
│               │  conclusions only         │    │                         │ │
│               └───────────────────────────┘    │  • Continue enrollment  │ │
│                                                │    to selected N_max    │ │
│                                                │  • Monitor between-arm  │ │
│                                                │    efficacy/futility    │ │
│                                                └────────────┬────────────┘ │
│                                                             │              │
│                                                             ▼              │
│                                                ┌─────────────────────────┐ │
│                                                │  STATE: STOP            │ │
│                                                │  Report single-arm AND  │ │
│                                                │  between-arm conclusions│ │
│                                                └─────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Mathematical Framework

### 2.1 Notation

| Symbol | Description |
|--------|-------------|
| $k = 1, \ldots, K$ | Arm index ($K \geq 2$, all experimental) |
| $\lambda_k$ | Hazard rate for arm $k$ (random, to be estimated) |
| $\lambda_{\text{hist},k}$ | Historical control hazard for arm $k$ (fixed, pre-specified) |
| $n_k$ | Number of events observed in arm $k$ |
| $T_k$ | Total exposure time (person-time) in arm $k$ |
| $c_k$ | Target hazard ratio threshold vs historical (e.g., 0.8) |
| $\text{HR}_k^{(\text{hist})}$ | Hazard ratio of arm $k$ vs historical: $\lambda_k / \lambda_{\text{hist},k}$ |
| $\text{HR}_{jl}$ | Hazard ratio between arms: $\lambda_j / \lambda_l$ |

### 2.2 Survival Model: Exponential with Gamma Prior

**Assumption**: Survival times follow an exponential distribution with arm-specific hazard rates.

$$T_{ik} \mid \lambda_k \sim \text{Exponential}(\lambda_k)$$

where $T_{ik}$ is the survival time for patient $i$ in arm $k$.

**Prior**: Gamma prior on each hazard rate (conjugate):

$$\lambda_k \sim \text{Gamma}(a_0, b_0)$$

where $a_0, b_0$ are typically small (e.g., $a_0 = b_0 = 0.001$) for a vague prior.

**Likelihood**: For $n_k$ events and total exposure $T_k$:

$$L(\lambda_k \mid \text{data}) \propto \lambda_k^{n_k} \exp(-\lambda_k T_k)$$

**Posterior**: By conjugacy:

$$\lambda_k \mid \text{data} \sim \text{Gamma}(a_k, b_k)$$

where:
- $a_k = a_0 + n_k$ (shape parameter)
- $b_k = b_0 + T_k$ (rate parameter)

**Sufficient statistics**: The posterior is fully characterized by $(a_k, b_k)$.

### 2.3 Historical Control Specification

The historical control hazard $\lambda_{\text{hist},k}$ is a **fixed constant**, typically derived from:

$$\lambda_{\text{hist},k} = \frac{\log(2)}{m_{\text{hist},k}}$$

where $m_{\text{hist},k}$ is the historical median survival time for the population relevant to arm $k$.

**Key point**: We do NOT place a prior on $\lambda_{\text{hist},k}$. It enters calculations only as a fixed benchmark.

### 2.4 Within-Arm Comparison (vs Fixed Historical)

**Hazard ratio**:

$$\text{HR}_k^{(\text{hist})} = \frac{\lambda_k}{\lambda_{\text{hist},k}}$$

**Decision quantity** — probability that HR is below target threshold $c_k$:

$$p_k^{\text{single}} = P\left(\text{HR}_k^{(\text{hist})} < c_k \mid \text{data}\right) = P\left(\lambda_k < c_k \cdot \lambda_{\text{hist},k} \mid \text{data}\right)$$

**Closed-form computation**:

Since $\lambda_k \mid \text{data} \sim \text{Gamma}(a_k, b_k)$:

$$p_k^{\text{single}} = F_{\text{Gamma}}\left(c_k \cdot \lambda_{\text{hist},k}; a_k, b_k\right)$$

where $F_{\text{Gamma}}(x; a, b)$ is the CDF of the Gamma distribution with shape $a$ and rate $b$.

**In R**:
```r
p_single <- pgamma(c_k * lambda_hist_k, shape = a_k, rate = b_k)
```

### 2.5 Between-Arm Comparison

**Hazard ratio** between arms $j$ and $l$:

$$\text{HR}_{jl} = \frac{\lambda_j}{\lambda_l}$$

**Decision quantity**:

$$p_{jl}^{\text{between}} = P\left(\text{HR}_{jl} < 1 \mid \text{data}\right) = P\left(\frac{\lambda_j}{\lambda_l} < 1 \mid \text{data}\right)$$

**Distribution of the ratio**:

If $\lambda_j \sim \text{Gamma}(a_j, b_j)$ and $\lambda_l \sim \text{Gamma}(a_l, b_l)$ independently, then:

$$\frac{\lambda_j / b_j}{\lambda_l / b_l} \sim \text{BetaPrime}(a_j, a_l)$$

Or equivalently, using the F-distribution relationship:

$$\frac{a_l}{a_j} \cdot \frac{b_j}{b_l} \cdot \frac{\lambda_j}{\lambda_l} \sim F(2a_j, 2a_l)$$

**Closed-form computation**:

$$P\left(\frac{\lambda_j}{\lambda_l} < c\right) = F_F\left(c \cdot \frac{b_j}{b_l} \cdot \frac{a_l}{a_j}; 2a_j, 2a_l\right)$$

where $F_F(x; d_1, d_2)$ is the CDF of the F-distribution with degrees of freedom $d_1$ and $d_2$.

For $c = 1$ (testing whether arm $j$ is superior to arm $l$):

$$p_{jl}^{\text{between}} = F_F\left(\frac{b_j}{b_l} \cdot \frac{a_l}{a_j}; 2a_j, 2a_l\right)$$

**In R**:
```r
p_between <- pf(
  q = (b_j / b_l) * (a_l / a_j),
  df1 = 2 * a_j,
  df2 = 2 * a_l
)
```

**More general form** (for HR < c):
```r
compute_p_hr_less_than_c <- function(a_j, b_j, a_l, b_l, c = 1) {
  pf(
    q = c * (b_j / b_l) * (a_l / a_j),
    df1 = 2 * a_j,
    df2 = 2 * a_l
  )
}
```

### 2.6 Summary of Posterior Computations

| Quantity | Formula | R Implementation |
|----------|---------|------------------|
| Posterior shape | $a_k = a_0 + n_k$ | `a_k <- a_0 + n_events_k` |
| Posterior rate | $b_k = b_0 + T_k$ | `b_k <- b_0 + exposure_k` |
| $P(\text{HR}_k^{\text{hist}} < c_k)$ | $F_{\text{Gamma}}(c_k \lambda_{\text{hist},k}; a_k, b_k)$ | `pgamma(c_k * lambda_hist_k, a_k, b_k)` |
| $P(\text{HR}_{jl} < c)$ | $F_F(c \cdot \frac{b_j a_l}{b_l a_j}; 2a_j, 2a_l)$ | `pf(c * b_j/b_l * a_l/a_j, 2*a_j, 2*a_l)` |

---

## 3. Decision Rules

### 3.1 Single-Arm Decision Rules

At each interim analysis, for each active arm $k$:

**Efficacy**: Arm $k$ demonstrates efficacy vs historical if:

$$p_k^{\text{single}} = P\left(\text{HR}_k^{(\text{hist})} < c_k \mid \text{data}\right) > \gamma_{\text{eff}}^{\text{single}}$$

**Futility**: Arm $k$ is dropped for futility if:

$$p_k^{\text{single}} = P\left(\text{HR}_k^{(\text{hist})} < c_k \mid \text{data}\right) < \gamma_{\text{fut}}^{\text{single}}$$

**Typical values**:
- $c_k = 0.8$ (target 20% improvement over historical)
- $\gamma_{\text{eff}}^{\text{single}} = 0.90$ (90% posterior probability for efficacy)
- $\gamma_{\text{fut}}^{\text{single}} = 0.10$ (10% posterior probability for futility)

### 3.2 Between-Arm Decision Rules

At each interim analysis (primarily in STATE_BETWEEN, but computed throughout):

**Efficacy**: Arm $j$ is superior to arm $l$ if:

$$p_{jl}^{\text{between}} = P\left(\text{HR}_{jl} < 1 \mid \text{data}\right) > \gamma_{\text{eff}}^{\text{between}}$$

**Futility**: Stop between-arm comparison for futility if:

$$p_{jl}^{\text{between}} < \gamma_{\text{fut}}^{\text{between}}$$

**Typical values**:
- $\gamma_{\text{eff}}^{\text{between}} = 0.975$ (more stringent for confirmatory)
- $\gamma_{\text{fut}}^{\text{between}} = 0.05$

### 3.3 Transition Trigger Rules

The trial moves from STATE_SINGLE to STATE_CONSIDER_CONVERSION when a trigger condition is met. Configurable options:

| Trigger Type | Condition |
|--------------|-----------|
| `"any_single_success"` | At least one arm $k$ has $p_k^{\text{single}} > \gamma_{\text{eff}}^{\text{single}}$ |
| `"all_single_success"` | All active arms have $p_k^{\text{single}} > \gamma_{\text{eff}}^{\text{single}}$ |
| `"k_of_K"` | At least $k$ of $K$ arms meet single-arm efficacy |

### 3.4 Conversion Decision Rules

At STATE_CONSIDER_CONVERSION, compute predictive probability $\pi_{\text{pred}}(N_{\text{add}})$ for candidate additional sample sizes.

**Go decision** (proceed to between-arm phase):

$$\exists N_{\text{add}} \in \mathcal{N}_{\text{candidates}}: \pi_{\text{pred}}(N_{\text{add}}) \geq \pi_{\text{go}}$$

If multiple $N_{\text{add}}$ satisfy this, select the smallest.

**No-go decision** (stop with single-arm conclusions only):

$$\max_{N_{\text{add}} \in \mathcal{N}_{\text{candidates}}} \pi_{\text{pred}}(N_{\text{add}}) < \pi_{\text{nogo}}$$

**Typical values**:
- $\pi_{\text{go}} = 0.70$ (70% chance of between-arm success to proceed)
- $\pi_{\text{nogo}} = 0.20$ (below 20% not worth continuing)
- $\mathcal{N}_{\text{candidates}} = \{30, 40, 50, \ldots, 100\}$ per arm

---

## 4. State Machine Architecture

### 4.1 State Definitions

```
typedef enum {
  STATE_SINGLE,              // Primary monitoring on within-arm comparisons
  STATE_CONSIDER_CONVERSION, // Evaluating whether to extend to between-arm
  STATE_BETWEEN,             // Full between-arm comparison phase
  STATE_STOP                 // Trial complete
} TrialState;
```

### 4.2 State Transition Logic

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        STATE TRANSITION DIAGRAM                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│                          ┌─────────────────┐                                │
│                          │  STATE_SINGLE   │                                │
│                          └────────┬────────┘                                │
│                                   │                                         │
│            ┌──────────────────────┼──────────────────────┐                  │
│            │                      │                      │                  │
│            ▼                      ▼                      ▼                  │
│   ┌────────────────┐    ┌────────────────┐    ┌─────────────────────┐      │
│   │ All arms hit   │    │ Max N reached  │    │ Transition trigger  │      │
│   │ futility       │    │ (no trigger)   │    │ met                 │      │
│   └───────┬────────┘    └───────┬────────┘    └──────────┬──────────┘      │
│           │                     │                        │                  │
│           │                     │                        ▼                  │
│           │                     │             ┌─────────────────────────┐   │
│           │                     │             │ STATE_CONSIDER_CONVERSION│   │
│           │                     │             └──────────┬──────────────┘   │
│           │                     │                        │                  │
│           │                     │         ┌──────────────┴──────────────┐   │
│           │                     │         │                             │   │
│           │                     │         ▼                             ▼   │
│           │                     │  ┌─────────────┐              ┌──────────┐│
│           │                     │  │ PP ≥ pp_go  │              │PP < nogo ││
│           │                     │  └──────┬──────┘              └────┬─────┘│
│           │                     │         │                          │      │
│           │                     │         ▼                          │      │
│           │                     │  ┌─────────────────┐               │      │
│           │                     │  │  STATE_BETWEEN  │               │      │
│           │                     │  └────────┬────────┘               │      │
│           │                     │           │                        │      │
│           │                     │    ┌──────┴──────┐                 │      │
│           │                     │    │             │                 │      │
│           │                     │    ▼             ▼                 │      │
│           │                     │ ┌──────┐    ┌────────┐             │      │
│           │                     │ │Eff/  │    │ Max N  │             │      │
│           │                     │ │Fut   │    │reached │             │      │
│           │                     │ └──┬───┘    └───┬────┘             │      │
│           │                     │    │            │                  │      │
│           ▼                     ▼    ▼            ▼                  ▼      │
│        ┌────────────────────────────────────────────────────────────────┐   │
│        │                        STATE_STOP                              │   │
│        │  • Compile final conclusions (single-arm and/or between-arm)   │   │
│        └────────────────────────────────────────────────────────────────┘   │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### 4.3 Pseudo-code for State Transitions

```r
update_trial_state <- function(trial) {
  
  switch(trial$state,
    
    "STATE_SINGLE" = {
      # Update posteriors
      trial <- update_posteriors(trial)
      
      # Check single-arm futility (drop arms)
      for (k in trial$active_arms) {
        if (trial$p_single[k] < trial$design$gamma_single_fut) {
          trial <- drop_arm(trial, k, reason = "futility")
        }
      }
      
      # Check if all arms dropped
      if (length(trial$active_arms) == 0) {
        trial$state <- "STATE_STOP"
        trial$conclusion <- "all_arms_futile"
        return(trial)
      }
      
      # Check transition trigger
      trigger_met <- check_transition_trigger(trial)
      if (trigger_met) {
        trial$state <- "STATE_CONSIDER_CONVERSION"
        return(trial)
      }
      
      # Check max sample size
      if (trial$n_enrolled >= trial$design$n_max_single_phase) {
        trial$state <- "STATE_STOP"
        trial$conclusion <- "max_n_single_phase"
        return(trial)
      }
      
      # Continue enrollment
      trial
    },
    
    "STATE_CONSIDER_CONVERSION" = {
      # Compute predictive probability for each candidate N_add
      pp_results <- compute_pp_curve(trial)
      
      # Find minimum N_add achieving pp_go (if any)
      viable_n <- pp_results$n_add[pp_results$pp >= trial$design$pp_go]
      
      if (length(viable_n) > 0) {
        # Go decision: proceed to between-arm phase
        trial$n_add_selected <- min(viable_n)
        trial$n_max_between <- trial$n_enrolled + trial$n_add_selected * length(trial$active_arms)
        trial$state <- "STATE_BETWEEN"
      } else if (max(pp_results$pp) < trial$design$pp_nogo) {
        # No-go decision: stop with single-arm conclusions
        trial$state <- "STATE_STOP"
        trial$conclusion <- "conversion_nogo"
      } else {
        # Ambiguous: could implement additional logic or default to no-go
        trial$state <- "STATE_STOP"
        trial$conclusion <- "conversion_ambiguous"
      }
      
      trial
    },
    
    "STATE_BETWEEN" = {
      # Update posteriors
      trial <- update_posteriors(trial)
      
      # Check between-arm efficacy
      if (trial$p_between > trial$design$gamma_between_eff) {
        trial$state <- "STATE_STOP"
        trial$conclusion <- "between_arm_efficacy"
        return(trial)
      }
      
      # Check between-arm futility
      if (trial$p_between < trial$design$gamma_between_fut) {
        trial$state <- "STATE_STOP"
        trial$conclusion <- "between_arm_futility"
        return(trial)
      }
      
      # Check max sample size
      if (trial$n_enrolled >= trial$n_max_between) {
        trial$state <- "STATE_STOP"
        trial$conclusion <- "max_n_between_phase"
        return(trial)
      }
      
      # Continue enrollment
      trial
    },
    
    "STATE_STOP" = {
      # Terminal state - no transitions
      trial
    }
  )
}
```

---

## 5. Predictive Probability Engine

### 5.1 Core Algorithm

The predictive probability $\pi_{\text{pred}}(N_{\text{add}})$ answers: "Given current data, if we enroll $N_{\text{add}}$ additional patients per arm, what is the probability that the final between-arm comparison will meet the efficacy criterion?"

**Algorithm**:

```
Input:
  - Current posterior parameters: (a_A, b_A), (a_B, b_B)
  - Additional patients per arm: N_add
  - Accrual rate, follow-up time
  - Success criterion: γ_between_eff
  - Number of Monte Carlo samples: n_outer

Output:
  - π_pred: Predictive probability of success

Algorithm:
  success_count ← 0
  
  FOR i = 1 TO n_outer:
    # Step 1: Draw "true" hazards from current posterior
    λ_A^(i) ~ Gamma(a_A, b_A)
    λ_B^(i) ~ Gamma(a_B, b_B)
    
    # Step 2: Simulate future events and exposure under these "true" hazards
    (n_A^future, T_A^future) ← simulate_future(λ_A^(i), N_add, accrual, followup)
    (n_B^future, T_B^future) ← simulate_future(λ_B^(i), N_add, accrual, followup)
    
    # Step 3: Compute final posterior parameters (accumulate sufficient stats)
    a_A^final ← a_A + n_A^future
    b_A^final ← b_A + T_A^future
    a_B^final ← a_B + n_B^future
    b_B^final ← b_B + T_B^future
    
    # Step 4: Compute P(HR_AB < 1 | final data) using closed form
    p_between^(i) ← F_F(b_A^final/b_B^final × a_B^final/a_A^final; 2a_A^final, 2a_B^final)
    
    # Step 5: Check success criterion
    IF p_between^(i) > γ_between_eff THEN
      success_count ← success_count + 1
  
  RETURN success_count / n_outer
```

### 5.2 Simulating Future Events

For exponential survival with hazard $\lambda$, accrual rate $r$ patients/month, and analysis at time $\tau$ after last patient enrolled:

```r
simulate_future_arm <- function(lambda, n_patients, accrual_rate, followup_months) {
  # Enrollment times (uniform accrual)
  enrollment_duration <- n_patients / accrual_rate
  enrollment_times <- sort(runif(n_patients, 0, enrollment_duration))
  
  # Analysis time
  analysis_time <- enrollment_duration + followup_months
  
  # Simulate survival times
  survival_times <- rexp(n_patients, rate = lambda)
  
  # Calendar time of event
  event_times <- enrollment_times + survival_times
  
  # Observed time (censored at analysis)
  observed_times <- pmin(survival_times, analysis_time - enrollment_times)
  event_indicators <- (event_times <= analysis_time)
  
  # Sufficient statistics
  n_events <- sum(event_indicators)
  total_exposure <- sum(observed_times)
  
  list(
    n_events = n_events,
    total_exposure = total_exposure
  )
}
```

### 5.3 Complete Predictive Probability Function

```r
compute_pp_between_success <- function(
  a_post,                    # Named vector: c(A = ..., B = ...)
  b_post,                    # Named vector: c(A = ..., B = ...)
  n_additional_per_arm,      # Scalar or vector to evaluate
  accrual_rate,              # Patients per month per arm
  followup_months,           # Follow-up after last patient
  gamma_between_eff,         # Success threshold (e.g., 0.975)
  n_outer = 1000             # Monte Carlo samples
) {
  
  # Handle vector of candidate N values
  if (length(n_additional_per_arm) > 1) {
    return(sapply(n_additional_per_arm, function(n) {
      compute_pp_between_success(
        a_post, b_post, n, accrual_rate, followup_months, gamma_between_eff, n_outer
      )
    }))
  }
  
  success_count <- 0
  
  for (i in seq_len(n_outer)) {
    
    # Step 1: Draw "true" hazards from current posterior
    lambda_A <- rgamma(1, shape = a_post["A"], rate = b_post["A"])
    lambda_B <- rgamma(1, shape = a_post["B"], rate = b_post["B"])
    
    # Step 2: Simulate future events and exposure
    fut_A <- simulate_future_arm(lambda_A, n_additional_per_arm, accrual_rate, followup_months)
    fut_B <- simulate_future_arm(lambda_B, n_additional_per_arm, accrual_rate, followup_months)
    
    # Step 3: Update posterior (accumulate sufficient stats)
    a_final <- c(
      A = a_post["A"] + fut_A$n_events,
      B = a_post["B"] + fut_B$n_events
    )
    b_final <- c(
      A = b_post["A"] + fut_A$total_exposure,
      B = b_post["B"] + fut_B$total_exposure
    )
    
    # Step 4: Compute P(HR_AB < 1 | final data) - closed form
    p_between <- pf(
      q = (b_final["A"] / b_final["B"]) * (a_final["B"] / a_final["A"]),
      df1 = 2 * a_final["A"],
      df2 = 2 * a_final["B"]
    )
    
    # Step 5: Check success criterion
    if (p_between > gamma_between_eff) {
      success_count <- success_count + 1
    }
  }
  
  success_count / n_outer
}
```

### 5.4 Finding Required Sample Size

```r
find_n_for_target_pp <- function(
  a_post, b_post,
  target_pp,                  # Target predictive probability (e.g., 0.80)
  accrual_rate,
  followup_months,
  gamma_between_eff,
  n_range = c(10, 200),       # Search range for N per arm
  tolerance = 5,              # Tolerance for binary search
  n_outer = 1000
) {
  
  # Binary search for minimum N achieving target PP
  lower <- n_range[1]
  upper <- n_range[2]
  
  # Check boundaries
  pp_lower <- compute_pp_between_success(
    a_post, b_post, lower, accrual_rate, followup_months, gamma_between_eff, n_outer
  )
  pp_upper <- compute_pp_between_success(
    a_post, b_post, upper, accrual_rate, followup_months, gamma_between_eff, n_outer
  )
  
  if (pp_lower >= target_pp) {
    return(list(n_required = lower, pp_achieved = pp_lower, converged = TRUE))
  }
  if (pp_upper < target_pp) {
    return(list(n_required = NA, pp_achieved = pp_upper, converged = FALSE,
                message = "Target PP not achievable within range"))
  }
  
  # Binary search
  while ((upper - lower) > tolerance) {
    mid <- round((lower + upper) / 2)
    pp_mid <- compute_pp_between_success(
      a_post, b_post, mid, accrual_rate, followup_months, gamma_between_eff, n_outer
    )
    
    if (pp_mid >= target_pp) {
      upper <- mid
    } else {
      lower <- mid
    }
  }
  
  list(
    n_required = upper,
    pp_achieved = compute_pp_between_success(
      a_post, b_post, upper, accrual_rate, followup_months, gamma_between_eff, n_outer * 2
    ),
    converged = TRUE
  )
}
```

### 5.5 Computing the Full PP Curve

```r
compute_pp_curve <- function(trial) {
  # Extract current posterior parameters
  a_post <- trial$posterior$a
  b_post <- trial$posterior$b
  
  # Candidate additional sample sizes
  n_candidates <- trial$design$conversion$n_max_candidates
  
  # Compute PP for each candidate
  pp_values <- compute_pp_between_success(
    a_post = a_post,
    b_post = b_post,
    n_additional_per_arm = n_candidates,
    accrual_rate = trial$design$accrual_rate,
    followup_months = trial$design$followup_months,
    gamma_between_eff = trial$design$between_arm_criteria$eff_prob,
    n_outer = trial$design$conversion$n_outer_pp
  )
  
  data.frame(
    n_add = n_candidates,
    pp = pp_values
  )
}
```

---

## 6. Software Architecture

### 6.1 Package Structure

```
evolveTrial/
├── R/
│   ├── hybrid_design.R           # Design specification
│   ├── hybrid_trial.R            # Trial state and simulation
│   ├── hybrid_posterior.R        # Posterior computations (closed-form)
│   ├── hybrid_decisions.R        # Decision rule implementations
│   ├── hybrid_predictive.R       # Predictive probability engine
│   ├── hybrid_simulate.R         # Monte Carlo simulation wrapper
│   ├── hybrid_summary.R          # Operating characteristics summary
│   ├── single_arm_trial.R        # Existing (minimal changes)
│   ├── multiarm_trial.R          # Existing (minimal changes)
│   └── utils.R                   # Shared utilities
├── tests/
│   ├── test_hybrid_posterior.R   # Unit tests for posterior computations
│   ├── test_hybrid_decisions.R   # Unit tests for decision rules
│   ├── test_hybrid_predictive.R  # Unit tests for PP engine
│   ├── test_hybrid_simulate.R    # Integration tests
│   └── test_hybrid_edge_cases.R  # Edge case tests
├── vignettes/
│   └── hybrid_design_workflow.Rmd
└── man/
    └── ...
```

### 6.2 Core S3 Classes

```r
# =============================================================================
# DESIGN SPECIFICATION CLASS
# =============================================================================

#' Create a Hybrid Single-to-Between Arm Survival Design
#'
#' @param arms Character vector of arm names
#' @param hist_median_surv Named numeric vector of historical median survival (months)
#' @param single_arm_criteria List with hr_threshold, eff_prob, fut_prob
#' @param between_arm_criteria List with eff_prob, fut_prob
#' @param conversion List with trigger, pp_go, pp_nogo, n_max_candidates, etc.
#' @param constraints List with n_max_total, max_duration_months, etc.
#' @param interims List with timing, events_per_interim, etc.
#' @param randomization List with scheme, initial_equal_n
#' @param prior List with a0, b0 (gamma prior parameters)
#' 
#' @return Object of class "hybrid_surv_design"
#' @export
create_hybrid_surv_design <- function(
  arms = c("A", "B"),
  
  hist_median_surv = c(A = 12, B = 12),
  
  single_arm_criteria = list(
    hr_threshold = 0.8,
    eff_prob = 0.90,
    fut_prob = 0.10
  ),
  
  between_arm_criteria = list(
    eff_prob = 0.975,
    fut_prob = 0.05
  ),
  
  conversion = list(
    trigger = "any_single_success",
    pp_go = 0.70,
    pp_nogo = 0.20,
    n_max_candidates = seq(30, 100, by = 10),
    n_outer_pp = 1000
  ),
  
  constraints = list(
    n_max_total = 200,
    n_max_single_phase = 100,
    max_duration_months = 48
  ),
  
  interims = list(
    timing = "event_driven",
    events_per_interim = 20,
    max_interims = 10
  ),
  
  randomization = list(
    scheme = "equal",
    initial_equal_n = 20
  ),
  
  prior = list(
    a0 = 0.001,
    b0 = 0.001
  ),
  
  accrual_rate = 5,
  followup_months = 12
) {
  

  # Validate inputs
  stopifnot(length(arms) >= 2)
  stopifnot(all(arms %in% names(hist_median_surv)))
  stopifnot(single_arm_criteria$hr_threshold > 0 && single_arm_criteria$hr_threshold < 1)
  stopifnot(conversion$trigger %in% c("any_single_success", "all_single_success", "k_of_K"))
  
  # Compute historical hazards from median survival
  hist_hazard <- setNames(log(2) / hist_median_surv[arms], arms)
  
  # Construct design object
  design <- list(
    arms = arms,
    K = length(arms),
    hist_median_surv = hist_median_surv[arms],
    hist_hazard = hist_hazard,
    single_arm_criteria = single_arm_criteria,
    between_arm_criteria = between_arm_criteria,
    conversion = conversion,
    constraints = constraints,
    interims = interims,
    randomization = randomization,
    prior = prior,
    accrual_rate = accrual_rate,
    followup_months = followup_months
  )
  
  class(design) <- c("hybrid_surv_design", "list")
  design
}

# Print method
print.hybrid_surv_design <- function(x, ...) {
  cat("Hybrid Single-to-Between Arm Survival Design\n")
  cat("=============================================\n\n")
  cat("Arms:", paste(x$arms, collapse = ", "), "\n")
  cat("Historical median survival:", paste(x$hist_median_surv, "months", collapse = ", "), "\n\n")
  cat("Single-arm criteria:\n")
  cat("  HR threshold:", x$single_arm_criteria$hr_threshold, "\n")
  cat("  Efficacy prob:", x$single_arm_criteria$eff_prob, "\n")
  cat("  Futility prob:", x$single_arm_criteria$fut_prob, "\n\n")
  cat("Between-arm criteria:\n")
  cat("  Efficacy prob:", x$between_arm_criteria$eff_prob, "\n")
  cat("  Futility prob:", x$between_arm_criteria$fut_prob, "\n\n")
  cat("Conversion:\n")
  cat("  Trigger:", x$conversion$trigger, "\n")
  cat("  PP go threshold:", x$conversion$pp_go, "\n")
  cat("  PP no-go threshold:", x$conversion$pp_nogo, "\n")
  invisible(x)
}

# =============================================================================
# TRIAL STATE CLASS
# =============================================================================

#' Initialize a Hybrid Trial
#'
#' @param design Object of class "hybrid_surv_design"
#' @return Object of class "hybrid_trial"
initialize_hybrid_trial <- function(design) {
  
  trial <- list(
    design = design,
    state = "STATE_SINGLE",
    
    # Data
    patients = data.frame(
      id = integer(),
      arm = character(),
      enrollment_time = numeric(),
      event_time = numeric(),
      observed_time = numeric(),
      event = logical()
    ),
    n_enrolled = setNames(rep(0L, design$K), design$arms),
    
    # Sufficient statistics for posterior
    posterior = list(
      a = setNames(rep(design$prior$a0, design$K), design$arms),
      b = setNames(rep(design$prior$b0, design$K), design$arms)
    ),
    
    # Current probabilities
    p_single = setNames(rep(NA_real_, design$K), design$arms),
    p_between = NA_real_,
    
    # Arm status
    active_arms = design$arms,
    dropped_arms = character(),
    successful_arms = character(),  # Single-arm efficacy achieved
    
    # Conversion info
    conversion_triggered = FALSE,
    n_add_selected = NA_integer_,
    n_max_between = NA_integer_,
    pp_at_conversion = NA_real_,
    
    # Timeline
    current_time = 0,
    interim_count = 0,
    
    # Final conclusions
    conclusion = NULL,
    final_p_single = NULL,
    final_p_between = NULL
  )
  
  class(trial) <- c("hybrid_trial", "list")
  trial
}
```

### 6.3 Core Functions

```r
# =============================================================================
# POSTERIOR COMPUTATIONS
# =============================================================================

#' Update Posterior Parameters from Trial Data
#'
#' @param trial Object of class "hybrid_trial"
#' @param analysis_time Time at which to compute posterior (for censoring)
#' @return Updated trial object
update_posteriors <- function(trial, analysis_time = NULL) {
  
  if (is.null(analysis_time)) {
    analysis_time <- trial$current_time
  }
  
  for (k in trial$design$arms) {
    # Get patients in this arm
    arm_data <- trial$patients[trial$patients$arm == k, ]
    
    if (nrow(arm_data) == 0) {
      trial$posterior$a[k] <- trial$design$prior$a0
      trial$posterior$b[k] <- trial$design$prior$b0
    } else {
      # Compute observed times (censored at analysis_time)
      time_on_study <- pmin(
        arm_data$event_time - arm_data$enrollment_time,
        analysis_time - arm_data$enrollment_time
      )
      time_on_study <- pmax(time_on_study, 0)  # Ensure non-negative
      
      events <- arm_data$event & (arm_data$event_time <= analysis_time)
      
      # Update sufficient statistics
      n_events <- sum(events)
      total_exposure <- sum(time_on_study)
      
      trial$posterior$a[k] <- trial$design$prior$a0 + n_events
      trial$posterior$b[k] <- trial$design$prior$b0 + total_exposure
    }
  }
  
  # Compute decision quantities
  trial <- compute_decision_quantities(trial)
  
  trial
}

#' Compute Decision Quantities (Single-arm and Between-arm Probabilities)
#'
#' @param trial Object of class "hybrid_trial"
#' @return Updated trial object with p_single and p_between
compute_decision_quantities <- function(trial) {
  
  # Single-arm probabilities: P(HR_k^hist < c_k | data)
  for (k in trial$design$arms) {
    c_k <- trial$design$single_arm_criteria$hr_threshold
    lambda_hist_k <- trial$design$hist_hazard[k]
    
    # Threshold on lambda scale
    lambda_threshold <- c_k * lambda_hist_k
    
    # P(lambda_k < threshold) = Gamma CDF
    trial$p_single[k] <- pgamma(
      q = lambda_threshold,
      shape = trial$posterior$a[k],
      rate = trial$posterior$b[k]
    )
  }
  
  # Between-arm probability (for two-arm case: P(HR_AB < 1))
  # For K > 2, could compute pairwise or vs "best"
  if (trial$design$K == 2) {
    arms <- trial$design$arms
    a_A <- trial$posterior$a[arms[1]]
    b_A <- trial$posterior$b[arms[1]]
    a_B <- trial$posterior$a[arms[2]]
    b_B <- trial$posterior$b[arms[2]]
    
    # P(lambda_A / lambda_B < 1) using F distribution
    trial$p_between <- pf(
      q = (b_A / b_B) * (a_B / a_A),
      df1 = 2 * a_A,
      df2 = 2 * a_B
    )
  } else {
    # For K > 2, implement pairwise or other comparison strategy
    # Placeholder: compare first arm vs pooled others, or implement max comparison
    trial$p_between <- NA_real_
  }
  
  trial
}

# =============================================================================
# DECISION RULES
# =============================================================================

#' Check Single-Arm Efficacy for an Arm
#' @return Logical
check_single_arm_efficacy <- function(trial, arm) {
  trial$p_single[arm] > trial$design$single_arm_criteria$eff_prob
}

#' Check Single-Arm Futility for an Arm
#' @return Logical
check_single_arm_futility <- function(trial, arm) {
  trial$p_single[arm] < trial$design$single_arm_criteria$fut_prob
}

#' Check Between-Arm Efficacy
#' @return Logical
check_between_arm_efficacy <- function(trial) {
  trial$p_between > trial$design$between_arm_criteria$eff_prob
}

#' Check Between-Arm Futility
#' @return Logical
check_between_arm_futility <- function(trial) {
  trial$p_between < trial$design$between_arm_criteria$fut_prob
}

#' Check Transition Trigger
#' @return Logical
check_transition_trigger <- function(trial) {
  trigger <- trial$design$conversion$trigger
  
  switch(trigger,
    "any_single_success" = {
      any(sapply(trial$active_arms, function(k) check_single_arm_efficacy(trial, k)))
    },
    "all_single_success" = {
      all(sapply(trial$active_arms, function(k) check_single_arm_efficacy(trial, k)))
    },
    "k_of_K" = {
      k_required <- trial$design$conversion$k_required
      sum(sapply(trial$active_arms, function(k) check_single_arm_efficacy(trial, k))) >= k_required
    },
    FALSE
  )
}

#' Drop an Arm
#' @return Updated trial object
drop_arm <- function(trial, arm, reason = "futility") {
  trial$active_arms <- setdiff(trial$active_arms, arm)
  trial$dropped_arms <- c(trial$dropped_arms, arm)
  attr(trial$dropped_arms, "reasons") <- c(
    attr(trial$dropped_arms, "reasons"),
    setNames(reason, arm)
  )
  trial
}
```

### 6.4 Simulation Engine

```r
# =============================================================================
# SINGLE TRIAL SIMULATION
# =============================================================================

#' Simulate a Single Hybrid Trial
#'
#' @param design Object of class "hybrid_surv_design"
#' @param scenario List with true_hazard (named vector) and other simulation params
#' @param seed Random seed
#' @return Object of class "hybrid_trial" with complete trajectory
simulate_hybrid_trial <- function(design, scenario, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Initialize trial
  trial <- initialize_hybrid_trial(design)
  
  # Run until terminal state
 while (trial$state != "STATE_STOP") {
    
    # Enroll patients until next interim
    trial <- enroll_patients_until_interim(trial, scenario)
    
    # Update posteriors at interim
    trial <- update_posteriors(trial)
    trial$interim_count <- trial$interim_count + 1
    
    # Apply state-specific logic and transitions
    trial <- update_trial_state(trial)
  }
  
  # Finalize
  trial$final_p_single <- trial$p_single
  trial$final_p_between <- trial$p_between
  
  trial
}

#' Enroll Patients Until Next Interim Analysis
#'
#' @param trial Object of class "hybrid_trial"
#' @param scenario Simulation scenario with true hazards
#' @return Updated trial object
enroll_patients_until_interim <- function(trial, scenario) {
  
  # Determine target for next interim
  if (trial$design$interims$timing == "event_driven") {
    # Enroll until we have enough events
    target_events <- trial$design$interims$events_per_interim * (trial$interim_count + 1)
    
    while (count_events(trial) < target_events && 
           sum(trial$n_enrolled) < trial$design$constraints$n_max_total) {
      trial <- enroll_one_patient(trial, scenario)
      trial$current_time <- trial$current_time + 1 / trial$design$accrual_rate
    }
    
  } else if (trial$design$interims$timing == "calendar") {
    # Enroll until calendar time
    target_time <- trial$design$interims$months_per_interim * (trial$interim_count + 1)
    
    while (trial$current_time < target_time && 
           sum(trial$n_enrolled) < trial$design$constraints$n_max_total) {
      trial <- enroll_one_patient(trial, scenario)
      trial$current_time <- trial$current_time + 1 / trial$design$accrual_rate
    }
  }
  
  trial
}

#' Enroll One Patient
#'
#' @param trial Object of class "hybrid_trial"
#' @param scenario Simulation scenario
#' @return Updated trial object
enroll_one_patient <- function(trial, scenario) {
  
  # Randomize to an active arm
  if (trial$design$randomization$scheme == "equal") {
    arm <- sample(trial$active_arms, 1)
  } else {
    # Implement RAR if needed
    arm <- sample(trial$active_arms, 1)
  }
  
  # Generate survival time under true hazard
  true_hazard <- scenario$true_hazard[arm]
  survival_time <- rexp(1, rate = true_hazard)
  
  # Create patient record
  new_patient <- data.frame(
    id = nrow(trial$patients) + 1,
    arm = arm,
    enrollment_time = trial$current_time,
    event_time = trial$current_time + survival_time,
    observed_time = NA,  # Will be computed at analysis
    event = TRUE,  # Will be updated for censoring
    stringsAsFactors = FALSE
  )
  
  trial$patients <- rbind(trial$patients, new_patient)
  trial$n_enrolled[arm] <- trial$n_enrolled[arm] + 1
  
  trial
}

#' Count Total Events Observed by Current Time
count_events <- function(trial) {
  sum(trial$patients$event_time <= trial$current_time)
}

# =============================================================================
# MULTIPLE TRIAL SIMULATION
# =============================================================================

#' Simulate Multiple Hybrid Trials and Compute Operating Characteristics
#'
#' @param design Object of class "hybrid_surv_design"
#' @param scenario Simulation scenario
#' @param n_sims Number of simulations
#' @param parallel Logical, use parallel processing
#' @param n_cores Number of cores (NULL for auto-detect)
#' @return Object of class "hybrid_sim_results"
#' @export
simulate_hybrid_trials <- function(
  design,
  scenario,
  n_sims = 1000,
  parallel = TRUE,
  n_cores = NULL
) {
  
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    # Export required objects
    parallel::clusterExport(cl, c("design", "scenario"), envir = environment())
    parallel::clusterEvalQ(cl, library(evolveTrial))
    
    results <- parallel::parLapply(cl, 1:n_sims, function(i) {
      simulate_hybrid_trial(design, scenario, seed = i)
    })
    
  } else {
    results <- lapply(1:n_sims, function(i) {
      simulate_hybrid_trial(design, scenario, seed = i)
    })
  }
  
  # Compile operating characteristics
  ocs <- compile_operating_characteristics(results, design, scenario)
  
  structure(
    list(
      design = design,
      scenario = scenario,
      n_sims = n_sims,
      trials = results,
      ocs = ocs
    ),
    class = c("hybrid_sim_results", "list")
  )
}

#' Compile Operating Characteristics from Simulation Results
compile_operating_characteristics <- function(results, design, scenario) {
  
  n_sims <- length(results)
  
  # Extract outcomes
  conclusions <- sapply(results, function(x) x$conclusion)
  states_reached <- sapply(results, function(x) x$conversion_triggered)
  n_enrolled_total <- sapply(results, function(x) sum(x$n_enrolled))
  n_interims <- sapply(results, function(x) x$interim_count)
  
  # Single-arm outcomes
  single_arm_efficacy <- sapply(design$arms, function(k) {
    mean(sapply(results, function(x) k %in% x$successful_arms))
  })
  
  # Between-arm outcomes
  reached_between <- mean(states_reached)
  between_efficacy_given_reached <- mean(
    sapply(results[states_reached], function(x) {
      x$conclusion == "between_arm_efficacy"
    })
  )
  between_efficacy_overall <- mean(conclusions == "between_arm_efficacy")
  
  # Type I error (depends on scenario)
  # If scenario is null (HR = 1), then between_efficacy is type I error
  
  list(
    # Sample size
    mean_n = mean(n_enrolled_total),
    sd_n = sd(n_enrolled_total),
    median_n = median(n_enrolled_total),
    
    # Single-arm
    prob_single_arm_efficacy = single_arm_efficacy,
    
    # Conversion
    prob_conversion = reached_between,
    mean_n_at_conversion = mean(sapply(results[states_reached], function(x) sum(x$n_enrolled))),
    
    # Between-arm
    prob_between_efficacy_given_conversion = between_efficacy_given_reached,
    prob_between_efficacy_overall = between_efficacy_overall,
    
    # By conclusion
    conclusion_table = table(conclusions) / n_sims,
    
    # Duration
    mean_interims = mean(n_interims)
  )
}
```

---

## 7. BATON Calibration

### 7.1 Design Parameter Vector (φ)

The parameters to be optimized by BATON:

| Parameter | Symbol | Range | Description |
|-----------|--------|-------|-------------|
| Single-arm efficacy threshold | $\gamma_{\text{eff}}^{\text{single}}$ | [0.80, 0.99] | Posterior prob for single-arm efficacy |
| Single-arm futility threshold | $\gamma_{\text{fut}}^{\text{single}}$ | [0.01, 0.20] | Posterior prob for single-arm futility |
| HR threshold vs historical | $c$ | [0.60, 0.90] | Target HR improvement |
| Between-arm efficacy threshold | $\gamma_{\text{eff}}^{\text{between}}$ | [0.95, 0.999] | Posterior prob for between-arm efficacy |
| Between-arm futility threshold | $\gamma_{\text{fut}}^{\text{between}}$ | [0.01, 0.10] | Posterior prob for between-arm futility |
| Conversion go threshold | $\pi_{\text{go}}$ | [0.50, 0.90] | PP threshold to proceed |
| Conversion no-go threshold | $\pi_{\text{nogo}}$ | [0.10, 0.40] | PP threshold to stop |

**Total dimensions**: 7 continuous parameters

### 7.2 Scenarios for Calibration

```r
# Standard scenario set for calibration
calibration_scenarios <- list(
  
  # Null scenarios (for Type I error control)
  null_global = list(
    name = "Global null",
    true_hazard = c(A = 0.0578, B = 0.0578),  # Both equal to historical (12 mo median)
    true_hr_vs_hist = c(A = 1.0, B = 1.0),
    true_hr_between = 1.0
  ),
  
  null_between = list(
    name = "Both beat historical, no between-arm difference",
    true_hazard = c(A = 0.0462, B = 0.0462),  # Both 20% better than historical (15 mo median)
    true_hr_vs_hist = c(A = 0.8, B = 0.8),
    true_hr_between = 1.0
  ),
  
  # Alternative scenarios (for power)
  alt_both_different = list(
    name = "A better than B, both beat historical",
    true_hazard = c(A = 0.0385, B = 0.0513),  # A: 18mo median, B: 13.5mo median
    true_hr_vs_hist = c(A = 0.67, B = 0.89),
    true_hr_between = 0.75  # A is 25% better than B
  ),
  
  alt_strong_difference = list(
    name = "Strong between-arm difference",
    true_hazard = c(A = 0.0347, B = 0.0578),  # A: 20mo median, B: 12mo median
    true_hr_vs_hist = c(A = 0.6, B = 1.0),
    true_hr_between = 0.6
  ),
  
  # Edge cases
  one_arm_futile = list(
    name = "One arm futile",
    true_hazard = c(A = 0.0385, B = 0.0694),  # A: 18mo, B: 10mo (worse than hist)
    true_hr_vs_hist = c(A = 0.67, B = 1.2),
    true_hr_between = 0.55
  )
)
```

### 7.3 Objective Function

```r
#' BATON Objective Function for Hybrid Design Calibration
#'
#' @param phi Named vector of design parameters
#' @param design_template Base design object
#' @param scenarios List of calibration scenarios
#' @param targets List of target operating characteristics
#' @param n_sims Simulations per scenario
#' @return List with objective value and constraint violations
baton_objective <- function(
  phi,
  design_template,
  scenarios,
  targets,
  n_sims = 5000
) {
  
  # Build design from phi
  design <- build_design_from_phi(design_template, phi)
  
  # Evaluate across scenarios
  ocs_by_scenario <- lapply(scenarios, function(scen) {
    sim_results <- simulate_hybrid_trials(design, scen, n_sims = n_sims, parallel = TRUE)
    sim_results$ocs
  })
  names(ocs_by_scenario) <- names(scenarios)
  
  # Compute constraint violations
  constraints <- list()
  
  # Type I error for single-arm (under global null)
  constraints$single_type1 <- max(ocs_by_scenario$null_global$prob_single_arm_efficacy) - targets$single_type1_max
  
  # Type I error for between-arm (under null_between where single-arm might succeed)
  constraints$between_type1 <- ocs_by_scenario$null_between$prob_between_efficacy_overall - targets$between_type1_max
  
  # Power for single-arm (under alt_both_different)
  constraints$single_power <- targets$single_power_min - min(ocs_by_scenario$alt_both_different$prob_single_arm_efficacy)
  
  # Power for between-arm (under alt_both_different)
  constraints$between_power <- targets$between_power_min - ocs_by_scenario$alt_both_different$prob_between_efficacy_overall
  
  # Compute objective (minimize expected sample size, subject to constraints)
  # Weighted average across scenarios
  expected_n <- weighted.mean(
    sapply(ocs_by_scenario, function(x) x$mean_n),
    w = targets$scenario_weights
  )
  
  # Add penalty for unnecessary conversions under null_between
  conversion_penalty <- ocs_by_scenario$null_between$prob_conversion * targets$conversion_penalty_weight
  
  objective_value <- expected_n + conversion_penalty
  
  list(
    value = objective_value,
    constraints = constraints,
    feasible = all(sapply(constraints, function(x) x <= 0)),
    ocs = ocs_by_scenario
  )
}

#' Build Design Object from Parameter Vector
build_design_from_phi <- function(design_template, phi) {
  
  design <- design_template
  
  design$single_arm_criteria$eff_prob <- phi["gamma_single_eff"]
  design$single_arm_criteria$fut_prob <- phi["gamma_single_fut"]
  design$single_arm_criteria$hr_threshold <- phi["hr_threshold"]
  design$between_arm_criteria$eff_prob <- phi["gamma_between_eff"]
  design$between_arm_criteria$fut_prob <- phi["gamma_between_fut"]
  design$conversion$pp_go <- phi["pp_go"]
  design$conversion$pp_nogo <- phi["pp_nogo"]
  
  design
}
```

### 7.4 Calibration Wrapper

```r
#' Calibrate Hybrid Design Using BATON
#'
#' @param design_template Base design object
#' @param scenarios Calibration scenarios
#' @param targets Target operating characteristics
#' @param baton_config BATON configuration
#' @return Calibration results
#' @export
calibrate_hybrid_design <- function(
  design_template,
  
  scenarios = calibration_scenarios,
  
  targets = list(
    single_type1_max = 0.10,
    between_type1_max = 0.05,
    single_power_min = 0.80,
    between_power_min = 0.80,
    scenario_weights = c(0.2, 0.2, 0.3, 0.2, 0.1),
    conversion_penalty_weight = 10
  ),
  
  baton_config = list(
    n_initial = 20,
    n_iterations = 150,
    n_sims_per_design = 3000,
    n_validation_sims = 20000
  )
) {
  
  # Define parameter bounds
  bounds <- list(
    gamma_single_eff = c(0.80, 0.99),
    gamma_single_fut = c(0.01, 0.20),
    hr_threshold = c(0.60, 0.90),
    gamma_between_eff = c(0.95, 0.999),
    gamma_between_fut = c(0.01, 0.10),
    pp_go = c(0.50, 0.90),
    pp_nogo = c(0.10, 0.40)
  )
  
  # Create objective wrapper for BATON
  objective_wrapper <- function(phi) {
    baton_objective(
      phi = phi,
      design_template = design_template,
      scenarios = scenarios,
      targets = targets,
      n_sims = baton_config$n_sims_per_design
    )
  }
  
  # Run BATON optimization
  # (Assuming baton_optimize is available from BATON package)
  result <- baton_optimize(
    objective = objective_wrapper,
    bounds = bounds,
    n_initial = baton_config$n_initial,
    n_iterations = baton_config$n_iterations,
    constrained = TRUE
  )
  
  # Validate final design
  final_design <- build_design_from_phi(design_template, result$best_phi)
  
  validation <- lapply(scenarios, function(scen) {
    sim_results <- simulate_hybrid_trials(
      final_design, scen, 
      n_sims = baton_config$n_validation_sims, 
      parallel = TRUE
    )
    sim_results$ocs
  })
  
  # Local sensitivity analysis
  sensitivity <- validate_local_neighborhood(
    result$best_phi, design_template, scenarios, targets
  )
  
  list(
    optimal_phi = result$best_phi,
    optimal_design = final_design,
    calibration_trace = result$trace,
    validation = validation,
    sensitivity = sensitivity,
    targets = targets
  )
}
```

### 7.5 Local Validation

```r
#' Validate Design in Local Neighborhood of Optimal
#'
#' @param phi Optimal parameter vector
#' @param design_template Base design
#' @param scenarios Calibration scenarios
#' @param targets Target OCs
#' @param perturbation Relative perturbation (default 5%)
#' @param n_sims Simulations per perturbed design
validate_local_neighborhood <- function(
  phi,
  design_template,
  scenarios,
  targets,
  perturbation = 0.05,
  n_sims = 10000
) {
  
  results <- list()
  
  for (param in names(phi)) {
    # Perturb down
    phi_low <- phi
    phi_low[param] <- phi[param] * (1 - perturbation)
    
    # Perturb up
    phi_high <- phi
    phi_high[param] <- phi[param] * (1 + perturbation)
    
    # Evaluate both
    oc_low <- baton_objective(phi_low, design_template, scenarios, targets, n_sims)
    oc_high <- baton_objective(phi_high, design_template, scenarios, targets, n_sims)
    
    results[[param]] <- list(
      phi_low = phi_low[param],
      phi_optimal = phi[param],
      phi_high = phi_high[param],
      feasible_low = oc_low$feasible,
      feasible_high = oc_high$feasible,
      objective_low = oc_low$value,
      objective_high = oc_high$value,
      sensitivity = (oc_high$value - oc_low$value) / (2 * perturbation * phi[param])
    )
  }
  
  # Summary: is the design robust?
  all_feasible <- all(sapply(results, function(x) x$feasible_low && x$feasible_high))
  max_sensitivity <- max(abs(sapply(results, function(x) x$sensitivity)))
  
  list(
    by_parameter = results,
    all_neighbors_feasible = all_feasible,
    max_sensitivity = max_sensitivity,
    robust = all_feasible && max_sensitivity < 100  # Threshold for "robust"
  )
}
```

---

## 8. Implementation Phases

### Phase 1: Core Infrastructure (Week 1-2)

**Objectives**:
- Implement posterior computation functions
- Implement decision rule functions
- Create design and trial S3 classes

**Deliverables**:
```
✓ create_hybrid_surv_design()
✓ initialize_hybrid_trial()
✓ update_posteriors()
✓ compute_decision_quantities()
✓ compute_p_hr_less_than_c()
✓ check_single_arm_efficacy()
✓ check_single_arm_futility()
✓ check_between_arm_efficacy()
✓ check_between_arm_futility()
✓ check_transition_trigger()
```

**Unit tests**:
- Verify posterior computations match Monte Carlo for small examples
- Verify F-distribution formula for HR comparisons
- Verify gamma CDF formula for single-arm comparisons

### Phase 2: Predictive Probability Engine (Week 2-3)

**Objectives**:
- Implement future data simulation
- Implement PP computation
- Implement sample size finding

**Deliverables**:
```
✓ simulate_future_arm()
✓ compute_pp_between_success()
✓ find_n_for_target_pp()
✓ compute_pp_curve()
```

**Unit tests**:
- Verify PP increases with sample size
- Verify PP approaches 1 when true effect is large
- Verify PP approaches 0 when true effect is null
- Benchmark computation time

### Phase 3: State Machine and Simulation (Week 3-4)

**Objectives**:
- Implement state transition logic
- Implement single trial simulation
- Implement multi-trial simulation

**Deliverables**:
```
✓ update_trial_state()
✓ enroll_patients_until_interim()
✓ simulate_hybrid_trial()
✓ simulate_hybrid_trials()
✓ compile_operating_characteristics()
```

**Integration tests**:
- Verify state transitions occur correctly
- Verify that with extreme thresholds, design reduces to pure single-arm or pure between-arm
- Verify parallelization produces consistent results

### Phase 4: BATON Integration (Week 4-5)

**Objectives**:
- Define calibration scenarios
- Implement objective function
- Implement calibration wrapper
- Implement local validation

**Deliverables**:
```
✓ baton_objective()
✓ build_design_from_phi()
✓ calibrate_hybrid_design()
✓ validate_local_neighborhood()
```

**Validation**:
- Run calibration on test problem
- Verify type I error control
- Verify power targets achieved
- Document calibration trace

### Phase 5: Documentation and Testing (Week 5-6)

**Objectives**:
- Write vignette
- Complete unit test coverage
- Performance optimization

**Deliverables**:
```
✓ Vignette: hybrid_design_workflow.Rmd
✓ Full test coverage (>90%)
✓ Performance benchmarks
✓ Documentation for all exported functions
```

---

## 9. Testing Strategy

### 9.1 Unit Tests

```r
# tests/test_hybrid_posterior.R

test_that("Gamma CDF formula for single-arm is correct", {
  # Compare closed form to Monte Carlo
  a <- 10
  b <- 5
  threshold <- 1.5
  
  # Closed form
  p_closed <- pgamma(threshold, shape = a, rate = b)
  
  # Monte Carlo
  set.seed(123)
  samples <- rgamma(100000, shape = a, rate = b)
  p_mc <- mean(samples < threshold)
  
  expect_equal(p_closed, p_mc, tolerance = 0.01)
})

test_that("F distribution formula for between-arm HR is correct", {
  # Compare closed form to Monte Carlo
  a_A <- 15
  b_A <- 8
  a_B <- 12
  b_B <- 6
  
  # Closed form
  p_closed <- pf(
    q = (b_A / b_B) * (a_B / a_A),
    df1 = 2 * a_A,
    df2 = 2 * a_B
  )
  
  # Monte Carlo
  set.seed(123)
  lambda_A <- rgamma(100000, shape = a_A, rate = b_A)
  lambda_B <- rgamma(100000, shape = a_B, rate = b_B)
  p_mc <- mean(lambda_A / lambda_B < 1)
  
  expect_equal(p_closed, p_mc, tolerance = 0.01)
})

test_that("Posterior update accumulates sufficient statistics correctly", {
  design <- create_hybrid_surv_design()
  trial <- initialize_hybrid_trial(design)
  
  # Add some patients manually
  trial$patients <- data.frame(
    id = 1:10,
    arm = rep(c("A", "B"), each = 5),
    enrollment_time = 0,
    event_time = c(2, 3, 5, 8, 10, 1, 4, 6, 7, 9),
    observed_time = NA,
    event = TRUE
  )
  trial$current_time <- 12
  
  trial <- update_posteriors(trial)
  
  # Check sufficient statistics
  expect_equal(trial$posterior$a["A"], design$prior$a0 + 5)
  expect_equal(trial$posterior$a["B"], design$prior$a0 + 5)
})
```

### 9.2 Integration Tests

```r
# tests/test_hybrid_simulate.R

test_that("Design reduces to single-arm when between-arm threshold is impossible", {
  design <- create_hybrid_surv_design(
    between_arm_criteria = list(eff_prob = 1.0, fut_prob = 0.0)  # Never achievable
  )
  
  scenario <- list(true_hazard = c(A = 0.04, B = 0.06))
  
  results <- simulate_hybrid_trials(design, scenario, n_sims = 100)
  
  # Should never reach between-arm efficacy
  expect_equal(results$ocs$prob_between_efficacy_overall, 0)
})

test_that("Type I error is controlled under null scenario", {
  skip_on_cran()  # Long-running test
  
  design <- create_hybrid_surv_design()  # Default calibrated design
  
  null_scenario <- list(true_hazard = c(A = 0.0578, B = 0.0578))  # HR = 1
  
  results <- simulate_hybrid_trials(design, null_scenario, n_sims = 10000)
  
  # Between-arm type I error should be < 0.05 (or whatever target)
  expect_lt(results$ocs$prob_between_efficacy_overall, 0.06)  # Allow small margin
})
```

### 9.3 Edge Case Tests

```r
# tests/test_hybrid_edge_cases.R

test_that("Trial handles all arms dropped for futility", {
  design <- create_hybrid_surv_design(
    single_arm_criteria = list(hr_threshold = 0.5, eff_prob = 0.9, fut_prob = 0.5)
  )
  
  # Scenario where both arms are worse than historical
  bad_scenario <- list(true_hazard = c(A = 0.08, B = 0.08))  # Worse than 0.0578
  
  trial <- simulate_hybrid_trial(design, bad_scenario, seed = 123)
  
  expect_equal(trial$state, "STATE_STOP")
  expect_true(trial$conclusion %in% c("all_arms_futile", "max_n_single_phase"))
})

test_that("Trial handles conversion at first interim", {
  design <- create_hybrid_surv_design(
    single_arm_criteria = list(hr_threshold = 0.95, eff_prob = 0.5, fut_prob = 0.01)
  )
  
  # Scenario with strong effect
  strong_scenario <- list(true_hazard = c(A = 0.03, B = 0.04))
  
  trial <- simulate_hybrid_trial(design, strong_scenario, seed = 123)
  
  # Should trigger conversion early
  expect_true(trial$conversion_triggered)
})
```

---

## Appendix A: Mathematical Derivations

### A.1 Derivation of P(HR < c) Using F Distribution

Let $\lambda_A \sim \text{Gamma}(a_A, b_A)$ and $\lambda_B \sim \text{Gamma}(a_B, b_B)$ be independent.

We want $P(\lambda_A / \lambda_B < c)$.

**Step 1**: Standardize by rates.

Let $X = \lambda_A / b_A \sim \text{Gamma}(a_A, 1)$ and $Y = \lambda_B / b_B \sim \text{Gamma}(a_B, 1)$.

Then $\lambda_A / \lambda_B = (b_A / b_B) \cdot (X / Y)$.

**Step 2**: Use the Gamma-to-F relationship.

If $X \sim \text{Gamma}(a, 1)$ and $Y \sim \text{Gamma}(b, 1)$ are independent, then:

$$\frac{X/a}{Y/b} \sim F(2a, 2b)$$

So $X/Y = (a/b) \cdot F$ where $F \sim F(2a, 2b)$.

**Step 3**: Combine.

$$\frac{\lambda_A}{\lambda_B} = \frac{b_A}{b_B} \cdot \frac{X}{Y} = \frac{b_A}{b_B} \cdot \frac{a_A}{a_B} \cdot F$$

where $F \sim F(2a_A, 2a_B)$.

**Step 4**: Compute probability.

$$P\left(\frac{\lambda_A}{\lambda_B} < c\right) = P\left(F < c \cdot \frac{b_B}{b_A} \cdot \frac{a_B}{a_A}\right) = F_F\left(c \cdot \frac{b_B \cdot a_B}{b_A \cdot a_A}; 2a_A, 2a_B\right)$$

Wait, let me recheck this. We have:

$$\frac{\lambda_A}{\lambda_B} < c \iff \frac{b_A}{b_B} \cdot \frac{a_A}{a_B} \cdot F < c \iff F < c \cdot \frac{b_B}{b_A} \cdot \frac{a_B}{a_A}$$

Hmm, that's different from what I wrote earlier. Let me recalculate.

Actually, the correct relationship is:

If $X \sim \text{Gamma}(a, 1)$, then $2X \sim \chi^2(2a)$.

So if $\lambda_A \sim \text{Gamma}(a_A, b_A)$, then $2 b_A \lambda_A \sim \chi^2(2a_A)$.

The ratio:
$$\frac{\lambda_A / a_A}{\lambda_B / a_B} \cdot \frac{b_A}{b_B} = \frac{2 b_A \lambda_A / (2 a_A)}{2 b_B \lambda_B / (2 a_B)} \sim F(2a_A, 2a_B)$$

Therefore:
$$\frac{\lambda_A}{\lambda_B} = \frac{a_A}{a_B} \cdot \frac{b_B}{b_A} \cdot F$$

And:
$$P\left(\frac{\lambda_A}{\lambda_B} < c\right) = P\left(F < c \cdot \frac{a_B}{a_A} \cdot \frac{b_A}{b_B}\right)$$

So the correct formula is:

$$P(\text{HR}_{AB} < c) = F_F\left(c \cdot \frac{a_B \cdot b_A}{a_A \cdot b_B}; 2a_A, 2a_B\right)$$

For $c = 1$:

$$P(\text{HR}_{AB} < 1) = F_F\left(\frac{a_B \cdot b_A}{a_A \cdot b_B}; 2a_A, 2a_B\right)$$

**Verification in R**:
```r
# Verify with Monte Carlo
a_A <- 15; b_A <- 8
a_B <- 12; b_B <- 6

# Closed form
p_closed <- pf(
  q = (a_B * b_A) / (a_A * b_B),
  df1 = 2 * a_A,
  df2 = 2 * a_B
)

# Monte Carlo
set.seed(123)
lambda_A <- rgamma(1e6, shape = a_A, rate = b_A)
lambda_B <- rgamma(1e6, shape = a_B, rate = b_B)
p_mc <- mean(lambda_A / lambda_B < 1)

c(closed = p_closed, mc = p_mc)
# Should be approximately equal
```

### A.2 Corrected Implementation

```r
#' Compute P(lambda_j / lambda_l < c) for independent Gamma posteriors
#' 
#' @param a_j Shape parameter for numerator (arm j)
#' @param b_j Rate parameter for numerator (arm j)
#' @param a_l Shape parameter for denominator (arm l)
#' @param b_l Rate parameter for denominator (arm l)
#' @param c Threshold for hazard ratio (default 1)
#' @return Probability P(HR < c)
compute_p_hr_less_than_c <- function(a_j, b_j, a_l, b_l, c = 1) {
  pf(
    q = c * (a_l * b_j) / (a_j * b_l),
    df1 = 2 * a_j,
    df2 = 2 * a_l
  )
}
```

---

## Appendix B: Glossary

| Term | Definition |
|------|------------|
| **HR** | Hazard ratio |
| **PP** | Predictive probability |
| **BATON** | Bayesian Adaptive Trial Optimization Navigator |
| **OC** | Operating characteristic |
| **RAR** | Response-adaptive randomization |
| **Sufficient statistics** | $(n_k, T_k)$ = (number of events, total exposure) for arm $k$ |
| **Conjugate prior** | Gamma prior for exponential likelihood |
| **Gatekeeping** | Sequential testing where one hypothesis gates another |

---

*End of specification document.*
