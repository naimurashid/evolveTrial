# Hybrid and seamless engines: how the phase 1 and phase 2 lessons transfer

Date 2026-11-04. Read-only research note for the PI. No code was edited and nothing was committed. Every number is
marked (measured) or (estimated). Timing scripts are in the session scratchpad, not in the repo.

## 1. What the two families are in the code

The seamless designs of the paper (Scenario B, `cohortB_hybrid_*`, 13 searched parameters) and the hybrid designs run
on one engine. That engine is `src/hybrid_trial_sim.cpp`. `compute_hybrid_oc_cpp()` (line 1385) calls
`run_hybrid_simulations_cpp()` (line 1321), which calls `simulate_hybrid_trial_cpp()` (line 652) once per trial. The R
entry point is `compute_hybrid_oc_rcpp()` in `R/hybrid_sim_rcpp.R`, and the paper repo reaches it through
`sim_fun_hybrid()` in `R/warmstart_wrappers.R` (lines 1386 to 1440). The engine is already C++, so there is no R trial
loop to port. The design document's phrase "seamless C++ engine" refers to this file.

The engine has four modes (line 38). The seamless designs use `hybrid` mode (single-arm monitoring, a
predictive-probability conversion decision, then between-arm monitoring) with posterior single-arm decisions. The other
modes (`single_arm`, `between_arm`, `dual_single_arm`) and the predictive-probability stopping options
(`efficacy_method` or `futility_method` set to `predictive`) are the other family. I call them "hybrid variants" below.
None of the production seamless designs uses them. A pure-R reference engine (`R/hybrid_trial.R`, `R/hybrid_ssr.R`)
exists as a fallback in `sim_fun_hybrid()`. It is not on the production path and is not worth porting.

How an evaluation is simulated. Each call runs three scenarios of R trials each (alternative, global null, and the
between-arm null that gives `type1_between`). Event times are piecewise exponential with constant hazards, so the
Weibull shape in the scenario CSV is not used unless the M1 switch is injected (atbo-longleaf
`pilot/verify/R/verify_m1_switch.R`). Looks fall every 3 months with a fixed number of patients per look. The wrapper
sets `n_outer = 1000` for the conversion predictive probability. With `parallel_replicates`, the R wrapper splits the
trials across `mclapply` workers with seeds drawn by `sample.int` after `set.seed`, so results depend on the worker
count, as phase 1 results did.

## 2. Where the time goes (profiled)

`perf` is blocked on this machine (`perf_event_paranoid` 4), so I profiled by switching components off. The design was
`cohortB_hybrid_balanced_bidir` (stored best vector, `n_intervals` 3), library `dev_p2_final`, one worker, load about 2.
Its operating characteristics matched the stored row (power 0.81, `EN_null` 85.9, `EN_alt` 92.1 at R = 1,000, measured).

| Run (1 worker, seconds per 1,000 trials per scenario, 3 scenarios) | Time |
|---|---|
| Production settings, `n_outer` 1,000 | 14.0 (measured) |
| `n_outer` 200 | 11.7 (measured) |
| `n_outer` 2 (conversion predictive probability nearly free) | 10.4 (measured) |
| `n_intervals` 1, `n_outer` 2 (single-arm posterior becomes one `pgamma`) | 0.06 (measured) |
| `n_intervals` 8, `n_outer` 2 | 28.8 (measured) |

The time splits as follows at the production design.

| Bucket | Share |
|---|---|
| Single-arm posterior by Monte Carlo (`compute_p_single_arm_internal`, 2,000 gamma draws per interval and a piecewise median per draw, at every look) | about 74% (measured by difference) |
| Conversion predictive probability (`compute_pp_predictive_internal`, up to 1,000 antithetic outer draws, each with two `qgamma` per interval, a simulated future cohort for both arms and a Newton solve) | about 26% (measured by difference) |
| Trial loop, patient draws, look slicing, between-arm Laplace posterior | under 1% (measured, 0.06 s) |

The `n_intervals` 1 run changes the design a little (power 0.80 against 0.81, measured), so the last row is a proxy
for the loop cost, not an exact share. The single-arm share grows with `n_intervals`, which BATON searches.

Hybrid variants (R = 300, `n_intervals` 3 unless stated, measured). `single_arm` mode with posterior decisions costs
10.3 s per 1,000 trials per scenario, almost all of it the same 2,000-draw posterior. Predictive-probability stopping
costs 45.6 s in `single_arm` mode and 49.6 s in `hybrid` mode. At `n_intervals` 1, predictive single-arm stopping still
costs 19.6 s, so the single-arm predictive probability (`compute_pp_efficacy_sa_cpp`, `compute_pp_futility_sa_cpp` in
`src/predictive_probability.cpp`) dominates there. It runs 1,000 outer draws at every look after the event gate, with no
early stop.

Parallel behavior today (production design, R = 3,000, measured). 16 workers took 3.85 s against an estimated 42 s at
one worker, an efficiency of about 0.68. At `n_intervals` 1, a stand-in for an engine without the single-arm draws,
efficiency was 0.63 (7.10 s at one worker, 0.70 s at 16). The 1-worker and 16-worker runs gave different type I error
estimates from the same seed, which confirms the worker-count dependence.

## 3. What transfers

**Exact single-arm posterior (phase 2, highest payoff).** The hybrid single-arm decision is P(median of the arm's
piecewise exponential survival exceeds m*), with m* the historical median divided by the hazard-ratio threshold. The
median exceeds m* exactly when the cumulative hazard at m* is below log 2. That is the weighted-sum-of-gamma-hazards
tail that phase 2 already computes (`gamma_sum_upper_tail` in `src/posterior_exact.h`), with weights equal to the
exposure of each interval up to m*. For the production seamless designs m* is 6 months and every interval is at least
12 months wide (follow-up 120 months, at most 10 intervals), so the exact answer is a single `pgamma` on the first
interval. The engine now spends three quarters of its time estimating that number with 2,000 draws. Replacing it removes
about 10.3 of 14.0 s (estimated from the measured split), a 3.8x gain at `n_intervals` 3 and about 8x at 8 (estimated).
The code change is a few lines, because the function already exists.

**RNG streams (phase 2).** Per-trial xoshiro256++ streams seeded from (scenario seed, trial index), with sums in
trial-index order, make results independent of the worker count and allow dynamic dispatch. They also give common
random numbers across candidate designs for free (Section 6).

**R and Rcpp bookkeeping (phase 1 spirit).** The trial loop is C++, but it is not lean. Each trial re-reads its
parameters from Rcpp lists and compares strings (`conversion_trigger`, `futility_action`), returns a named `List`, and
the driver converts the result field by field. Inside the predictive-probability loop, `simulate_future_arm_internal`
returns an Rcpp `List` and `ph_beta_mode_var_cpp` returns one, so every outer draw allocates R objects. Replacing them
with plain structs and parsing the parameters once per call keeps R's random stream, so it can be bit-identical. It also
removes the R API calls that block threading. The payoff is modest, an estimated 1.2 to 2x on the predictive-probability
part only.

**Faster draws inside the remaining Monte Carlo (phase 2 future extension).** The predictive probability keeps its Monte
Carlo. Its antithetic outer draws use `qgamma`, which costs about 0.53 s per million calls against 0.03 s for `rgamma`
in R (measured, vectorized). Direct gamma draws (Marsaglia-Tsang on xoshiro) or a cheaper inverse would cut the draw
cost about tenfold (estimated), and the future-cohort simulation would gain from the faster uniform. Combined with the
struct change, the predictive probability should run 3 to 10x faster (estimated).

## 4. What does not transfer, and why

- **A closed form for the predictive probability.** Phase 2 replaced draws by the infinite-draw limit of the same rule.
  The predictive probability's infinite-draw limit is an integral over future event counts and exposures by interval and
  arm, with no closed form. A normal approximation to the predictive distribution of the log hazard ratio gives a
  closed form, but it is a different decision rule, not a faster computation of the current one. It would change the
  designs and needs a separate PI ruling. Quasi-Monte Carlo (Sobol) outer draws are a middle path that keeps the rule
  and lowers noise per draw, but it is research effort for a 26% bucket.
- **Look slicing.** `compute_interval_metrics` rescans every patient at every look (and twice for the evaluated arm),
  but the whole loop costs under 1% (measured). Incremental bookkeeping would also change the summation order, so it
  could not sit in a bit-identical layer. It is not worth doing.
- **The between-arm posterior.** It is already the closed-form Laplace normal tail (`compute_ba_ph_posterior_internal`),
  so phase 2's between-arm gain is already in this engine.
- **A phase 1 style port.** There is no R loop to port. The pure-R hybrid engine is a fallback only.
- **The single-arm predictive probability (hybrid variants).** Its inner step is already a `pgamma` on pooled counts, so
  the exact posterior does not help it. Only the draw, struct and early-stopping changes apply. Its inner step pools
  the intervals into one rate, which is a different model from the piecewise-median posterior used for posterior
  decisions. That is a modeling note, not a speed item, but a gate on this path should not assume the two agree.

## 5. Estimated speed-up and verification

| Step | Kind | Seamless gain at 1 worker | Basis |
|---|---|---|---|
| Exact single-arm posterior | statistical | 3.8x at `n_intervals` 3, about 8x at 8 | estimated from the measured split |
| Structs and parse once, in the predictive probability and the trial wrapper | bit-identical | 1.2 to 2x on the 26% bucket | estimated |
| xoshiro streams, direct gamma draws in the predictive probability | statistical | 3 to 10x on the 26% bucket | estimated from the `qgamma` timing |
| All three | | 10 to 30x per evaluation (14.0 s to roughly 0.5 to 1.4 s per 1,000 trials per scenario) | estimated |
| Hybrid variants with predictive stopping | | 3 to 10x | estimated, same draw arguments |

This is well below phase 2's 190x for the single-arm and between-arm engine, because the predictive-probability Monte
Carlo stays. After the first step it is about 97% of the remaining time (estimated). Once it shrinks, the `mclapply`
start-up per call will exceed the work, as in phase 2's G6, so in-process C++ threads become necessary rather than more R
workers. Threads need the struct change first (no R API calls inside the loop) and per-trial streams.

Gates, in the phase 2 pattern.

1. **Bit-identical layer** for the struct and parse-once refactor only. `identical()` on the per-trial data frame from
   `run_hybrid_simulations_cpp` at fixed seeds, for the four seamless vectors and the hybrid variants, at `n_intervals`
   1, 3 and 8, at 1 and 16 workers. A bisection tool over trial index, as in phase 1, finds the first divergence.
2. **Exact probability check** (phase 2 G1). The exact single-arm probability against 10^7-draw Monte Carlo on random
   posterior shapes, rates and thresholds, including thresholds beyond the first cutpoint and beyond the last.
3. **Determinism** (phase 2 G3). Identical results at 1, 4, 16 and 19 workers and two block sizes.
4. **Statistical equivalence against fresh runs**, per the PI's phase 2 ruling, not against stored rows. Pooled
   R = 100,000 per design for the current engine at fresh seeds and for the new one, two one-sided tests on power, type I
   error, `type1_between` and conversion rate, and expected sample sizes within their margins. Phase 2 showed that the
   exact rule is the infinite-draw limit and that finite draws add stops. A third arm at a high draw count (for example
   20,000 single-arm draws) separates a draw-count shift from a defect.
5. **Speed**, with the known `per_1000_s` convention stated (each call runs three scenarios here).

## 6. Order by payoff per effort

1. Exact single-arm posterior in `hybrid_trial_sim.cpp`, reusing `gamma_sum_upper_tail`. Small change, largest gain,
   needs gates 2 and 4.
2. Structs and parse-once in the trial wrapper and the predictive probability. Low effort, bit-identical (gate 1), and
   the prerequisite for threads.
3. Per-trial xoshiro streams with in-process threads in C++. Moderate effort. Removes the worker-count dependence and
   the start-up cost (gates 3 and 4).
4. Direct gamma draws in the predictive probability, and an early stop for the single-arm predictive probability like
   the conversion one already has. Moderate effort, statistical, gate 4.
5. Quasi-Monte Carlo for the predictive probability. Research effort, deferred.
6. Incremental look slicing. Not worth doing.

## 7. How BATON could use the faster engine

- **More replications per evaluation.** At an estimated 0.5 to 1.4 s per 1,000 trials per scenario, an R = 10,000
  seamless evaluation costs about 15 to 40 core-seconds (estimated) against about 420 today (estimated from the measured
  14.0 s). Searching at R = 10,000 instead of 3,000 cuts the Monte Carlo standard error of a type I error near 0.08 from
  about 0.005 to 0.0027 (computed), which tightens the feasibility boundary the surrogate has to learn. The separate
  top-candidate re-screen at 10,000 could then be dropped.
- **Larger budgets and screens.** The measured seamless feasibility screen cost 19.96 core-hours (atbo-longleaf
  `docs/reviews/jasa_round1/2026-10-26-reopen-analyses-plan.md`, line 401). At 10 to 30x it would cost about 0.7 to
  2 core-hours (estimated), so the 300-point screen could become a much larger initial design.
- **In-process batching.** BATON proposes a batch of three candidates per iteration. One C++ call that runs all
  candidates and all trials on one thread pool avoids a fork per candidate. Per-trial streams then give common random
  numbers across candidates, so differences between designs carry less noise, which helps both the surrogate and the
  final comparison.
- **The surrogate becomes visible.** GP fitting and acquisition took 2 to 4% of iteration wall time
  (`docs/reviews/jasa_round1/t1-timing-preflight.md`, verdict item 2, measured). With a 10 to 30x faster simulator that
  share rises to perhaps 20 to 50% (estimated), so the design document's decision not to pursue faster surrogates
  (Section 5) should be revisited once the engine lands.
- **`n_outer` is a measurement setting.** Like the old posterior draw count, it should not become a searched
  coordinate. If it stays finite, it belongs with R as a fidelity setting (lower in search, higher at certification).

## 8. Side observations for the gates

- `ss_method`, `n_add_candidates` and `max_additional_n` are passed to the C++ engine but not read there. The C++ engine
  sets the added sample size to `nmax_ba` minus the enrolled count. Sample-size re-estimation exists only in the R
  engine.
- Event times ignore `weibull_shape` without the M1 switch (known, PI ruling 2026-10-21). Any gate run should state
  which mode it used.
