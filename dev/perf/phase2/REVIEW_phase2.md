# Review of phase 2 (exact posteriors, counter-based RNG), commits 0854a58..3bb0694

Reviewer: independent Opus review, fresh context, 2026-11-04. Read-only. Files read from the commits (`git show 3bb0694:...`),
not the working tree. Checks run against `dev/perf/lib/dev_p2` (and `dev_l2` for the default-mode check), one core.

Counts: critical 0, high 1, medium 2, low 5. No defect found that would make a wrong result pass. The high finding is
about how gate G5 will read a correct engine, not a code defect.

## High

**H1. The exact engine is the infinite-draw limit of the phase 1 rule, and the phase 1 Monte Carlo noise biases
operating characteristics by about the size of the G5 power margin.** `src/simulate_chunk_exact.cpp:555-556, 576-578,
596-598`; plan Section 3 G5.
Phase 1 compares a Monte Carlo estimate (800 draws for R23, 5,000 for R01) with the threshold at every look. Noise gives
each look an extra chance to cross the futility threshold, so phase 1 stops more often for futility and has lower power
than the noise-free rule. This is not a defect, but it means a correct phase 2 engine differs systematically from the
reference rows. Confirmed by a run (`.simulate_chunk_exact_replay`, which is phase 1's loop with phase 1's MC posteriors,
with the draw counts set in cfg), R23 alternative, R = 3,000:

| Rule | Triplet power | Triplet futility stop | E[N] doublet / triplet |
|---|---|---|---|
| exact | 0.828 | 0.127 | 44.38 / 47.46 |
| MC 800 draws (phase 1 setting) | 0.813 | 0.134 | 44.18 / 47.29 |
| MC 100 draws | 0.784 | 0.153 | 44.08 / 46.98 |
| MC 30 draws | 0.765 | 0.176 | 43.45 / 46.31 |

The trend is monotone in the draw count and has the same sign as the P2-2 sanity pattern (exact higher power, fewer
futility stops, larger E[N], z 1.1 to 1.6 in both designs). Extrapolating, the 800-draw shift is on the order of 0.01 in
power, which is the G5 margin. So G5 can fail on R23-type designs with a correct engine, and the classification check
(power on the same side of 0.80) is exposed for reference rows near 0.80. Recommendation: before reading G5, decide how a
failure in this direction is attributed. The paired diagnostic the plan already prescribes can separate the two causes
if it also runs the replay instantiation at a large draw count (for example 20,000) on the same slices: an engine defect
shows up as exact vs MC-large disagreement, the expected shift only as MC-800 vs MC-large disagreement.

## Medium

**M1. No deterministic test ties the ExactPost wiring to the phase 1 decision rule.** `src/simulate_chunk_exact.cpp:515-557`
(policy), `:576, 596, 766-767, 787, 791` (call sites).
The replay identity (`tests/testthat/test-simulate-chunk-exact.R:48-60`) proves the template loop is phase 1's, but only
with `McPost`. Which threshold goes to which call (interim `futility_thr_interim` with `below`, final `success_thr` with
`above`, `futility_thr` with `below`; PH `ge` at the futility margins, `lt` at `log_eff_final`) is checked only by reading
and by the statistical sanity run. I read every call site against `src/simulate_chunk.cpp:375, 514, 741-742, 763, 767`
and found them correct, so this is a test gap. A swap (for example `ph_ge` and `ph_lt`, or the two final HC thresholds)
would pass the package tests and G3, and be caught only by G5 or the sanity table. Suggested test: a hook that records
each look's (slice, threshold, probability) from `ExactPost` and checks it against `.exact_pr_median_*` and `pnorm`
called directly, and against `McPost` with 10^5 draws within 4 standard errors.

**M2. No distributional test of the CounterRng samplers.** `src/simulate_chunk_exact.cpp:403-409`.
Tests cover the raw xoshiro outputs, the seeding and the uniforms (`test-simulate-chunk-exact.R:25-45`), but not the
exponential gap, the Weibull event time or the arm-assignment probabilities produced from them. I checked by reading: the
gap is `-log(U) * (1 / rate)` (the scale phase 1 passes to `R::rexp`), the event time is `scale * (-log U)^(1 / shape)`
(R's `rweibull` formula), and arm assignment reuses the `ProbSampleNoReplace` port with `rT = U`, which selects arm j with
its normalized probability; skipping the draw when one arm is eligible does not change any probability. All correct. A
parameterization slip (shape and scale swapped, rate used as scale) would again be caught only by G5. A cheap hook
returning per-trial patient draws, compared with `pweibull` and `pexp` by KS, would close this.

## Low

**L1. Quadrature non-convergence is silent in the engine.** `src/posterior_exact.cpp:246-268`, `src/simulate_chunk_exact.cpp:522-529`.
`contour_tail` returns its value when the 2,000-panel budget runs out, with `converged = false` only in `TailInfo`, and
`ExactPost` passes no `TailInfo`. Not observed: a stress run of 20,000 random cases (2 to 8 terms, shapes 10^-4 to 5,000,
scales 10^-8 to 10^3, threshold log 2) gave 0 non-converged, 0 non-finite, error estimate never above 1e-11, rounding
estimate never above 1e-10, at most 735 integrand evaluations. Extreme cases checked against independent references: shapes
1e-6 to 0.05 (singular densities), shapes 5,000 / 3,000 with tails of 1e-98 and 1 - 1e-14, scale ratios 1e-16 and 1e12,
equal and nearly equal scales; all agree to 1e-12 or better (one apparent 1.1e-6 discrepancy was my reference's error,
resolved by conditioning on the other gamma). Suggest a hard stop (or at least a counted warning) when `converged` is false,
so a future input outside the tested region cannot pass silently.

**L2. The engine option is silently ignored outside `run_simulation_pure`.** `R/simulate_chunk_exact.R:19-25`,
`R/simulation_driver.R:969-974`.
With `use_cpp_engine = FALSE` "exact" runs the R body (documented and tested), and other simulators (the hybrid path) never
read the option. A caller who sets `evolveTrial.engine = "exact"` gets phase 1 there with no message. Suggest a one-time
message or a note in the option's documentation.

**L3. The PH error conditions of phase 1 are not evaluated.** `src/simulate_chunk_exact.cpp:540-554` vs
`src/simulate_chunk.cpp:455-472`.
Phase 1 draws the gammas and evaluates two medians per draw only so that a NaN stops the run. The exact engine skips this,
so it would continue where phase 1 stops. Reachable only with a non-positive posterior rate or an infinite `exp(logHR)`;
the supported set has positive prior rates. Note only.

**L4. PSOCK workers load the package from their own library paths.** `R/simulation_driver.R` (cluster setup, unchanged),
`R/simulate_chunk_exact.R:63`.
When the parent loads evolveTrial with `lib.loc = dev/perf/lib/dev_p2` and the worker library paths do not include it, the
workers load whichever build they find. Without an exact-capable build installed elsewhere this fails loudly (no
`.simulate_chunk_exact`), so it cannot pass wrongly today; it could once a different exact build is installed. Same
exposure as phase 1. The gate harness should set `R_LIBS` for workers (it may already; I did not check `p2_cert.R`).

**L5. Cosmetic.** `ExactPost::warned` warns once per block, phase 1 once per call (unreachable in the set,
`simulate_chunk_exact.cpp:524-527`). `contour_tail` can return `-0` (`posterior_exact.cpp:256`); harmless in every
comparison.

## Checked and found correct

- Trial loop: line-by-line comparison of `run_trial` (`simulate_chunk_exact.cpp:620-805`) with `simulate_chunk.cpp:595-777`
  (arrival, phantom-arrival rules, rebalancing, look schedule, post-enrollment looks, final time, vs-ref and HC finals, NA
  handling, final N). Identical apart from the two policies.
- Single-arm reduction: `med < thr` iff cumulative hazard at `thr` exceeds log 2 for the piecewise median of
  `posterior_sampling.cpp:98-133`; weights `min(thr, end) - start`, open-ended last interval, threshold on a cutpoint gives
  the full interval and zero to the next, interval starts by repeated addition as phase 1; shape-0 terms are identically
  zero as `R::rgamma(0, .)` is; no atom at `med == thr`, so `above = 1 - below`. NaN threshold gives 0 / 0 as phase 1's
  tally does. Strict versus non-strict comparisons do not matter (continuous) except the `>=` against the decision
  thresholds, which is kept.
- Contour method: residue bookkeeping for the left vertex (`I = Pr - 1`), saddle brackets, tail bound on `[Y, inf)`
  (each factor bound holds), branch cuts never crossed. `n == 1` uses `pgamma(c, shape, scale, upper)` correctly.
- Between-arm: `ge` is the upper `pnorm` at the futility margins, `lt` the lower at `log_eff_final`, R's `rnorm` edge cases
  mirrored.
- Generator: seed64 from two 31-bit draws (62 bits), mix64 is a bijection so trial keys are distinct under one seed,
  splitmix64 seeding of xoshiro256++, U in [2^-53, 1 - 2^-53] (never 0 or 1, and the top value is representable).
  Cross-scenario collisions need a 64-bit coincidence.
- Driver and reducer: dynamic dispatch (`clusterApplyLB`), reducer sums in trial-index order with phase 1 accumulation,
  checks every index once, raw vectors at position = index.
- Default mode: confirmed by a run that `run_scenarios` on R01 (R = 200, sequential, seed 11) gives an identical table and
  identical `.Random.seed` with `dev_l2` and `dev_p2`. `.exact_engine_selected()` only reads options.

## Reports versus tests

- P2-2's "replay identity on all 20 scenario frames of the 10 designs" is run evidence; the package tests cover R01 and
  R23 at 120 trials. Accurate as stated.
- P2-2's statement that the replay "checks that the copied trial logic is the phase 1 logic" is accurate. Its scope (it
  cannot check the policies) is stated in the report and is the gap in M1 and M2.
- P2-1's G1 claims match the tests and my own reference checks above.
