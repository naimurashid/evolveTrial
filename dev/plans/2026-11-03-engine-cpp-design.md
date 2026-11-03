# evolveTrial engine speed-up: design

Date: 2026-11-03. Approved by the PI section by section in the session of 2026-11-03.
Branch `perf/engine-cpp`, cut from 58a4850 (tag `jasa-engine-p6-nodropout`, the engine of record for the
JASA single-arm and between-arm designs). Evidence files are in `dev/perf/evidence/`.

## 1. Scope, baseline and targets

The goal is to make the single-arm and between-arm simulator as fast as possible. BATON's surrogate is out of
scope, because fitting and acquisition take 2 to 4% of a search iteration's wall time (atbo-longleaf
`docs/reviews/jasa_round1/t1-timing-preflight.md`, lines 16 and 100 to 104). The seamless C++ engine is also
out of scope, as are the frozen JASA results and the reproducibility package. Main does not contain 58a4850,
so reconciling this branch with main is a separate, later step.

Baseline on the engine of record, P6 options on, dropout off (`dev/perf/evidence/profile_report.md`).

| Design | s per 1,000 trials, 1 worker | Share of time in R code a port can replace | 16-worker efficiency |
|---|---|---|---|
| R01, single-arm, Scenario A | 13.4 | 36% | 0.64 |
| R23, between-arm, Scenario B | 75.1 | 99% | 0.56 |

Where the time goes at one worker.

| Bucket | R01 | R23 |
|---|---|---|
| Posterior sampling already in C++ (gamma draws) | 64% | 0% |
| Posterior sampling in R (median solver loop, rgamma and rnorm, sampler loop) | 1% | 67% |
| Data slicing at each look | 17% | 20% |
| Trial loop bookkeeping, data frames, gates | 17% | 12% |
| Patient-level random draws | 0.4% | 0.5% |

The between-arm posterior runs in R because `posterior_helpers.R:105` redefines `sample_vs_ref_medians_ph`
after `posterior_cpp_dispatchers.R:241` (files are sourced alphabetically and DESCRIPTION has no Collate
field). The C++ version it shadows uses `arma::randn`, which does not match R's `rnorm`. Swapping it in moved
R23 power from 0.79 to 0.75 at R = 300, so it is dead code and must not be reused.

The work runs in two phases.

- **Phase 1, bit-identical.** Port the body of `simulate_chunk` to C++. Every reported number of the engine of
  record must reproduce exactly. Target is 10x or better for between-arm designs at one worker (75 s to under
  8 s per 1,000 trials) and about 1.5x for single-arm designs, where 64% of the time is the fixed count of
  gamma draws that bit-identity forbids changing. Parallel efficiency does not improve in this phase, because
  chunk boundaries are part of the seeding.
- **Phase 2, statistically equivalent.** Exact posterior probabilities by one-dimensional integration, which
  removes the gamma draws and the Monte Carlo noise in each decision. Per-trial counter-based seeding, so
  results stop depending on the worker count and chunks can be balanced dynamically. Incremental event
  bookkeeping in place of re-slicing the registry at every look.

## 2. Architecture of the bit-identical port

One new Rcpp entry point, `simulate_chunk_cpp(sim_indices, args, seed)`, replaces the body of
`simulate_chunk` and returns the same plain list (`sum_final_n`, `sum_final_events`, `sum_stop_eff`,
`sum_stop_fut`, `sum_final_eff`, `sum_final_fut`, `sum_final_inc`, `sum_final_time`, `n_sims`, and the raw
per-trial vectors when `return_variance` is on). `aggregate_results`, `run_scenarios`, the cluster code and
all argument parsing stay in R, so callers see no change. An option `evolveTrial.use_cpp_engine`, default
TRUE, keeps the R body available for the gate.

The R driver is unchanged (`splitIndices`, the `sample.int` seed list, `set.seed` once per chunk, FORK
workers). R's random number generator is not thread-safe and each chunk owns one stream, so the C++ chunk
runs single-threaded with one `RNGScope` and no callbacks into R that draw.

Per trial, the C++ chunk does the following, in the order of `dev/perf/evidence/rng_driver_inventory.md` Section 1.

- Arrival loop with `R::rexp(1 / rate)` drawn before every break test. The look comes between the arrival
  draw and the arm draw. The `looks_stop_at_final` phantom-arrival rules are copied as written.
- Arm assignment by a port of R's `ProbSampleNoReplace` (revsort, then one `unif_rand`). It draws one uniform
  even when only one arm is eligible.
- `R::rweibull` and `R::runif` per patient. The uniform is drawn and then discarded when dropout is off.
- Registry as preallocated arrays per arm (enrollment time, event time, censoring time), capacity
  `ceiling(1.2 * Nmax)` as now.
- Look slicing, interval metrics, `findInterval` semantics and the rebalanced cutpoints (type 7 quantile once
  50 or more events are observed) in C++. Long-double accumulation wherever R uses `sum` or `mean`, `mean`
  with R's second correction pass, and `median` as the mean of the two middle values.
- Single-arm posterior through the existing `R::rgamma` loop (interval outer, sample inner).
- Between-arm posterior through a new C++ copy of the live R sampler (`R::rnorm`, then `R::rgamma` for each
  sample, in the same order). The vs-reference look draws `num_posterior_draws` (800), not the interim count,
  and the control-cache draws that are spent and discarded when more than one experimental arm is evaluated
  are kept.
- The hazard-ratio mode and variance keep the existing C++ Newton solver `ph_beta_mode_var_cpp`.

The package `Makevars` adds `-ffp-contract=off` and nothing else. Local R and Longleaf's R both compile with
plain `-O2` and no `-march`.

Order of work. The single-arm path is ported and gated first, then the between-arm path. Each lands as its own
commit, with the gate result in the message.

## 3. Verification gate and tests

Three layers, in order of strength. Nothing in phase 1 is committed with layer 1 or layer 2 failing.

1. **Chunk-level identity test** (testthat, seconds). For each of the 10 acceptance designs, run the R body
   and the C++ body on the same chunk indices and seed at R = 50 and require `identical()` on the returned
   list, raw vectors included. Repeat at three seeds and at 16 and 19 workers, so the chunk boundaries
   differ. On a mismatch, a bisection tool over `sim_indices` reports the first trial that diverges. This tool
   is built first, before any porting, because it is the main debugging instrument.
2. **Reference-row gate**, the phase 1 acceptance criterion. Reproduce all 75 stored certification rows
   (`dev/perf/evidence/acceptance_set_reference_rows.csv`, 10 designs, 5 or 10 seeds each, R = 10,000 per seed, 19
   workers, P6 options on, dropout off). Pass means `power_count` and `type1_count` are equal as integers and
   `EN_null` and `EN_alt` are equal as doubles. The rows were cross-checked against the committed Tier 1
   readout with no mismatches.
3. **Search-side check.** Two or three E2c BATON evaluations from the stored ledgers, reproduced at R = 3,000
   and 16 workers from their recorded seeds, to cover the 16-chunk case that certification does not.

Each floating-point reduction (long-double `sum` and `mean`, `median`, type 7 quantile, `findInterval`) gets
a unit test against the R function on random inputs, since a one-ulp difference there flips a gate and shifts
the random stream for the rest of the chunk.

Phase 2 reuses the cloud prototype's equivalence gate (atbo-longleaf `pilot/cloud/compare_equiv.R`) in its
statistical mode. That means pooled R = 100,000 per design, two one-sided tests on power and type I error
with bounds at the certification margins, and expected sample size within its Monte Carlo standard error.

## 4. Rollout, compute and risks

Work happens in the worktree `~/research/evolveTrial-perf`, so the installed engine and the frozen JASA
stack stay untouched. Implementation follows the lab delegation practice. An Opus agent does the C++ port,
a Sonnet agent builds the bisection tool and the floating-point unit tests from a written specification, and
the coordinator runs the gate itself and makes every commit.

Estimated compute for the gate (only the new engine runs, because the reference rows already exist).

| Layer | Trials | Cost after the port (estimate) |
|---|---|---|
| 1, chunk identity | about 6,000 | minutes, local |
| 2, reference rows | 1.5 million | 3 to 5 core-hours, Longleaf |
| 3, E2c search evaluations | about 20,000 | under 1 core-hour |

Risks.

1. **R version.** The reference rows came from Longleaf's R 4.5.0 and this machine runs R 4.6.0. Layer 1
   compares the two bodies on one machine, so it does not depend on the version. Layer 2 runs on Longleaf
   with the original modules (gcc/12.2.0, r/4.5.0), so a libR difference cannot pass for a port defect.
2. **A one-ulp floating-point difference flips a gate.** The bisection tool finds the first diverging trial,
   and the unit tests check each reduction against R.
3. **Posterior draw counts differ by path.** The between-arm look draws 800 samples and the single-arm look
   draws 5,000. This is undocumented, and phase 1 keeps it exactly. Phase 2's exact integration removes it.
4. **The JASA submission.** Nothing on this branch feeds the frozen text, the tables or the reproducibility
   package, and it is not cited in the resubmission.

## 5. Decisions recorded in the design session

- Equivalence standard. Bit-identical first, then statistically equivalent (PI choice of the recommended
  option).
- The number of posterior draws is not a searched BATON coordinate. The objective and constraints depend on
  it only through noise, so an optimizer would learn to pick values whose noise happens to satisfy a
  constraint. It is a measurement setting. If it stays finite after phase 2, it belongs with R as a fidelity
  setting (low in search, high at certification).
- BoTorch or GPyTorch are not pursued for speed, because the surrogate is 2 to 4% of wall time.

## 6. Future extensions (PI 2026-11-03, keep in mind, not scheduled)

These apply only once differences from the engine of record are allowed, and each needs the phase 2
statistical gate.

- **Faster generator family.** A counter-based or xoshiro stream with a ziggurat normal and Marsaglia-Tsang
  gammas, in place of R's Mersenne Twister and its gamma method. Estimated 3 to 5x on the draw portion only.
  It is a fallback for any path that still samples after exact integration, since integration removes the
  draws altogether.
- **The dead C++ hazard-ratio sampler** (`sample_vs_ref_medians_ph_cpp`, `arma::randn`). Its 2.5x speed-up
  came from compiled code, not from its normal generator, so phase 1 already captures that gain. The R23
  power shift it showed (0.79 to 0.75 at R = 300) is within noise for a difference of two runs on different
  streams (standard error about 0.033), but that is not evidence of equivalence. It has never been tested
  against the live R sampler.
- **Posterior draw count as a fidelity setting.** If any finite draw count survives phase 2, it sits with R
  (low in search, high at certification), never as a searched design coordinate.
