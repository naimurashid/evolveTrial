# evolveTrial 58a4850: RNG inventory, seeding, chunking, and C++ port feasibility

Source: `~/research/evolveTrial` at 58a4850 (branch verify/p6-nodropout), read through a temporary
detached worktree (removed after the read). File:line references are to that commit. Empirical checks
were run in plain R 4.6.0 plus one scratch `Rcpp::sourceCpp` file (no package install). The script and
C++ file are in `/tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/537dc1ca-2d2f-4a33-9259-7d770244bb0a/scratchpad/perf/rngcheck/`.

Scope note. Inside `run_simulation_pure`, the "between-arm" path is `compare_arms_option = TRUE`, the
vs-reference comparison. The hybrid SA-to-BA engine (`compute_hybrid_oc_cpp`, `src/hybrid_trial_sim.cpp`)
is already C++ and is out of scope here. `R/design_analysis.R` has no RNG calls. It only calls
`run_scenarios(..., seed = ...)` (lines 208, 290, 375).

Frozen configuration that decides which branches run (F0 repo commit 1863cecb, `R/warmstart_wrappers.R`):

| Setting | Single-arm (lines 380-480) | Multi-arm (lines 742-842) |
|---|---|---|
| `compare_arms_option` | FALSE | TRUE (2 arms, reference plus 1 experimental) |
| `use_ph_model_vs_ref` | FALSE | TRUE |
| `num_posterior_draws` (interim not set, so equal) | 5000 | 800 |
| `randomization_probs` | 1.0 | 0.5 / 0.5 |
| `censor_max_time_sim` | `fu_time` (finite) | `fu_time` (finite) |
| `rebalance_after_events` | 50 | 50 |
| `person_time_milestones` | NULL (calendar beats) | NULL |
| Cluster | `num_workers = cached_num_workers` (SLURM_CPUS_PER_TASK), `cluster_type = "auto"` (FORK on Linux), no external cluster | same |

## 1. RNG inventory per simulated trial, in execution order

No file calls `RNGkind`, `.Random.seed`, L'Ecuyer streams or `clusterSetRNGStream` on this path. That
holds for evolveTrial at 58a4850, F0 1863cecb `R/` and `scripts/` (the only `.Random.seed` uses are in
the unrelated `biomarker_sim_harness.R` and `bo_parameter_importance.R`), BATON v0.4.0 `R/`, and the
pilot `R/`, `p6recert/R` and `verify/R` directories. The pilot's `opt_common.R:299-308` saves and
restores `.Random.seed` around an eval, outside the simulator. The kind is therefore the default
Mersenne-Twister / Inversion / Rejection, which is what R reports locally.

### 1a. Accrual loop (one iteration per arrival), `R/simulation_driver.R:542-679`

| Order | Call | file:line | Length | Loop | Notes |
|---|---|---|---|---|---|
| A1 | `rexp(1, rate = overall_accrual_rate)` | simulation_driver.R:543 | 1 | yes, every arrival | Drawn before every break test. The arrival that ends accrual (the `looks_stop_at_final` cap at 552-555, the phantom vs-ref break at 561-569, and the breaks at 635/638/645) has already used its draw. `exp_rand()` uses a variable number of uniforms internally. |
| (A1b) | rebalance check `extract_events_fast` | 501-520, 572-603 | 0 | yes | No RNG. It moves cutpoints (deterministic `quantile` type 7, line 1337) once there are 50 or more events, and every later posterior draw count and gate depends on this. |
| A2 | interim look, if `current_time >= next_calendar_look` | 629-632 -> `interim_check` | see 1b | at most one look per arrival | Posterior draws happen here, *between* A1 and A3 of the same iteration. |
| A3 | `sample(elig_arms, size = 1, prob = probs)` | 652 | 1 | yes, per enrolled patient | **sample() flag.** It goes through `sample.int(length(x), 1, FALSE, prob)`, then `.Internal(sample)`, then `ProbSampleNoReplace`: `revsort(p, perm, n)`, then `rT = totalmass * unif_rand()`, which is exactly one `unif_rand` and does *not* use `R_unif_index`, so `sample.kind` does not matter. When there is a single eligible arm (single-arm, or multi-arm after one arm stops or fills) a uniform is still used (checked: `sample("E", 1, prob = 1)` advances the stream by one uniform). With equal probabilities, `revsort` puts the second arm first, so u <= 0.5 picks arm 2 and u > 0.5 picks arm 1 (checked: u = 0.989 chose "D"). |
| A4 | `rweibull(1, shape, scale)` | 655-656 | 1 | yes, per enrolled patient | `scale * pow(-log(unif_rand()), 1/shape)` in libR. |
| A5 | `runif(1, 0, censor_max_time_sim)` | 657 | 1 | yes, per enrolled patient | Rejection loop `while (u <= 0 \|\| u >= 1)` in libR. It uses **zero** uniforms if the max is non-finite (checked: `runif(1, 0, Inf)` returns NaN and leaves the stream unchanged) or if min == max. Frozen `fu_time` is finite, so one uniform. Line 658 then overwrites the value when `random_dropout = FALSE` (Section 5). |

### 1b. Interim look, single-arm (HC), `R/interim_logic.R:171-263`

For each arm still recruiting (172-173), the slice (no RNG, `state_management.R:115-149`) is followed by
the gates at 218-223 (n, events, `stats::median` follow-up, person-time fraction). If the gates pass:

| Order | Call | file:line | Length | Loop |
|---|---|---|---|---|
| H1 | `draw_posterior_hazard_samples` -> `draw_posterior_hazard_samples_cpp` | interim_logic.R:47-54 -> posterior_cpp_dispatchers.R:37-46 -> posterior_sampling.cpp:57-65 | K x S = n_intervals x `num_posterior_draws_interim` (frozen 8 x 5000 = 40,000 by default) | C++ double loop: interval j outer, sample i inner, `R::rgamma(alpha_j, 1/beta_j)` |

The order (column by column) is the same as the R fallback `rgamma(num_samples, ...)` per column
(posterior_cpp_dispatchers.R:69-79). `R::rgamma` uses R's RNG through the `Rcpp::RNGScope` that the
generated wrapper creates (`src/RcppExports.cpp:61`), which calls GetRNGstate and PutRNGstate. Gamma
sampling uses a data-dependent number of `norm_rand`/`unif_rand`/`exp_rand` draws (rejection), but the same
libR function runs either way.

### 1c. Interim look, between-arm (vs-ref), `R/interim_logic.R:84-169`

| Order | Call | file:line | Length | Condition |
|---|---|---|---|---|
| V0 | `precompute_ctrl_posteriors` -> C++ gamma loop | interim_logic.R:116-119 -> posterior_helpers.R:79-91 | K x `num_posterior_draws` | only when `!use_ph_model_vs_ref && length(eval_arms) > 1`. **The result is never used**: neither the C++ dispatcher (posterior_cpp_dispatchers.R:171-191) nor the R fallback reads `ctrl_cache`, so these draws are spent and thrown away. Not reached in the frozen 2-arm PH configuration. |
| per experimental arm that is recruiting and passes `gates_pass_for_both_arms` (127): | | | | |
| V1 (PH, frozen) | `rnorm(num_samples, mean, sd)` | posterior_helpers.R:124 | S = **`args$num_posterior_draws`** (800), *not* the interim count (interim_logic.R:7) | once |
| V2 (PH, frozen) | `rgamma(length(shapes), shape = shapes, rate = rates)` | posterior_helpers.R:130 | K per sample | R-level `for (i in seq_len(num_samples))` loop, posterior_helpers.R:128-133 |
| V1' (independent model) | `sample_vs_ref_medians_independent_cpp`: lamC (K x S) then lamT (K x S) | posterior_cpp_dispatchers.R:182-191 -> posterior_sampling.cpp:357-381 | 2 x K x S | not frozen |

**Shadowing (important).** `sample_vs_ref_medians_ph` is defined twice:
`posterior_cpp_dispatchers.R:241` (it dispatches to `sample_vs_ref_medians_ph_cpp`) and
`posterior_helpers.R:105` (a pure-R loop). DESCRIPTION has no `Collate:` field, so the files are sourced
in alphabetical order and `posterior_helpers.R` comes later and wins. The installed namespace confirms it:
`deparse(evolveTrial:::sample_vs_ref_medians_ph)` shows the `coalesce_num` / pure-R version. The PH draws
are therefore R `rnorm` plus R `rgamma`. The only C++ call in that path is the deterministic Newton
solver `ph_beta_mode_var_cpp` (posterior_sampling.cpp:163-222). `sample_vs_ref_medians_ph_cpp` is dead
code, and it would **not** match: it uses `arma::randn` (posterior_sampling.cpp:301), which under
RcppArmadillo's `ARMA_RNG_ALT` (`RcppArmadillo/rng/Alt_R_RNG.h:82-101`) is a Marsaglia polar method built
on `Rf_runif`, not R's inversion `norm_rand`. Checked: `arma::randn(800) * sd + m` is not identical to
`rnorm(800, m, sd)`.

### 1d. Post-enrollment looks and final analysis

- Post-enrollment looks: simulation_driver.R:681-696 (`looks_stop_at_final = FALSE`: one look at the end
  of accrual, then a look every beat until `max_trial_time`) or 697-716 (TRUE, the P6/D1 setting: looks
  at `last_enroll + j*beat` strictly before T\*). Each look draws as in 1b or 1c. There is no other RNG.
- Final, single-arm: simulation_driver.R:759-766, the same C++ gamma loop, K x `num_posterior_draws`, for
  each arm still recruiting (751).
- Final, vs-ref: simulation_driver.R:798-803 -> `sample_vs_ref_medians` (posterior_helpers.R:93-103)
  -> V1/V2 (PH), for each recruiting non-reference arm. `ctrl_cache` is NULL.

### 1e. Not on the path

- `R/predictive_probabilities.R` (rgamma/rpois at 281-285, 367-368) and its C++ counterparts: none of the
  `calculate_predicted_*` functions has a caller in `R/`. The `pred_*` and `num_posterior_draws_pred`
  arguments are packed into `args` (408-413) but never used.
- `src/predictive_probability.cpp`: only `compute_pp_predictive_cpp` is called, from `hybrid_ssr.R:231`
  (hybrid path).
- `data_generation.R:86,125` (`simulate_piecewise_exponential_data`): called only from
  predictive_probabilities.R.
- `hybrid_sim_rcpp.R:57-58`: the hybrid engine has its own scheme (`worker_seeds <- sample.int(...)`
  after `set.seed`). Not this driver.

### 1f. Count

Per trial there are **6 RNG call sites** on each frozen path. These are 4 accrual sites (A1 rexp, A3
sample, A4 rweibull, A5 runif) and 2 posterior sites. Single-arm: the posterior sites are the C++
`R::rgamma` loop reached from the HC interim and the HC final. Multi-arm PH: R `rnorm` and R `rgamma`
reached from the vs-ref interim and the final. The independent vs-ref model would add V0 and the two
C++ loops instead. Per driver call (not per trial) there are 2 more: `sample.int` (line 1036) and
`set.seed` (line 470).

Every count is data-dependent. Arrivals = enrolled + 1 (the last arrival breaks). A3-A5 run once per
enrolled patient. Each look draws only for arms still recruiting that pass the gates, and stopping
removes later looks. The final analysis draws only for arms still recruiting. Rebalancing moves the
cutpoints, which moves the gates. The gamma and exponential rejection samplers and the runif boundary
loop add variable internal draw counts (these are inside libR and match automatically).

## 2. Seeding and chunking

`run_simulation_pure` has no `seed` argument and never calls `set.seed` on its own stream. Seeds come
from the caller's global stream:

```r
# run_scenarios, simulation_driver.R:1272
if (!is.null(seed)) set.seed(seed)
# run_simulation_pure, simulation_driver.R:976-977, 1036, 1042-1043
chunks <- parallel::splitIndices(num_simulations, workers)
chunks <- chunks[lengths(chunks) > 0]
seed_list <- replicate(length(chunks), sample.int(.Machine$integer.max, 1L), simplify = TRUE)
  function(idx, chunk_indices, seed_vals) {
    simulate_chunk(chunk_indices[[idx]], seed = seed_vals[[idx]], tick = function() {})
# simulate_chunk, simulation_driver.R:469-471
if (!is.null(seed)) { set.seed(seed) }
```

- Workers: `workers <- max(1L, min(workers, num_simulations))` (959), where `workers` is
  `length(cluster)`, `num_workers`, or `detectCores() - 1` (951-958). One chunk per worker. Within a
  chunk, trials run in order on one stream seeded once.
- If `num_simulations < getOption("evolveTrial.sequential_threshold", 100)`, or if parallel is off, one
  chunk runs with `seed = NULL` on the caller's stream (968-974, 1054).
- Chunk seeds use `sample.int(2^31 - 1, 1)`, which is `R_unif_index` under "Rejection". That is in R and
  stays in R. A check showed it equals `R_unif_index + 1`.
- The FORK cluster is created fresh in each `run_simulation_pure` call (1024). `newForkNode` does not
  reseed (checked by deparse), and `set.seed` in the chunk resets the state anyway.
- The frozen call chain is: wrapper `set.seed(seed)` (warmstart_wrappers.R:521/885), then
  `run_scenarios(..., parallel = FALSE, seed = seed)` (sim_fun_evolveTrial_unified.R:380-386), which calls
  `set.seed(seed)` again. Scenario 1 then uses seed draws 1..k, and scenario 2 (the next scenario, e.g.
  alt then null) uses draws k+1..2k from the same parent stream. Scenario order is therefore part of the
  seeding.
- **Results depend on the worker count.** The number of chunks sets both how many seeds are drawn and
  which trials share a stream, and a trial's draws depend on how many draws came before it in its chunk.
  `splitIndices` sizes: n=10,000 with 19 workers gives 518x2, 527x12, 528x5. n=10,000 with 16 workers
  gives 617x2, 626x12, 627x2. n=3,000 with 16 workers gives 185x2, 187x2, 188x12.
- **Workers used in the frozen runs (conflicts with the brief).** The P6 recert/rep and verify runs used
  **19**, not 16:
  - `docs/reviews/jasa_round1/e2-p6-recert-harness-report.md:145` and `p6-reported-harness-report.md:93,246`
    give `--cpus-per-task=19`.
  - `pilot/verify/freeze/verify_tasks.tsv` sets cpus = 19 for every stageA/stageB task (2 for build/equiv).
  - The P6rep rows in a previous session's scratchpad show `cpus = 19`, `n_workers_seen = 38`
    (19 FORK children x 2 scenario calls).
  - `p6recert.sbatch` and `p6rep.sbatch` default to 4 and 2 CPUs and are overridden at submission.
  - `workers` = `SLURM_CPUS_PER_TASK` (p6recert.R:55, problems.R:246-271).

  16 CPUs applies to the E2/E2c search tranches (`pilot/e2c/longleaf/submit_e2c.sh:176`,
  `pilot/longleaf/e2/submit_bench.sh:17`). E2c certification jobs used 19 (`submit_e2c.sh:179`). So
  bit identity must be to the 16-chunk result for search evaluations and to the 19-chunk result for
  certification, P6 and verify. It is not worker-count invariant.

## 3. What a chunk returns and how results are merged

`simulate_chunk` returns a plain list (simulation_driver.R:885-898) holding:
- named numeric vectors per arm: `sum_final_n`, `sum_final_events`, `sum_stop_eff`, `sum_stop_fut`,
  `sum_final_eff`, `sum_final_fut`, `sum_final_inc`;
- the scalar `sum_final_time` and `n_sims`;
- optionally, with `return_variance = TRUE` (set by the frozen sim_fun), the per-trial vectors
  `final_n_raw` and `final_events_raw` (a list of numerics per arm, in chunk order) and `final_time_raw`.

`aggregate_results` (901-944) adds the chunk sums with `<<-` in chunk order and copies the raw vectors
into position by a running offset. There is no rbind or data.frame. Floating-point detail: the double
accumulators (`final_time`, `final_n`, `final_events`) use plain double `+`, first per trial in chunk
order (856-863) and then across chunks in chunk order. A port must keep both orders. The `var()`,
`quantile()` and data.frame steps (1060-1135) stay in R.

Per-trial and per-look structures a C++ port needs:
- `state` (`make_state`, state_management.R:70-105): `arm_status` (recruiting / stopped_efficacy /
  stopped_futility), `enrolled_counts`, per-arm registries (id, enroll_time, true_event_time,
  random_censor_time, capacity `ceiling(1.2 * Nmax)`), `registry_row_idx`, `stop_efficacy_per_sim_row`,
  `stop_futility_per_sim_row`, `sim_final_n_current_run`, `stop_time`, `.__backstop_fired__` (milestones
  only).
- Loop scalars: `current_time`, `next_calendar_look`, `last_enroll_time`, `patient_id`, and the
  rebalance state (`interval_cutpoints_current`, `rebalance_done`).
- Slice (state_management.R:115-149): `observed_time`, `event_status` (rows with time_available > 0), and
  metrics `events_per_interval` (findInterval left-closed plus `pmin(idx, K)`), `person_time_per_interval`,
  `events_total`, `person_time_total`, `median_followup`.
- Per-trial final vectors: eff/fut/inc/n/events per arm, and `final_time`.

## 4. Feasibility of a draw-order-preserving C++ port

| Obstacle | Exact match possible? | Evidence |
|---|---|---|
| `rexp(1, rate)` | Yes: `R::rexp(1.0 / rate)`. Rcpp's `R::rexp` takes a **scale**. | checked identical over 2,000 seeds |
| `sample(arms, 1, prob)` (ProbSampleNoReplace, revsort ties) | Yes: `Rcpp::sample(x, 1, false, probs)` (a port of R's FixupProb/revsort/loop), or call `revsort` + `unif_rand()` directly. It is not `R_unif_index`, so `sample.kind` does not matter. A length-1 arm set must still use one uniform. | checked identical with 1 and 2 arms, equal probabilities. Unequal probabilities not checked (not used in frozen designs). |
| `rweibull(1, k, scale)`, `runif(1, 0, c)` | Yes: `R::rweibull`, `R::runif`, the same libR functions (including the non-finite, zero-draw case). | checked identical |
| `rgamma` (C++ loop and R-level vectorized) | Yes: `R::rgamma(shape, 1.0 / rate)`. The R vectorized call with vector shapes/rates draws element by element in the same order. | checked identical for a >= 1, a < 1, and vector shapes |
| `rnorm(S, m, s)` (PH) | Yes with `R::rnorm(m, s)` in a loop. **No** with `arma::randn` (polar method). The existing `sample_vs_ref_medians_ph_cpp` must not be reused. | checked: R::rnorm identical, arma::randn differs |
| Vectorized vs scalar order | Not a problem. Every vector draw in this path is 1-D and consumed in index order. The C++ gamma matrix is filled column by column (interval outer), the same as the R fallback. | code |
| Data-dependent draw counts | Yes, but only if the control flow is copied exactly: one look per arrival at most, the arrival that ends accrual still uses its rexp, the look comes between rexp and sample, the D1 phantom-arrival rules, gates, rebalancing at >= 50 events, vs-ref interim using `num_posterior_draws` rather than the interim count, and the discarded V0 control-cache draws when there are more than 1 experimental arm. | Section 1 |
| Floating-point reductions that drive branches | Yes, but each must be copied by hand. R's `sum()` on doubles accumulates in long double. `mean()` uses long double with a second correction pass (logical/integer: long double sum / n). `stats::median` is `mean()` of the two middle values. `quantile` type 7 uses `(1-h)*x[lo] + h*x[hi]`. `seq(..., length.out)` and `findInterval` also apply. A single flipped comparison (gate, `p >= thr`) changes the later draw count and desynchronizes the rest of the chunk. Also compile with `-ffp-contract=off`. The local build has plain `-O2`, no `-march`, so no FMA. The Longleaf `Makeconf` flags were not checked. | code; R semantics |
| `.Random.seed` / RNGkind manipulation | None in the path. A C++ chunk function needs one RNGScope (automatic with Rcpp attributes). It must avoid R callbacks that draw. | grep |
| Parallelism | R's RNG is not thread-safe and the stream is per chunk. Keep the R driver (`splitIndices`, `seed_list`, `set.seed` per chunk, FORK workers) and port only the body of `simulate_chunk` as one single-threaded C++ call per chunk. No OpenMP inside a chunk. | Section 2 |

**Verdict: feasible.** Every RNG primitive on both frozen paths has a C++ call that is bit-identical in
value and stream position. The port must keep the R-level seeding and chunking unchanged.

**Hardest obstacle.** The RNG calls themselves are not the problem. The hard part is copying R's
floating-point reductions exactly (long-double `sum`/`mean`, `median`, quantile type 7) together with the
exact branch order. These decide which gates pass and when arms stop, which sets how many gamma/normal
draws are consumed. One mismatched comparison shifts the stream for every later trial in the chunk. The
PH shadowing trap (the R loop is live, the arma-based C++ version is dead and does not match) is the main
RNG-specific trap.

## 5. Dropout with `random_dropout = FALSE`

Confirmed. `simulation_driver.R:657` always runs `t_random_censor <- runif(1, min = 0, max = censor_max_time_sim)`,
and line 658 `if (!random_dropout) t_random_censor <- Inf` overwrites it afterwards. The uniform is
therefore used in both modes, and the stream is unchanged as long as `censor_max_time_sim` is finite and
positive (frozen: `fu_time`). The roxygen at 269-271 says the same. Edge case: if `censor_max_time_sim`
were `Inf`, `runif` would use zero draws and return NaN in *both* modes. The stream would still match
between the modes, but not with a port that draws unconditionally.
