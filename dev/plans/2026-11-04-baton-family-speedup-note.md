# BATON-D, BATON-C and the program (phaseless) project: how the engine speed-up lessons transfer

Date 2026-11-04. Read-only research note for the PI, companion to `2026-11-04-hybrid-seamless-speedup-note.md`. No code
was edited and nothing was committed. One short profile was run on a scratch copy of BATOND (Section 2.2). Its script
and output are in the session scratchpad, not in any repository. Every number is marked (measured), (estimated) or
(computed). Measured numbers cite a file and line or a commit.

## 1. Where the project descriptions conflict with what I found

1. **BATON-D speed-up is not only in worktrees.** It is merged to BATOND `main` (head 7a9a56f) and has been in
   production since 2026-10-29. The 193 remaining Variant B Layer 1 searches run on 7a9a56f with
   `BATOND.sim_engine = "cpp"` in strict mode (baton-dosefinding `docs/DECISIONS-2026-10-29.md` item 2,
   `config/variant_b_layer1.yml` lines 21 to 25), and so does the null-1 Layer 2 work (`config/null1_layer2_part2.yml`
   line 32, BATOND 5bcdba8, which contains 7a9a56f). The port plan's sentence that Layer 1 runs on the R engine
   (`docs/plans/2026-10-28-rcpp-simulator-port.md` line 8) was superseded by that decision.
2. **The `~/research/BATOND` checkout is on `search-reliability`, which does not contain the C++ engine.** Its merge base
   with `main` is 6346b43. Two of its commits change ported R code (d8d7748 BOIN end-of-trial selection honours `n_min`,
   in `R/sim_boin.R`, and 40c83d3 BOIN12's terminal step, in `R/sim_boin12.R`). When that branch meets `main`, the C++
   simulators (`src/select_mtd.cpp`, `src/sim_boin*.cpp`) must get the same changes, or the strict C++ engine will
   disagree with the R reference. The existing fidelity gate catches this, but only if it is re-run.
3. **BATON-C is not on BATON branches.** `baton-c/CLAUDE.md` lines 26 to 31 say the code lives on BATON branches
   (`BATON-C-dev`, `contextual-parameters`, `robust-api`). Those branches still exist, but the live code is the standalone
   package `~/research/BATONC-certified` (origin `ayoung31/BATON-C`, version 0.2.0 at 673e6e3). It imports BATON 0.9.2 or
   later and hetGP, and it does not depend on or call evolveTrial (`DESCRIPTION` lines 18 to 29, and no `evolveTrial`
   string in its `R/`, `tests/`, `vignettes/` or `docs/`). `baton-c/CLAUDE.md` line 45 also lists evolveTrial as a
   package. The `baton-c` repository holds only `CLAUDE.md`, so there are no BATON-C studies to profile.
4. **"Phaseless designs" is now a program-calibration paper.** Since the 2026-10-15 rulings the project calibrates a
   fixed architecture (phase I dose finding, randomized phase II screening, independent phase III confirmation) and
   compares A, A′ and B (`baton-phaseless/CLAUDE.md` lines 10 to 27). The joint-allocation direction of 2026-10-19 is
   recorded but not examined, and no simulation campaign is authorized. It is docs only, as the PI expected.
5. **The exact engine as default is an uncommitted working-tree change** in `evolveTrial-perf`
   (`R/simulate_chunk_exact.R` lines 4 to 9 say "the default since 0.0.0.9002, PI ruling 2026-11-04", and
   `git status` shows that file, `R/simulation_driver.R`, `DESCRIPTION` and `NEWS.md` modified). "Becoming the default"
   is accurate.

## 2. BATON-D (BATOND)

### 2.1 What was already done and measured

The 2026-10-28 performance survey (`baton-dosefinding/docs/notes/2026-10-28-batond-performance-survey.md`) measured the R
engine on the Variant A traces. A candidate evaluation simulates 18,400 phase I trials (computed from lines 23 to 27).
These are 1,600 on the nominal truth and S0, 4,800 on the six nulls of the search stream, and 12,000 for the calibration
of `c` on the six nulls of the calibrate stream. The median cost was 33.5 to 36.3 s per candidate for BOIN12 and 16.4 to
17.3 s for U-BOIN under 14 concurrent workers, and the 160 Variant A searches cost 143.5 core-hours (measured, lines 62
to 67). Everything outside candidate evaluation was 4.0 to 8.4% of wall (measured, lines 69 to 71). The surrogate and
acquisition were 4 to 8% of search wall (measured, line 90). So BATON-D, like evolveTrial, is simulator-bound.

Then four commits, all bit-identical, all on `main`.

| Commit | What | Gate | Speed |
|---|---|---|---|
| 12a86a7 | Vectorised common random number stream (cell table cached, no S3 dispatch in `crn_draw()`) | 65,937 stream objects and 225,978 draws identical | about 2x per stream build, 3x per draw (measured, commit message) |
| 0d59077, 184a85e | C++ simulators for BOIN, BOIN12, U-BOIN and a reimplementation of `BOIN::select.mtd`, behind `BATOND.sim_engine`, `rng = false` | G1 to G7 of the port plan, `identical()` on every replicate, 0 mismatches (`results/rcpp_port/impl_report.md` lines 9 to 19) | BOIN12 23.5 to 5.1 s (4.4x), U-BOIN 9.0 to 4.6 s (1.9x) per candidate at `-O2`, per-trial simulator call 1.6 to 0.3 ms (measured, `impl_report.md` lines 341 to 348) |
| 2927a8b | Duration read from trace columns instead of a data frame, per-bundle argument checks, faster `crn_seed()` | 484,800 trials identical against 8be4349, 0 fallbacks | duration 720 to under 90 microseconds per timed trial, BOIN12 47.0 to 6.6 s (7.1x) and U-BOIN 13.3 to 5.7 s (2.3x) per candidate on the loaded machine (measured, commit message) |
| 7a9a56f | U-BOIN randomized stage 2 and both 3+3 rules in C++ | 2,200 cells, every 3+3 path for 1 to 5 doses exhaustive (6,636,320 calls) | not reported |

None of the survey's search-level options appears in any commit I found (`git log --all` on BATOND). These are O1
stage-ordered evaluation, O4 the exact disk cache, O5 scoring all five expansion sizes `m` per simulation and O6 pinning
workers to performance cores. The in-memory evaluation memo (601b5c0) predates the survey and is what its line 41
describes.

### 2.2 Where the time goes now (measured today)

I ran the port's own profile script (`results/rcpp_port/profile_cpp.R`) against a `git archive` copy of 7a9a56f in the
scratchpad, so no worktree was touched. One core, load about 2.5, `pkgload::load_all()` build (which compiles at `-O0`),
the first stored candidate of `search_c10_A_r1.rds` (BOIN12) and `search_c6_A_r1.rds` (U-BOIN), search replicates
900001 to 900800 and calibrate 1 to 2000, as in the port report. Total under one minute.

| Bucket (share of one candidate) | BOIN12, 4.43 s | U-BOIN, 3.88 s |
|---|---|---|
| C++ simulator loop (`boin12_core_cpp`, `uboin_core_cpp`) | 25% | 24% |
| Stream construction in R (`crn_stream()`, `set.seed`, draws, save and restore of `.Random.seed`) | 20% | 26% |
| R side of the dispatch (`*_cpp_dispatch` minus the core, including `cpp_crn_ok()` 8%) | 18% | 20% |
| `oc_record()` | 9% | 9% |
| `joint_record()` | 8% | 9% |
| Per-replicate closure, argument assembly, `rbind`, duration and the rest | about 20% | about 12% |
| Of the total, calibration of `c` (`calibrate_threshold()`) | 60% | 62% |
| Of the total, the joint step (`joint_step()`) | 86% | 85% |

All rows are measured (profile output in the scratchpad, `prof/profile_7a9a56f.rds`). That is about 0.24 ms per trial
(computed, 4.43 s over 18,400 trials). About three quarters of a candidate is still R bookkeeping around a C++ core.
Production loads the package through `pkgload::load_all()` (`scripts/null1_l2p2/n1b_common.R` line 1109,
`scripts/null1_l2p1/n1a_common.R` line 891), and `src/Makevars` sets only `-ffp-contract=off`. So production most likely
runs the `-O0` build. The port measured `-O2` at 15% faster per candidate for BOIN12 (6.0 to 5.1 s, measured,
`impl_report.md` line 343).

### 2.3 What transfers

- **Phase 1, the rest of the loop in C++ (bit-identical).** The phase 1 lesson was to port the whole trial loop and its
  bookkeeping, not only the hot kernel. In BATOND the replicate loop, stream build and records are still R. One C++ call
  per candidate, truth and replicate block would build each replicate's stream with R's own `set.seed`, `unif_rand` and
  `exp_rand` in the same order (the survey's O2 already proposed this, line 150), run the simulator, and build the
  `oc_record` and `joint_record` rows directly. This is also the in-process batching lesson, since it removes 18,400 R
  dispatches per candidate. A scratch measurement in the port session (C++ cell lookup in `crn_stream()` and a field-level
  trace read) already reached 3.8 and 3.3 s at `-O2` (measured, `impl_report.md` line 373), before the stream was moved.
  Estimated result 0.5 to 1.0 s per candidate, 4 to 8x over today (estimated from the table above, assuming the C++ core
  at `-O2` and the record building in C++ together cost about a third of today's core time plus the remaining calls).
- **An optimised build.** Install with `R CMD INSTALL` or `pkgbuild::compile_dll(debug = FALSE)` and keep
  `-ffp-contract=off`. About 15% for free (measured on the earlier build, above). Its fidelity was gated at `-O2`
  (`impl_report.md` line 389).
- **Look slicing has an analog that is done.** Trial duration was the BATOND equivalent and 2927a8b took it from 720 to
  under 90 microseconds.

### 2.4 What does not transfer, and why

- **Closed-form posteriors.** BATOND has no Monte Carlo posterior. BOIN, BOIN12 and U-BOIN decide on beta-binomial
  quantities computed by `pbeta` (`R::pbeta` in `src/`), which are already exact. Phase 2's main gain has no target here.
- **Per-trial xoshiro streams.** BATOND already seeds each replicate from (scenario, replicate, stream)
  (`crn_seed()`, `R/crn.R` line 63, `set.seed` at line 253), so results do not depend on worker count and common random
  numbers hold across candidates by design. Replacing R's generator would change every replayed record, the certify-block
  accounting in `docs/CERTIFY_BLOCKS.md` and the identity with frozen campaigns. The cost is real and the benefit is
  only the 20 to 26% stream bucket (BOIN12 to U-BOIN), which the bit-identical C++ stream also removes.
- **Worker start-up.** Production parallelism is one search per worker. `BATOND.cores` forks `mclapply` inside a candidate
  (`R/enumerate.R` line 741), but at 4.4 s per candidate start-up is not the bottleneck. It would become one after the
  loop port if `BATOND.cores` above 1 is used, so keep it at 1 and parallelise across searches.

### 2.5 Statistical options with the largest payoff (PI rulings needed)

These are not phase 1 or 2 techniques, but at the search level they pay more than any simulator change left.

- **O1, stage-ordered evaluation.** 62 to 97% of candidates fail a `c`-free constraint (measured, survey lines 132 to
  135), yet every candidate pays for the calibration of `c` and the null simulations, which are 85% of a candidate today
  (measured above). Skipping them for such candidates saves 55 to 84% of evaluation time (estimated in the survey, line
  141, and still valid since the joint share did not fall). It changes the surrogate's training data, so the search
  trajectory changes. A low-fidelity variant keeps the surrogate fed.
- **O5, all five `m` per simulation.** `joint_record()` already holds the expansion responses for every `m`, and scoring
  another `m` costs under 0.1% (survey line 180). About five times more scored programs per simulated trial (estimated in
  the survey). It changes the search.
- **O4, the exact disk cache.** 31.6% of all simulation repeated earlier work (measured, survey line 170). Bitwise-safe by
  its key, so it needs only concurrency-safe writes. Gate is `identical()` on rescored records.
- **Exact conditional expansion (the true phase 2 analog for the joint step).** Given a replicate's phase I end state at
  the selected dose, the expansion responses are binomial with the truth's rate, so the confirmation indicator can be
  replaced by its exact conditional probability (a finite binomial sum). This is phase 2's "infinite-draw limit" for the
  expansion part only. It is not a speed item, since expansion is under 0.1% of time, but it removes the expansion part
  of the Monte Carlo noise in `false_confirm`, which the survey found makes point-feasibility verdicts noisy (line 216).
  It changes the estimator, so it needs a ruling and a fresh-run equivalence gate.
- **Exact enumeration of the phase I outcome tree.** The 3+3 rules are already enumerated exhaustively in the tests
  (7a9a56f). For BOIN, a dynamic program over (current dose, per-dose counts) would give exact operating characteristics
  with no replicates (feasibility for BOIN12 and U-BOIN, with bivariate outcomes and stage 2 randomization, not
  estimated). It would remove Monte Carlo noise from every constraint but conflicts with the stream-based search,
  selection and certify design. Research effort, defer.

### 2.6 Gates

Bit-identical steps use the port plan's own ladder (G1 to G7, `identical()` and md5 on every replicate, no tolerance,
speed reported and never gated). The loop port must also show the final `.Random.seed` is unchanged after a call, as
184a85e and 7a9a56f did. The `search-reliability` merge (Section 1 item 2) re-runs G1, G3, G4 and G7 for BOIN and
BOIN12. O1, O5 and the exact expansion need a search-level gate fixed before any run. Point-feasible candidates scored
identically, then a paired search comparison (same seeds, both rules) on stored cells, judged on the selected design's
certified constraints, not on stored rows.

### 2.7 Order by payoff per effort

1. Optimised build for production (`-O2`). Trivial, about 15% (measured), no gate beyond the existing `-O2` fidelity runs.
2. O4 exact disk cache. Low effort, about 1.45x on campaign core-hours (computed from 31.6%), bitwise-safe.
3. Remaining loop in C++ (stream, records, replicate loop). Moderate effort, 4 to 8x per candidate (estimated),
   bit-identical.
4. O1 stage-ordered evaluation. Low to moderate effort, 2 to 6x per search (estimated from 55 to 84%), PI ruling.
5. O5 all `m`. Low effort, up to 5x scored programs per simulation (estimated), PI ruling.
6. Exact conditional expansion. Low effort, noise not speed, PI ruling.
7. Exact enumeration. Research effort, deferred.

Items 1 to 4 together give roughly 10 to 30x on campaign core-hours (estimated). Before item 3 is worth starting, the
`search-reliability` divergence should be resolved, since a second C++ edit on top of an unmerged R change doubles the
gate work.

## 3. BATON-C (BATONC-certified)

### 3.1 What it simulates

BATON-C takes a user `sim_fun(theta, n_rep, seed, xi_eval)` (`R/bo_calibrate_contextual.R`). The package helper
`make_simulator()` (`R/make_simulator.R` line 89) builds one from a `dgp_fun` and a `model_fun` that each handle one
trial, and runs them with an R `lapply` over replicates and scenarios (lines 123 to 185). The only simulators in the
package are toys. The platform-trial vignette is a vectorised two-stage binomial (`vignettes/platform-trial-workflow.Rmd`
lines 56 to 58). I found no caller of `bo_calibrate_contextual()` outside the package in `trials-research`, the BATON
worktrees or evolveTrial. **BATON-C does not call evolveTrial**, so no path benefits directly from the new exact engine
today.

If a future BATON-C study wraps evolveTrial, the paths that benefit are `run_simulation_pure()` and `run_scenarios()`
(which calls it) on frames inside the phase 1 supported set (`R/simulate_chunk_cpp.R` lines 188 to 216, the
single-arm and between-arm vs-reference checks). Under the default, a frame outside the set falls back silently to the
phase 1 engine (`R/simulate_chunk_exact.R` lines 7 to 9). The hybrid and seamless engine (`compute_hybrid_oc_rcpp()`) is
not affected. Such a `sim_fun` should call `run_simulation_pure()` once per (theta, xi) with all `n_rep` trials. Building
it through `make_simulator()` with one trial per call would pay the per-call cost on every trial, which the phaseless
Gate 1 measured at 0.038 s per call (Section 4.1) against about 0.4 ms per trial for the exact engine (phase 2 readout,
G6 row).

### 3.2 Where the time goes (measured in the package's own report)

The implementation report (`docs/certified-q-implementation-report.md` lines 174 to 189) timed the certified 0.2.0 path
on the toy, reference BLAS, one core per worker, while tests ran on other cores.

| Step | Time |
|---|---|
| P_R for one candidate (M = 1,024 nodes, N = 4,000 trajectories), dominated by `L %*% eps` and the 1,024 by 1,024 Cholesky (`R/acquisition.R` lines 229 to 250) | 2.41 s (measured) |
| E[R] for one candidate | 0.007 s (measured) |
| One certified selection, 26 of 30 candidates scanned before the first eligible one | 60.4 s at 1 core, 31.1 s at 4 (measured) |
| Adequacy check, M = 2,048 | 28.9 s per scanned candidate (measured) |
| Full toy run, B = 3, 15 warm-start rows | about 2.9 minutes (measured) |

So BATON-C's time is in the surrogate and certification, not the simulator. This is the reverse of BATON, BATON-D and
evolveTrial. This machine's R uses the reference BLAS (`libblas.so.3.12.0`, measured with `sessionInfo()`). With a real
evolveTrial simulator at the exact engine's speed, a 1,000-trial evaluation per context would cost about 0.4 s
(estimated), which is still below one P_R.

### 3.3 What transfers

- **A closed-form screen before Monte Carlo (the phase 2 idea).** The certified scan computes the expensive P_R for
  candidates in objective order until one is eligible. E[R] is the cheap closed-form part, and Markov's inequality gives
  P(R at least tau) at most E[R] / tau. Any candidate with E[R] / tau below gamma cannot be eligible, so it can be skipped
  without sampling, and the selected design is unchanged. It is exact, so it belongs in a bit-identical layer (the
  report's T1 test, `identical()` against brute force). Payoff depends on how many scanned candidates are clearly
  infeasible, up to 5 to 10x on a selection like the toy's (estimated, 26 of 30 scanned were infeasible).
- **Batching and threads.** One BLAS call for all scanned candidates' trajectories, and `n_cores` set with BLAS threads
  at 1 as the report advises. The 4-core selection reached an efficiency of about 0.49 (computed from 60.4 and 31.1 s),
  partly from block overshoot, which is the phase 2 G6 lesson again.
- **RNG discipline is already in place.** The certification set and QMC nodes are fixed per run and leave
  `.Random.seed` untouched, and T1 checks 1 and 2 cores give identical results (report lines 160 to 168).

### 3.4 What does not transfer, and why

- **The phase 1 port.** There is no R trial loop of BATON-C's own to port. The simulator is the user's.
- **A closed form for P_R.** P_R is the probability that a Gaussian process trajectory keeps a Q-mass of contexts
  feasible. It has no closed form. Faster ways to compute it are approximations (a low-rank or pivoted Cholesky of the
  posterior covariance at the nodes, or sequential stopping of trajectories once the interval clears gamma). Both change
  the certification rule's numbers and need a ruling and a tolerance gate against the full computation.
- **An optimised BLAS** (OpenBLAS or MKL, which the report expects to cut per-candidate time several-fold, line 187)
  is the cheapest gain, but it changes last-bit results, so T1 and T2's `identical()` comparisons would have to be
  regenerated under the new BLAS. It does not change the method.

### 3.5 Gates and order

1. Optimised BLAS on Longleaf. Trivial, several-fold (estimated in the report). Gate is P_R within its Monte Carlo
   standard error on a fixed corpus and the same selected design.
2. Markov-bound screen in the scan. Low effort, up to 5 to 10x per selection (estimated), bit-identical, gate T1-style
   `identical()` on selections.
3. Batched trajectories across candidates. Moderate effort, deterministic, gate as item 1.
4. Low-rank P_R or sequential trajectory stopping. Research effort, ruling and tolerance gate.
5. If BATON-C ever wraps evolveTrial, call it once per (theta, xi) with explicit `options(evolveTrial.engine = "exact")`
   so an unsupported frame stops instead of falling back.

## 4. The program (phaseless) project

### 4.1 Which simulator it would use

The Gate 1 reuse inspection (`baton-phaseless/docs/GATE1_REUSE_INSPECTION-2026-10-15.md` lines 22 to 41) found that no
stage needs a new simulator. Phase I is BATOND (BOIN, BOIN12, U-BOIN, 3+3, with its streams and timing). Phase II
randomized PFS screening and phase III run through evolveTrial's `run_simulation_pure()` on the vs-reference path, one
program replicate per call, measured at 0.038 s per call for a 120-patient three-arm trial on one core (line 26). Two
default-off evolveTrial changes are needed, E1 a per-replicate return and E2 an event-driven final analysis. The cost
estimate there is 0.1 to 0.2 s per program replicate (estimated, lines 180 to 182). Separately, the 2026-10-22 plan
(`PLAN-2026-10-22-phase2-confirmation-constructions.md` lines 18 to 25) already composes phase I selection probabilities
with the exact Simon two-stage probability per dose (construction B, done), and proposes a seamless reuse of the phase I
beta posterior (construction A, to specify).

### 4.2 How to build it fast from the start

1. **Do not call `run_simulation_pure()` once per replicate.** At 0.038 s per call against about 0.4 ms per trial in
   the exact engine, the call costs about 100 times the trial (estimated). Build E1 into the phase 2 chunk engine as a
   per-trial record, and let the orchestrator group program replicates by their retained-dose set so each group is one
   call with many trials. This is the in-process batching lesson.
2. **Key evolveTrial trials on the program replicate.** Phase 2 seeds each trial from (scenario seed, trial index). If the
   engine accepted explicit trial keys (program replicate, stage, arm) instead of 1 to n, results would not depend on
   grouping or worker count, and B against A′ would be paired beyond phase I. Gate 1 records that pairing is exact only
   for phase I today (line 164 to 166). Whether the exact engine's streams are per arm was not checked.
3. **Build E2 inside the exact C++ engine, not in the R body.** A frame outside the supported set falls back silently to
   the phase 1 engine under the default. Run the project with `options(evolveTrial.engine = "exact")` set explicitly so
   such a frame stops.
4. **Phase I on the strict C++ engine and an optimised build**, reusing stored certified locks and their replays as
   construction B did. The per-trial data construction A needs (patients, responses and toxicities at the selected
   dose) are not in the stored records (plan line 45), so the replay must write them, once.
5. **Exact composition wherever a stage is fixed-sample and binary.** Construction B already does this. For
   construction A, the expansion cohort's confirmation probability given the phase I end state is a finite binomial sum,
   so the expansion need not be simulated (as in Section 2.5). Phase III sized by Schoenfeld and tested by Wald could be
   computed from the normal approximation of the log hazard ratio given true HR and events, which is what the Kieser
   and drugdevelopR benchmark assumes. That is a modelling choice, not a speed-up of the simulated rule, so it needs a
   ruling.
6. **Threads last.** Parallelise across design points and scenarios first. C++ threads inside a program loop only after
   E1 and E2 are in C++ with no R calls inside the loop.

Estimated cost per program replicate with items 1 to 4, about 2 to 5 ms against the 0.1 to 0.2 s Gate 1 estimate, 30 to
100x (estimated from 0.24 ms per BATOND trial measured above and about 0.4 ms per evolveTrial trial, plus R
orchestration). 10,000 replicates per design point per scenario would fall from about 17 to 33 core-minutes to under
one (computed from those estimates).

### 4.3 Gates

1. E1 and E2 default-off, with `identical()` on every existing output (the bit-identical layer) and per-trial records
   that sum to the aggregates exactly.
2. E2 against an independent R implementation at fresh seeds, two one-sided tests on power, type I error and expected
   events, per the PI's phase 2 ruling (fresh runs, never stored rows).
3. External validation of the vs-reference path, which Gate 1 found missing (lines 43 to 46), for example against the
   drugdevelopR case 01.02 closed forms.
4. Determinism across worker counts and groupings once trials are keyed on the program replicate.
5. Any exact composition against brute-force simulation at 10^6 replicates.

## 5. Summary across the three projects

| Project | Where the time goes | Basis | Best next step | Estimated gain |
|---|---|---|---|---|
| BATON-D | R bookkeeping around the C++ core (about 75% of a 4.4 s candidate), calibration of `c` 60% | measured today | `-O2` build, disk cache, rest of loop in C++, then O1 by ruling | 10 to 30x on campaign core-hours |
| BATON-C | GP trajectory sampling for P_R in certification, simulator negligible on toys | measured in the package report | optimised BLAS, Markov-bound screen | several-fold to 10x per selection |
| Program (phaseless) | per-call overhead of one-replicate evolveTrial calls (not yet run at scale) | Gate 1 measurement and estimate | batch E1 into the exact engine, key trials on the program replicate | 30 to 100x per program replicate |
