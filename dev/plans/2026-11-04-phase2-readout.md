# Phase 2 readout (statistically equivalent engine)

Date 2026-11-04, run unattended overnight on the PI's authorization of 2026-11-03. Plan `2026-11-04-phase2-plan.md`
(gates fixed in 0854a58 before any run). Branch `perf/engine-cpp`, local only, not pushed. Final engine 669ca3e, library
`dev/perf/lib/dev_p2_final` (built by `layer2/build_dev_tree.sh`, lineage fields in DESCRIPTION).

## 1. What was built

- Exact posterior probabilities (`src/posterior_exact.{h,cpp}`). Between-arm decisions use the closed-form normal tail
  of the Laplace posterior on the log hazard ratio, which is what the phase 1 draws estimated. Single-arm decisions use
  the tail of a weighted sum of independent gamma hazards (median below the threshold exactly when the cumulative hazard
  at the threshold exceeds log 2), by `pgamma` for one interval and contour-shifted characteristic-function inversion
  otherwise. Plain Gil-Pelaez on the real line converged too slowly at the designs' small priors (P2-1 report).
- Phase 2 chunk engine (`src/simulate_chunk_exact.cpp`, `R/simulate_chunk_exact.R`). Phase 1 trial logic, per-trial
  xoshiro256++ streams seeded from (scenario seed, trial index), dynamic block dispatch on the existing cluster
  (`clusterApplyLB`), sums in trial-index order. Opt-in through `options(evolveTrial.engine = "exact")`. The default
  stays the phase 1 engine, and its results and R's random stream are unchanged (checked against `dev_l2` several times).
- Faster look slicing in the phase 2 engine only (P2-4), `identical()` in output to the engine before it.

## 2. Gates

| Gate | Result | Evidence |
|---|---|---|
| G1 exact probabilities | PASS. 2,551 deterministic rows, max error 1.0e-13. 200 Monte Carlo cases at 10^7 draws, max abs z 2.86. Between-arm 300 cases, max abs z 2.65. Coordinator's independent R check 12 of 12 | `phase2/P2-1_REPORT.md`, def87a3 |
| G3 determinism | PASS 100 of 100 (identical at 1, 4, 16, 19 workers and two block sizes, 10 designs) | `phase2/g3_results.csv`, 3bb0694, rerun at 923caae |
| G4 no regression | PASS. Package tests 1013 passed, 0 failed, 33 skipped. Phase 1 layer 1 0 mismatches in 2,120 chunks | 3bb0694, 669ca3e |
| G5 statistical equivalence | **FAIL, 9 of 10 designs.** All power, type I and classification checks pass on all 10. R01 `EN_null` 11.328 against reference 11.222 (+0.94%). Two one-sided tests need z above 1.645 on both sides, so the difference plus 1.645 standard errors must stay inside the 1% margin. The upper z was 0.2 | `evidence/phase2_gate/g5/`, bc47c0f |
| G5 diagnostic | At 200,000 fresh-seed trials each, exact minus phase 1 on R01 `EN_null` is -0.005 (z -0.29). Phase 1 at fresh seeds sits +0.068 (z 2.41) above the stored reference. The failure is the reference rows' Monte Carlo position | `evidence/phase2_gate/g5_diag/G5_DIAGNOSTIC.md`, d2cc56d |
| G6 speed | PASS. R23 at 1 worker 0.39 s per 1,000 trials elapsed (0.11 s CPU) against phase 1's 9.7 s on the same machine. R01 0.39 against 8.4 | `evidence/phase2_gate/timing/`, 0ae175f |
| G6 efficiency | **Not met.** 16-worker efficiency 0.25 (phase 1 0.45 to 0.50) at the step 7 protocol, 0.38 for R23 at 20,000 trials. Fixed parallel start-up now exceeds the work | same |

Timing note. In `timing_p2.csv` the column `per_1000_s` is the elapsed time of one `run_scenarios` call, which runs two
scenarios of 1,000 trials each, so per 1,000 trials it is half the column. The machine carried load 7 to 8 from another
session's job during G6, so absolute times are rough and the ratios are the useful numbers. The G5 run (1.5 million
trials, 75 rows) took 5 minutes locally. Phase 1's layer 2 took 7.1 core-hours on Longleaf.

## 3. Review

Independent Opus review of 0854a58..3bb0694 (`phase2/REVIEW_phase2.md`): 0 critical, 1 high, 2 medium, 5 low. The high
item is the G5 attribution question. The exact engine is the infinite-draw limit of the phase 1 rule, and finite draws
add futility stops (R23 power under the alternative at R 3,000: 0.828 exact, 0.813 at 800 draws, 0.784 at 100). M1, M2,
L1 and L2 are closed in 669ca3e (`phase2/REVIEW_FIXES.md`, nine reverted mutations all caught). L3 to L5 are notes. A
false positive from the commit hook (duplicate-index storage in the replay hook) was checked against phase 1 and
discarded.

## 4. For the PI

1. **G5 acceptance.** Accept phase 2 on the diagnostic's attribution, or re-specify G5 (for example, compare against fresh
   phase 1 runs at matched R instead of the stored rows). The recorded verdict stays FAIL until you rule.
2. **G6 efficiency.** The prespecified criterion was not met, for a reason that did not exist when it was written. If
   multi-worker throughput matters, the next step is in-process threads in C++ (each trial owns its stream, so there is no
   R call inside the loop), not more R workers.
3. **Margins** (0.01 power, 0.005 type I error, 1% expected sample size) were chosen without you (plan Section 6).
4. **The single-arm inversion method** replaces the Gil-Pelaez line the plan named (P2-1 report Section 2).
5. With exact probabilities no posterior draw count is left, so "draw count as fidelity" (design Section 6) is moot.
