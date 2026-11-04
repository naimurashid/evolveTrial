# EDGEFIX: the two Codex P2 defects in the exact engine

Date 2026-11-04. Worktree `~/research/evolveTrial-merge`, branch `perf/merge-main`, HEAD 8f03f69 (0.0.0.9003). Nothing
committed, pushed or installed into `~/R/library`. Source: `~/research/evolveTrial-perf/dev/perf/phase2/CODEX_REVIEW_branch.md`.

## Counts

| Item | Result |
|---|---|
| Full suite before, HEAD unmodified (`p30_edge_base`, `edge_base_tests.log`) | 1234 passed, 0 failed, 35 skipped, 0 errors, 0 warnings |
| Full suite after (`p30_edge`, `edge_tests.log`) | 1305 passed, 0 failed, 35 skipped, 0 errors, 0 warnings, 184 s |
| Reconciliation | 1234 + 55 (test-posterior-exact.R 60 to 115, 2 new tests, loops counted per iteration) + 16 (test-engine-equivalence.R 182 to 198, 1 new test) = 1305. Every other file unchanged per file. Skips unchanged |
| New tests against the old build (`p30_edge_base`) | 28 expectations fail (16 in posterior-exact, 12 in engine-equivalence; the remaining new expectations pass there because they test paths, such as the contour below 1e300, that did not change), so they detect both defects |
| Golden outputs (`edge_golden_identical.R`, `identical()`, no tolerance) | 6 of 6 tables and 6 of 6 RNG states identical. `engine_golden.rds` untouched (not in `git status`) |
| Layer 1, 10 designs (`edge_layer1.log`) | 2,120 chunks compared, 0 mismatches (0 result, 0 stream), 542 s, `dev/perf/harness/out/layer1_20261104_091147.csv` |
| `tools::codoc` | 0 entries (before: `random_dropout` missing from the Rd). `tools::checkRd` clean |
| Codex reproductions (`edge_codex_repro.R`, logs `edge_codex_repro_{edge_base,edge}.log`) | before: exact inconclusive in all three; after: R body's decision or fallback, see Section 1 and 2 |

## 1. Invalid priors (predicate)

Prior parameters the C++ engines read (packed by `.pack_chunk_cfg_supported()` and `.pack_vsref()`, read in
`src/simulate_chunk.cpp:177-178, 194-195` and `src/simulate_chunk_exact.cpp:221-222, 238-239`):

| Parameter | Read by | Status |
|---|---|---|
| `prior_alpha_params_model` (gamma shapes, K) | HC posterior (both engines), PH baseline hazards (`ph_beta_mode_var_cpp`, PH gamma draws) | now required finite and `> 0` |
| `prior_beta_params_model` (gamma rates, K) | same | now required finite and `> 0` |
| `ph_loghr_prior_sd` | PH sampler as `ph_sigma` | not added: packed as `max(1e-6, sd)`, the same clamp as the R body (`R/posterior_helpers.R:114`), so the engines never read a value `<= 0`; V10 keeps it finite |
| `ph_loghr_prior_mean` | PH sampler as `ph_mu` | a location, no sign condition; V10 keeps it finite |

The check sits inside item 13g ("well-formed design inputs"), which already required both prior vectors finite. A new
label would have broken the two tests that pin the exact item-label lists (`test-simulate-chunk-cpp-vsref.R:150`). The
failing-item message names 13g. `R/simulate_chunk_cpp.R:38-50`.

Codex reproduction (R01, zero caps, shapes -1, 10 trials): before, R body and cpp final efficacy 1, exact inconclusive 1.
After, the default runs the R body (efficacy 1, `identical()` table and RNG), and explicit `engine = "exact"` stops with
"outside the supported set ... Failing items: 13g well-formed design inputs".

Phase 1 users. The predicate is shared, so phase 1 now refuses these frames and runs the R body. Checked
(`edge_phase1_fallback.R`, 40 trials R01, 8 trials R23, seed 21) on R01 and R23, each with one shape -1, one shape 0, one
rate 0 and one rate -0.5. 7 of 8 frames: old cpp `identical()` to the old R body, new cpp and new default `identical()`
to old cpp. R01 with shape -1 and normal enrolment errors under every engine, old and new. Only the message differs: old
cpp said "missing value where TRUE/FALSE needed (interim futility probability is NA)", and now (R body) it says "missing
value where TRUE/FALSE needed". The zero-caps frame was also `identical()` under old cpp and old R body. So phase 1
results are unchanged wherever a run completes.

## 2. Overflow in weighted gamma scales

Root cause as Codex said: `median_probs()` forms `w / rate`, which is `Inf` for rate 1e-310, and
`gamma_sum_upper_tail()` rejects a non-finite scale as invalid and returns NaN with the default `converged = TRUE` and
`status = 0`. The final HC analysis reads a NaN as "not efficacy, not futility", which gives inconclusive.

A second instance of the same class turned up during the work. The contour itself goes wrong silently, with
`converged = TRUE`, when `tau = scale / c` reaches about 1e307 while still finite. `edge_contour_range.R` compares 429
two-term cases (shapes 0.01 to 5, tau1 1e2 to 1.6e308, tau2 1e-3 to 2) against an independent reference: `integrate()`
in both orders of the two terms, kept only where the two orders agree to 1e-12. On the old build every case up to
tau = 1e306 is within 1e-11, and 35 cases with tau >= 1e307 and shape 0.01 are off by up to 8e-4 (`edge_contour_range_base.log`).

Fix (`src/posterior_exact.cpp`, `posterior_exact.h`):
- `kHugeTau = 1e300`. Above it the contour is not used, which leaves six decades of margin below the observed breakdown.
- `median_probs()` notes when any `w / rate` overflows and passes that case to a new `gamma_sum_upper_tail_log()`,
  whose scales are `log w - log rate`. Nothing forms tau in this path. Non-overflow inputs take the old call unchanged,
  with the same arithmetic.
- `gamma_sum_upper_tail()` sends n >= 2 cases with largest tau above 1e300 to the same rule, using `log T - log c`.
  The n = 1 pgamma path is unchanged, since R's pgamma is accurate for any finite scale.
- `huge_scale_tail()` (method 3). With x = 1/tau <= 1e-300, `log P(a, x) = -a log tau - lgamma(a + 1)` holds to a
  relative error below 1e-300. One term gives the closed form `-expm1(log P)`. With several terms, the probability
  Pr(T <= 1) is at most B = prod over huge terms of P(a_k, x_k). If B <= 1e-12 (the quadrature's error target) the
  tail is 1, and `err_est = B` is recorded. Otherwise it returns NaN with `converged = FALSE`. That happens only when
  a term with posterior rate below about 1e-300 also has posterior shape below about 0.04.
- New `MedianProbs` status 3 means valid input but a non-finite result. `median_probs` sets it whenever the tail is
  not finite.
- `ExactPost::below` (`src/simulate_chunk_exact.cpp`, the L1 site) now stops on status 2 or 3 with "the exact posterior
  probability is not finite (status, method, threshold, bound) ... use options(evolveTrial.engine = \"cpp\")". Before,
  status 2 NaNs were passed on. Status 2 cannot occur in the supported set after Section 1.

Results. Codex case (shape 1, rate 1e-310, no events, threshold 1): tail 1, method 3, status 0, and `identical()` to
`pgamma(log(2) * 1e-310, 1, lower.tail = FALSE)`. Engine level (R01, zero caps, shape 1, rate 1e-310 in interval 1, and
also all intervals): exact now gives final futility 1, the same as the R body. Before it gave inconclusive 1. With
shape 0.01 and rate 1e-310 the exact engine now stops with the status 3 error (the R body decides, since `rgamma` with an
infinite scale gives `Inf`). On the new build the contour-range grid gives NaN (stop) in exactly the 63 tau > 1e300,
shape 0.01 cases, and every other case is within 1e-11 (`edge_contour_range_edge.log`).

Other exact entry points checked for NaN with success:

| Entry point | Finding |
|---|---|
| single-interval pgamma path (n = 1) | same defect. An overflowed scale was rejected as invalid (NaN, status 0). Fixed by the log path, closed form |
| contour, tau in (1e307, overflow) | finite but wrong, `converged = TRUE`. Fixed by `kHugeTau` |
| contour, NaN integrand | already safe. A NaN error estimate never meets the tolerance, so the panel budget runs out, `converged = FALSE`, and L1 stops |
| HC interim | already safe. `ISNAN(pr_fut)` stops (`simulate_chunk_exact.cpp:811`) |
| HC final | the silent site (`!ISNAN(p) && ...` reads a NaN as inconclusive). Now unreachable, because `below()` stops first |
| between-arm `loghr_probs` (pnorm) | NaN only for a NaN mean or an invalid sd from `ph_beta_mode_var_cpp`. Every caller already stops on NaN (interim `:839`, final `:1014, 1018`). With finite inputs pnorm is finite. No change |
| NaN threshold (`median_probs` returns 0, 0, status 0) | by design (phase 1 tally semantics). Unreachable, because item 10 requires finite thresholds |
| test hook `.exact_gamma_sum_tail` with an invalid scale | returns NaN with `converged = TRUE` (TailInfo has no status). Hook only, since the engine never passes invalid input. Left as is |

## 3. Tests

`tests/testthat/test-posterior-exact.R`, two new tests:
- "Codex overflow case". This covers the reproduction (tail 1, method 3), the one-term tiny-shape closed form against
  R's pgamma (1e-15), several terms with one overflowed (tail 1, B <= 1e-12, and the independent log-space integrate
  reference `ref_two_log` to 1e-12), and the tiny-shape case (status 3, NaN, `converged` FALSE).
- "extreme but valid scales". The contour at tau 1e2, 1e150 and 1e299, shapes 0.01 and 1, against `ref_two_log` (1e-10),
  and the huge-scale rule at 1e301 and 1e305 against the same reference (1e-12). The G1 rows do not record their
  scales, so they cannot be reused here, and an independent integrate reference is used instead.

`tests/testthat/test-engine-equivalence.R`, one new test, "Codex P2". (a) R01 with zero caps and shapes -1, shapes 0, or
rate 0: the default is `identical()` to the R body and to cpp, and explicit exact errors naming 13g. The R body gives
efficacy for shapes -1. (b) Rate 1e-310: the exact decisions equal the R body's (futility 1). (c) Shape 0.01 with rate
1e-310: exact errors "not finite" (the forced non-finite).

## 4. Leftovers

`random_dropout` added to the roxygen block (`R/simulation_driver.R:176-180`) and by hand to
`man/run_simulation_pure.Rd`, in `\usage` and `\arguments` (no roxygen run). The two "FALSE keeps 422eb40" comments
were left as they were. NEWS has one bullet under 0.0.0.9003.

## Conflicts with the brief

1. "Returns the correct tail (effectively 1)" holds for the Codex case (shape 1). For an overflowed term with a tiny
   posterior shape the true tail is not close to 1 (for shape 0.01 it is about 0.9992). There the one-term case returns
   the closed form, and the several-term case stops with an error, since nothing cheap resolves it.
2. "Between-arm log-HR prior variance if applicable": judged not applicable (clamped identically in the R body and the
   packer). No check was added.
3. Making the L1 site stop on status 2 as well as status 3 goes one step past "non-finite results set a non-zero
   status". Status 2 was already non-zero but was passed on silently. It is unreachable in the supported set now.
4. The tau = 1e300 cutoff turns the band 1e300 < tau <= 1e306 with shape below about 0.04, where the contour was still
   accurate, into an error. This needs a posterior rate below about 1e-300.

## Files

Changed: `R/simulate_chunk_cpp.R`, `R/simulation_driver.R` (roxygen only), `man/run_simulation_pure.Rd`, `NEWS.md`,
`src/posterior_exact.cpp`, `src/posterior_exact.h`, `src/simulate_chunk_exact.cpp`, `tests/testthat/test-posterior-exact.R`,
`tests/testthat/test-engine-equivalence.R`. New and untracked, in `dev/perf/phase3/`: this report,
`edge_codex_repro.R`, `edge_phase1_fallback.R`, `edge_contour_range.R`, `edge_golden_identical.R`, and logs/RDS
`edge_*`, `p30_build_edge*.log`. Libraries `dev/perf/lib/p30_edge` and `p30_edge_base` and their build trees can be
deleted.
