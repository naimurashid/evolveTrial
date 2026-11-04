# P2-1 report: exact posterior probabilities (phase 2, step P2-1)

Date 2026-11-04. Worktree `~/research/evolveTrial-perf`, branch `perf/engine-cpp`, nothing committed, installed only to
`dev/perf/lib/dev_p2`. Plan `dev/plans/2026-11-04-phase2-plan.md` Section 1 item 1, gate G1 (Section 3).

**Result: G1 PASS.** 2,400 single-arm cases, maximum absolute error 9.99e-14 against every deterministic reference
(2,551 reference rows), 200 Monte Carlo cases at 10^7 draws with max |z| 2.86, between-arm closed form identical to
`pnorm` on the phase 1 Laplace output (500 cases) and 300 Monte Carlo comparisons at 10^6 draws with max |z| 2.65.
Package tests 783 passed, 0 failed, 33 skipped (phase 1: 723 / 0 / 33; the 60 new expectations are the new file).

## 1. Verified parameterizations (phase 1 code, file:line)

| Item | Where | What the code does |
|---|---|---|
| Posterior shape | `src/posterior_sampling.cpp:49` | `prior_alpha(j) + events(j)` (double + int) |
| Posterior rate | `src/posterior_sampling.cpp:50, 53-57` | `prior_beta(j) + PT(j)`, replaced by `1e-6` (with an R warning) when `<= 0` |
| Draw | `src/posterior_sampling.cpp:61-63` | `R::rgamma(shape, 1 / rate)`, i.e. Gamma(shape, rate) |
| HC call | `src/simulate_chunk.cpp:330-338` | `posterior_medians` passes E, PT, prior alpha, beta, `lens` unchanged |
| HC median | `src/posterior_sampling.cpp:98-133` | intervals start at 0, `current_time += delta_t_j`; past the grid it extends with the last hazard; `Inf` if the last hazard is 0 |
| PH median | `src/simulate_chunk.cpp:397-418` (`median_piecewise_r`) | same algorithm, `std::log(2.0)` (same double as the literal `0.69314718055994530942`, checked `0x1.62e42fefa39efp-1`), stops on NaN instead of propagating |
| HC tallies | `src/simulate_chunk.cpp:340-354, 375, 763, 767` | interim futility `mean(med < thr)`, final efficacy `mean(med > thr)`, final futility `mean(med < thr)`; NA if any median is NaN |
| PH log-HR | `src/simulate_chunk.cpp:443-448` | `ph_beta_mode_var_cpp` then `R::rnorm(bmean, bsd)` |
| PH probabilities | `src/simulate_chunk.cpp:514, 741, 742` | `mean(logHR >= m)` (futility), `mean(logHR < m)` (efficacy) |
| Laplace sd | `src/posterior_sampling.cpp:232, 250-255` | `var = min(1e6, max(1e-6, -1/hess))` if `hess < 0`, else 1; Newton breaks before a non-finite step |
| Cutpoints | `src/r_semantics.h:148-190`, `src/simulate_chunk.cpp:652` | rebalanced cutpoints always start at 0, `lens = diff(cut)` |

The two median functions agree. Code reading: identical loops except NaN handling. Numerically (G1 check
`median_functions_agree`): 40,000 medians from `.chunk_ph_sample` (which uses `median_piecewise_r`) are bit-identical to
`calculate_median_survival_matrix_cpp` on replayed draws, including cases with a zero last hazard (open-ended tail
giving `Inf`). The `lambda == 0` branch inside the loop (`posterior_sampling.cpp:115`, `simulate_chunk.cpp:405`) is
unreachable: it needs `cum + 0 * dt >= log 2` with `cum < log 2` on entry (otherwise the previous iteration returned).
The identity `med < thr` iff `sum_k w_k lambda_k > log 2` held on 1.38 million phase 1 draws with zero hazards and
thresholds inside, at and beyond cutpoints, 0 mismatches outside near ties (`median_identity`).

## 2. Method

### Single-arm: Pr(med < thr)

Weights: `w_k = min(thr, start_k + lens_k) - start_k` for `start_k < thr` (last interval: `thr - start_k`), else 0,
with `start_k` by repeated addition as in the median code. With `theta_k = w_k / rate_k`, `S = sum theta_k G_k`,
`G_k ~ Gamma(a_k, 1)`, and Pr(med < thr) = Pr(S > log 2), Pr(med > thr) = 1 - that.

1. Terms with `a_k = 0` or `theta_k = 0` are identically 0 (phase 1: `R::rgamma(0, s) = 0`) and are dropped. Terms with
   exactly equal scale are merged (`Gamma(a1, t) + Gamma(a2, t) = Gamma(a1 + a2, t)`).
2. No term left: 0. One term: `R::pgamma(log 2, a, scale = theta, upper)`.
3. Two or more: Laplace inversion on a deformed Bromwich contour. Scale c = log 2 to 1 (`tau = theta / c`). With
   `L(s) = prod (1 + tau_k s)^(-a_k)`, analytic off the cuts `(-inf, -1/tau_k]`,
   `Pr(T <= 1) = (1/2 pi i) int_G e^s L(s)/s ds` for a contour G right of the pole at 0 and of every cut. Moving the
   vertex between `-1/tau_max` and 0 subtracts the residue 1 at the pole, so the same integral is `-Pr(T > 1)`. The
   vertex `gamma` is the real minimum of `|e^s L(s)/s|` on each side (both convex in s, unique, found by safeguarded
   Newton; the left side is parametrized by `v = 1 + tau_max s` so `1 + tau_k s` keeps full relative accuracy near the
   branch point). The side with the smaller saddle value is used, so the smaller of the two probabilities is computed
   directly. The contour is the hyperbola `s(y) = gamma + i y - (1/2)(sqrt(y^2 + b^2) - b)`, `b = 1/sqrt(g''(gamma))`:
   vertical at the vertex (steepest descent) and asymptotically at slope 1/2 into the left half plane, so the integrand
   decays like `e^(-y/2)` whatever the total shape. By conjugate symmetry the probability is
   `(1/pi) Im int_0^inf e^s L(s)/s s'(y) dy`.

Why not the vertical Gil-Pelaez line the plan names. On the line the integrand decays only like `y^-(1 + sum a)`; with
the acceptance designs' priors (alpha = beta = 0.22 to 0.48) and no events, `sum a` is about 0.4 to 1, and truncating at
1e-9 would need y near 10^10 to 10^25 with oscillation throughout. Bending the same inversion integral into the left half
plane (allowed because L is analytic there apart from the cuts on the negative real axis) removes that problem. The
Moschopoulos series is used only as a reference: it converges slowly when the scales differ by orders of magnitude (a
threshold just past a cutpoint gives a tiny weight), which happened in 136 of 2,400 gate cases within 4,000 terms.

Error control.
- Truncation of y at Y: rigorous. For `y >= Y`, `|e^s| <= e^(gamma + b/2 - y/2)`, `|s| >= y`, `|1 + tau s| >= tau y`,
  `|s'| <= sqrt(5)/2`, and `prod max(1, (tau_k y)^-a_k)/y` is non-increasing, so the neglected integral is at most
  `2 sqrt(5)/2 e^(gamma + b/2 - Y/2) prod max(1, (tau_k Y)^-a_k) / Y`. Y doubles from the initial panel width until this
  is below 1e-16 pi.
- [0, Y]: panels at 0, h0, 2 h0, 4 h0, ..., Y (h0 = min(b, distance from the vertex to the nearest singular point)),
  globally adaptive Gauss-Kronrod 10/21 (QUADPACK qk21 nodes), bisecting the worst panel until the summed `|K21 - G10|`
  over pi is at most 1e-12 (cap 2,000 panels; a `converged` flag reports the cap). `|K21 - G10|` overstates the K21
  error by orders of magnitude on analytic integrands, so 1e-12 is a conservative a posteriori estimate, not a proof.
- Rounding: `DBL_EPSILON * int |f| / pi` is reported. On the 2,094 contour cases of G1: max error estimate 9.7e-13,
  max rounding estimate 2.7e-16, all converged, about 126 to 190 integrand evaluations.

Edge cases, mirroring the phase 1 tallies: `thr <= 0` or `-Inf`: below 0, above 1. `thr = +Inf`: below 1 when the last
interval has positive shape, else the finite cumulative hazard over the grid. `thr` NaN: both 0 (every comparison in
`tally` is FALSE, so phase 1 returns 0 for both, not complementary). Rate `<= 0`: replaced by 1e-6, `status = 1` (no
warning emitted; see Section 6). Any non-finite E, PT, alpha, beta or length, a negative length, or a negative shape:
both NaN, `status = 2`.

### Between-arm

`loghr_probs(bmean, bsd, m)`: `ge = pnorm(m, bmean, bsd, upper)`, `lt = pnorm(m, bmean, bsd, lower)`, with R's
`rnorm` semantics (R 4.6 `src/nmath/rnorm.c`): NaN mean, or sd NaN, infinite or negative, gives NaN draws, so both
NaN (phase 1: `logical_mean` is NA and the engine stops at `simulate_chunk.cpp:515, 743, 747`); sd 0 or infinite mean
gives the constant `bmean`, so `ge = (bmean >= m)`, `lt = (bmean < m)`; NaN m gives NaN.
What `ph_beta_mode_var_cpp` can return: sd is always in [1e-3, 1e3] or exactly 1 (if `hess` is NaN, `hess < 0` is FALSE
and var = 1), so sd 0 and NaN sd do not occur. The mean is NaN only for NaN or negative inputs (the Newton loop breaks
before accepting a non-finite step, so a finite start stays finite).

## 3. Gate G1 (`dev/perf/phase2/g1_check.R`, results `g1_results.csv`, log `g1_check.log`)

Inputs: 2,400 cases, 300 per count of positive weights 1 to 8, thresholds 834 inside an interval, 855 exactly at a
cutpoint, 862 beyond the last cutpoint (counts among cases with a deterministic reference; 800 per mode generated);
events 0 to 300 per interval (30% zero), weak priors alpha = beta in 0.05 to 1 plus 5% of cases with alpha = 0 in one
interval; equal-length grids as in the designs (40%) and random lengths.

| Check | Rows | Max abs error | Max abs z | Failed |
|---|---|---|---|---|
| single-arm vs pgamma (one term left) | 304 | 0 | | 0 |
| single-arm vs R integrate (two terms) | 292 | 2.2e-14 | | 0 |
| single-arm vs Moschopoulos, long double (2 to 8 terms, where converged) | 1,955 | 1.0e-13 | | 0 |
| single-arm vs phase 1 sampler, 10^7 draws | 200 | 4.1e-4 | 2.86 | 0 |
| hook vs gamma-sum tail on weights recomputed in R | 2,400 | 0 | | 0 |
| median identity on phase 1 draws | 1.38e6 draws | 0 mismatches | | 0 |
| `median_piecewise_r` vs `calculate_median_survival_matrix_cpp` | 40,000 | 0 mismatches | | 0 |
| between-arm: identical to `pnorm` on `ph_beta_mode_var_cpp` output | 500 x 3 margins | 0 | | 0 |
| between-arm vs `.chunk_ph_sample`, 10^6 draws (both sides) | 300 | 1.2e-3 | 2.65 | 0 |
| between-arm rnorm edge semantics | 12 | | | 0 |

Per count of positive weights (deterministic max error / MC max |z|): 1: 0 / 2.00; 2: 1.0e-13 / 2.86; 3: 9.8e-14 / 2.36;
4: 9.8e-14 / 2.23; 5: 2.7e-14 / 2.56; 6: 1.8e-14 / 1.80; 7: 1.7e-14 / 2.37; 8: 2.0e-14 / 2.23. The 136 cases without a
deterministic reference (Moschopoulos not converged within 4,000 terms, three or more terms) are all in the Monte Carlo
set; the plan asked for at least 50 Monte Carlo cases, 200 were run (all 136 plus 64 stratified).

History, for the record. The first run reported FAIL on 5 of 290 integrate rows (errors 6e-7 to 1.7e-5). On all five the
long-double Moschopoulos series agreed with the exact value to 6e-15, and the cause was the reference: for a shape below
1 only the first integration piece used the substitution `v = g^a`, and QUADPACK misjudged the steep density on the next
piece. The substitution now covers every piece; tolerances and pass rules were not changed. The same reference flaw
appeared in my first unit-test helper and was fixed the same way.

## 4. Tests

New `tests/testthat/test-posterior-exact.R` (60 expectations, about 1 s): pgamma reduction, two weights against
integrate on design-like weak-prior inputs at 3.8, 4.527 and 5.999 months, a five-weight Monte Carlo check with the
phase 1 sampler, every edge case above, equal-scale merge, extreme scale ratios and shapes near 300, the between-arm
closed form, Monte Carlo and rnorm edge semantics. Suite against `dev_p2` (`dev/perf/harness/run_pkg_tests.R`):
**783 passed, 0 failed, 33 skipped, 0 errors, 0 warnings** (phase 1 count 723 / 0 / 33).

## 5. Timing (`dev/perf/phase2/p2_1_timing.R`, desktop, one process)

Microseconds per call through the R hook (about 2 us of that is `.Call` overhead), prior alpha = beta = 0.3, events 0 to
20 per interval, against the phase 1 Monte Carlo it replaces (5,000 draws plus medians, same inputs):

| Positive weights | Exact (us) | Integrand evaluations | Phase 1 MC (us) |
|---|---|---|---|
| 1 | 2.5 | 0 (pgamma) | 400 |
| 2 | 12.0 | 129 | 500 |
| 3 | 13.5 | 126 | 800 |
| 4 | 14.5 | 126 | 900 |
| 5 | 16.5 | 126 | 1,000 |
| 6 | 18.0 | 126 | 1,200 |
| 7 | 21.0 | 126 | 1,300 |
| 8 | 21.5 | 126 | 1,400 |

Between-arm `loghr_probs`: about 1 us per hook call including overhead (two `pnorm`). The machine's load average was
about 8 at the start of the timing run (the G1 run had just finished), so these are upper bounds.

## 6. Places where the plan, the brief or the code disagree (reported, not improvised around)

1. **Gil-Pelaez.** Plan Section 1 item 1 and Section 6 item 2 name Gil-Pelaez inversion. I used the same inversion
   integral on a bent contour (Section 2) because the vertical line is the slow case for the designs' small shapes. The
   plan left "method and error control chosen by the implementer"; this is the decision for PI review.
2. **Error bound.** The truncation bound is rigorous; the quadrature error is the conventional Kronrod estimate, not a
   proof. The 1e-9 target is met empirically by four orders of magnitude (max 1.0e-13 over 2,551 reference rows).
3. **R version.** The plan says the desktop runs R 4.6.0 and `src/r_semantics.h` cites R 4.6.0 sources; this machine runs
   **R 4.6.2** (2026-10-29), which built `dev_p2`. It does not affect P2-1 (no R-semantics helper is used), but the
   coordinator may want to know before G4 and G6.
4. **Rate guard warning.** `draw_posterior_hazard_samples_cpp` emits an R warning when a posterior rate is `<= 0`. The
   exact function applies the same replacement and returns `status = 1` but does not warn; P2-2 decides whether to warn.
   Unreachable in the supported set (beta > 0, PT >= 0).
5. **Invalid inputs.** Phase 1's response to a NaN input is draw-dependent (NaN medians give an NA tally, but a NaN in
   the last interval gives `Inf` medians, `posterior_sampling.cpp:129-133`). The exact function returns NaN for any
   non-finite input in any interval, a deliberate simplification; unreachable in the supported set.
6. **Draw-dependent PH errors not reproduced.** Phase 1 evaluates the gamma draws and both medians only for their error
   conditions (`simulate_chunk.cpp:458-472`). With valid inputs the only reachable one is a draw with
   `exp(logHR) = Inf` (logHR > 709.78), which gives NaN rates or `0 * Inf` and stops. With the prior sd of 1 in every
   design, the Laplace sd is at most 1, so this has probability `pnorm(709.78, bmean, bsd, upper)`, nil for any
   realistic `bmean`. The exact engine does not reproduce that error.
7. **`.exact_pr_loghr` ignores `interval_lengths`.** It keeps the `.chunk_ph_sample` signature as the brief asked; the
   log-HR posterior does not depend on the lengths.
8. **`ph_beta_mode_var_cpp` returns an `Rcpp::List`.** P2-2 calling it per look allocates an R list each time. A plain
   C++ variant would be faster but would duplicate `posterior_sampling.cpp` code, so I left it to P2-2.
9. **Brief's NaN `thr` rule.** The brief says Pr(med > thr) = 1 - Pr(med < thr). That holds for every non-NaN `thr`;
   for NaN `thr`, phase 1 gives 0 for both, and the exact function mirrors that.
10. **Other agents' files.** `dev/perf/phase2/` already held P2-3 files (`P2-3_REPORT.md`, `p2_*.R/sh`) and
    `dev/perf/evidence/phase2_gate/` exists untracked. I did not touch them.

## 7. Files changed

New: `src/posterior_exact.h`, `src/posterior_exact.cpp` (plain C++ in namespace `evt_exact`:
`gamma_sum_upper_tail`, `median_probs`, `loghr_probs`; hooks `.exact_gamma_sum_tail`, `.exact_pr_median`,
`.exact_pr_median_below`, `.exact_pr_median_above`, `.exact_loghr_probs`, `.exact_pr_loghr`, all `rng = false`),
`tests/testthat/test-posterior-exact.R`, `dev/perf/phase2/g1_check.R`, `dev/perf/phase2/ref_moschopoulos.cpp` (gate
reference, compiled by sourceCpp, not in the package), `dev/perf/phase2/p2_1_timing.R`, `dev/perf/phase2/g1_results.csv`,
`dev/perf/phase2/g1_check.log`, this report.
Regenerated by `Rcpp::compileAttributes()`: `R/RcppExports.R`, `src/RcppExports.cpp` (additions only).
Unchanged: `src/simulate_chunk.cpp`, `R/simulate_chunk_cpp.R`, `R/simulation_driver.R`, `src/posterior_sampling.cpp`,
`src/Makevars`. Library: `dev/perf/lib/dev_p2` (created).

Reproduce: `R CMD INSTALL -l dev/perf/lib/dev_p2 .`, then `Rscript dev/perf/phase2/g1_check.R dev/perf/lib/dev_p2 2400 64 12`
(about 70 s on 12 cores), `Rscript dev/perf/harness/run_pkg_tests.R dev/perf/lib/dev_p2 <out.rds>`,
`Rscript dev/perf/phase2/p2_1_timing.R`.
