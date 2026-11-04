# G5 failure diagnostic: R01 expected sample size under the null

Date 2026-11-04. Gate run bc47c0f (engine 3bb0694). Verdict under the prespecified rule: **FAIL**, 9 of 10 designs pass.
The one failing check is R01 `EN_null`, phase 2 11.3280 against the pooled reference 11.2223 (+0.94%, upper one-sided z
0.2 against the 1% margin). Every power, type I error and classification check passed on all 10 designs. The margins
were not changed and the gate was not re-run.

## Question

Is the gap a property of the exact engine, or of where the five reference seeds landed?

## Method

`dev/perf/phase2/g5_diag_draws.R` runs R01 through the captured certify arguments (`dev/perf/evidence/args/args_R01.rds`,
which reproduce the reference row exactly at seed 900001, R 10,000, 19 workers, `capture_report.md`), changing only the
seed, R, the worker count and the posterior draw count. Fresh seeds 772001 to 772010, R 20,000 each, 14 workers, so
200,000 trials per arm of the comparison. Phase 1 is `dev_l2` (engine 1c91871, bit identical to the engine of record), at
its production 5,000 draws and at 500. Phase 2 is `dev_p2` (3bb0694). A smaller first pass (seeds 771001 and 771002,
40,000 trials, draws 500, 5,000, 50,000) is in this directory as well.

## Result (`big/summary.csv`)

| engine | trials | EN_null (SE) | EN_alt (SE) | type I | power |
|---|---|---|---|---|---|
| phase 1, 500 draws | 200,000 | 11.293 (0.013) | 28.164 (0.019) | 0.01469 | 0.81819 |
| phase 1, 5,000 draws | 200,000 | 11.290 (0.013) | 28.206 (0.018) | 0.01526 | 0.82020 |
| phase 2, exact | 200,000 | 11.285 (0.012) | 28.203 (0.018) | 0.01481 | 0.82021 |

- Exact minus phase 1 at 5,000 draws, `EN_null`: -0.005 (SE 0.018, z -0.29, -0.05%). `EN_alt`: -0.003 (z -0.11).
- Phase 1 at fresh seeds minus the pooled reference (11.2223): +0.068, z 2.41. All five reference seeds (11.16 to
  11.28) lie below the phase 1 fresh-seed mean.

## Attribution

On R01 the exact engine and the engine of record agree on expected sample size to 0.05%. The G5 failure comes from the
reference rows' own Monte Carlo position (about 2.4 standard errors low), not from the engine. The independent review's
point H1 (finite posterior draws add futility stops) is real, but on R01 its size is below the noise at 200,000 trials.
On R23 the review measured it on power (0.828 exact, 0.813 at 800 draws, at R 3,000), and G5 passed R23 power at +0.006.

Whether to accept phase 2 on this attribution, or to re-specify G5 (for example, compare against fresh phase 1 runs
rather than stored rows), is the PI's decision. Nothing here changes the recorded verdict.
