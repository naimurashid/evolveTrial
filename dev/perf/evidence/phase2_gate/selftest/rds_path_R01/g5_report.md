# G5 statistical equivalence report

Verdict: **SUBSET-PASS (not a gate result)** (1 of 1 designs pass)

Input: 5 rds files under /tmp/claude-1000/-home-naimrashid-Downloads-adaptive-trial-bo-paper/b23583c2-0578-4df7-b82f-32abbed34f31/scratchpad/fake (mode p2 required; engine_ok, checks, R = 10000, no duplicates, variance present)
Engine identity: worktree 0854a589c47de31da8dcf0487d5a45d7cbf692cb P2_ENGINE_SHA NA
Reference: dev/perf/evidence/acceptance_set_reference_rows.csv
Margins: power 0.010, type I 0.005, EN 1% of the pooled reference mean; each one-sided test at level 0.05 (z > 1.6449).

Each z pair is (z_lo, z_hi) of the two one-sided tests; a test rejects when its z exceeds the critical value.

| design | seeds | rows ok | power ref | power p2 | diff | z (lo, hi) | ok | type I ref | type I p2 | diff | z (lo, hi) | ok | EN0 ref | EN0 p2 | z (lo, hi) | ok | EN1 ref | EN1 p2 | z (lo, hi) | ok | class | design |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R01 | 5 | yes | 0.8214 | 0.8214 | +0.0000 | 4.1, 4.1 | yes | 0.0141 | 0.0141 | +0.0000 | 6.7, 6.7 | yes | 11.2223 | 11.2223 | 3.2, 3.2 | yes | 28.2242 | 28.2242 | 5.4, 5.4 | yes | yes | PASS |
