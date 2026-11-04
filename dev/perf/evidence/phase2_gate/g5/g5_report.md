# G5 statistical equivalence report

Verdict: **FAIL** (9 of 10 designs pass)

Input: 75 rds files under dev/perf/evidence/phase2_gate/rows (mode p2 required; engine_ok, checks, R = 10000, no duplicates, variance present)
Engine identity: worktree 3bb06940d12a976c34692fd7a34e3b77b0063f70 P2_ENGINE_SHA 3bb06940d12a976c34692fd7a34e3b77b0063f70
Reference: dev/perf/evidence/acceptance_set_reference_rows.csv
Margins: power 0.010, type I 0.005, EN 1% of the pooled reference mean; each one-sided test at level 0.05 (z > 1.6449).

Each z pair is (z_lo, z_hi) of the two one-sided tests; a test rejects when its z exceeds the critical value.

| design | seeds | rows ok | power ref | power p2 | diff | z (lo, hi) | ok | type I ref | type I p2 | diff | z (lo, hi) | ok | EN0 ref | EN0 p2 | z (lo, hi) | ok | EN1 ref | EN1 p2 | z (lo, hi) | ok | class | design |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| R01 | 5 | yes | 0.8214 | 0.8204 | -0.0009 | 3.7, 4.5 | yes | 0.0141 | 0.0154 | +0.0013 | 8.3, 4.8 | yes | 11.2223 | 11.3280 | 6.2, 0.2 | NO | 28.2242 | 28.2107 | 5.2, 5.7 | yes | yes | FAIL |
| R02 | 5 | yes | 0.8609 | 0.8592 | -0.0017 | 3.8, 5.3 | yes | 0.0154 | 0.0153 | -0.0001 | 6.3, 6.5 | yes | 13.2804 | 13.2985 | 4.8, 3.6 | yes | 27.6693 | 27.6586 | 7.1, 7.6 | yes | yes | PASS |
| R03 | 5 | yes | 0.9605 | 0.9608 | +0.0003 | 8.4, 7.9 | yes | 0.0123 | 0.0125 | +0.0001 | 7.3, 6.9 | yes | 21.4467 | 21.4991 | 10.3, 6.3 | yes | 29.8521 | 29.8491 | 41.3, 42.1 | yes | yes | PASS |
| R04 | 5 | yes | 0.8251 | 0.8259 | +0.0008 | 4.5, 3.8 | yes | 0.0057 | 0.0060 | +0.0003 | 10.9, 9.8 | yes | 13.6777 | 13.6590 | 4.0, 5.3 | yes | 28.0653 | 28.0826 | 7.2, 6.3 | yes | yes | PASS |
| R23 | 10 | yes | 0.8045 | 0.8103 | +0.0058 | 9.0, 2.4 | yes | 0.0824 | 0.0827 | +0.0003 | 4.3, 3.8 | yes | 59.2282 | 59.2301 | 7.4, 7.3 | yes | 91.1411 | 91.3112 | 17.1, 11.7 | yes | yes | PASS |
| R24 | 10 | yes | 0.8101 | 0.8130 | +0.0029 | 7.4, 4.1 | yes | 0.0863 | 0.0862 | -0.0001 | 3.9, 4.0 | yes | 63.8694 | 63.9092 | 12.2, 10.7 | yes | 83.7993 | 83.8560 | 23.6, 20.6 | yes | yes | PASS |
| R25 | 10 | yes | 0.8167 | 0.8178 | +0.0011 | 6.4, 5.2 | yes | 0.0937 | 0.0925 | -0.0012 | 2.9, 4.8 | yes | 70.8033 | 70.8453 | 17.1, 15.2 | yes | 81.5041 | 81.4980 | 30.7, 31.2 | yes | yes | PASS |
| R26 | 5 | yes | 0.8126 | 0.8133 | +0.0007 | 4.3, 3.8 | yes | 0.0853 | 0.0847 | -0.0006 | 2.5, 3.2 | yes | 60.7887 | 60.8639 | 7.2, 5.6 | yes | 86.0825 | 86.1484 | 14.4, 12.3 | yes | yes | PASS |
| R27 | 10 | yes | 0.8145 | 0.8169 | +0.0024 | 7.2, 4.4 | yes | 0.0935 | 0.0908 | -0.0027 | 1.8, 5.9 | yes | 74.2236 | 74.2090 | 20.2, 21.0 | yes | 82.1388 | 82.1524 | 32.0, 31.0 | yes | yes | PASS |
| R28 | 10 | yes | 0.8117 | 0.8128 | +0.0012 | 6.4, 5.1 | yes | 0.0949 | 0.0932 | -0.0016 | 2.6, 5.1 | yes | 66.8567 | 66.7654 | 12.0, 15.7 | yes | 79.6104 | 79.6171 | 27.4, 26.9 | yes | yes | PASS |

Failing designs: 
- R01: complete=TRUE rows_ok=TRUE power_ok=TRUE type1_ok=TRUE ENnull_ok=FALSE ENalt_ok=TRUE power_class_ok=TRUE type1_class_ok=TRUE

Classification sensitivity (reference pooled value within 3 standard errors of its threshold): 
- R23: power is 2.5 SE from 0.80, type I is -14.3 SE from 0.10
