#!/bin/bash
# P3-2 mutation check of gate B's sensitivity (scratch build dev/perf/lib/p32_mut, working tree untouched): in the phase 1
# engine the HC efficacy check uses the final success value instead of the interim one (m) or the vs-ref efficacy check is skipped (v).
# Gate B on C03/C07 (HC) or D10/D03 (vs-ref) must then report mismatches.
set -e
cd ~/research/evolveTrial-merge
M=$1; B=dev/perf/build/p32_mut_$M; L=dev/perf/lib/p32_mut_$M
rm -rf "$B" "$L"; mkdir -p "$B" "$L"
git ls-files -co --exclude-standard | grep -v '^dev/' | tar -T - -c 2>/dev/null | tar -x -C "$B" || true
if [ "$M" = "m" ]; then
  python3 - "$B/src/simulate_chunk.cpp" <<'PY'
import sys; p=sys.argv[1]; s=open(p).read()
old="tally_mean(tally(med, c.success_thr_interim[a], false))"; assert s.count(old)==1
s=s.replace(old,"tally_mean(tally(med, c.success_thr[a], false))"); open(p,'w').write(s)
PY
else
  python3 - "$B/src/simulate_chunk.cpp" <<'PY'
import sys; p=sys.argv[1]; s=open(p).read()
old="if (pr_eff >= c.vr_eff_prob_interim) {"; assert s.count(old)==1
s=s.replace(old,"if (false && pr_eff >= c.vr_eff_prob_interim) {"); open(p,'w').write(s)
PY
fi
R CMD INSTALL --no-test-load -l "$L" "$B" > dev/perf/phase3/p32_build_mut_$M.log 2>&1
tail -1 dev/perf/phase3/p32_build_mut_$M.log
