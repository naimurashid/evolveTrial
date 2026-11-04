#!/bin/bash
# P3-2: install the merge working tree (tracked + untracked, non-ignored files outside dev/) into dev/perf/lib/p32_$1.
# With a second argument "head", installs a git archive of HEAD instead (the unmodified 8e1430f baseline).
# Same method as dev/perf/phase3/p30_install.sh.
set -e
cd ~/research/evolveTrial-merge
B=dev/perf/build/p32_$1; L=dev/perf/lib/p32_$1
rm -rf "$B" "$L"; mkdir -p "$B" "$L"
if [ "$2" = "head" ]; then
  git archive HEAD | tar -x -C "$B"; rm -rf "$B/dev"
else
  git ls-files -co --exclude-standard | grep -v '^dev/' | tar -T - -c | tar -x -C "$B"
fi
R CMD INSTALL --no-test-load -l "$L" "$B" > dev/perf/phase3/p32_build_$1.log 2>&1
tail -1 dev/perf/phase3/p32_build_$1.log
