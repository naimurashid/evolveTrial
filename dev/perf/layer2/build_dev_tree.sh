#!/bin/bash
# build_dev_tree.sh -- assemble the layer 2 engine tree from one commit of this branch (engine C++ port, step 5;
# dev/plans/2026-11-03-step5-layer2.md, Section 2).
#   build_dev_tree.sh <evolveTrial git dir> <commit> <out tree>
# Output: git archive of <commit> with two DESCRIPTION fields appended.
#   P6RecertEvolveTrialSHA  7565a17 (the P6 build field). The frozen P6 protocol (p6r_protocol.R, mode p6) refuses any
#                           library without it. The branch descends from 7565a17 through the M4 commits 0d5546f and
#                           58a4850, so the field states the build's lineage. Asserted below, not assumed.
#   EngineCppSHA            the full commit id. layer2_cert.R requires it to equal L2_ENGINE_SHA in dev mode.
# VerifyM4PatchMD5 and VerifyM4DriverMD5 are NOT written: this tree is not 7565a17 plus the patch, and its driver differs
# from the M4 driver by the seam and the dispatcher.
# Install: R CMD INSTALL --no-multiarch -l <lib> <out tree>. Writes only <out tree>.
set -euo pipefail
SRC=${1:?git dir}; REF=${2:?commit}; OUT=${3:?out tree}
P6_SHA=7565a17c625573d2a697d0d4a85efb2fa6c58ad9
M4_SHA=58a4850                                       # engine of record (verify/p6-nodropout)
[ ! -e "$OUT" ] || { echo "refusing: $OUT exists" >&2; exit 2; }
SHA=$(git -C "$SRC" rev-parse --verify "$REF^{commit}")
git -C "$SRC" merge-base --is-ancestor "$P6_SHA" "$SHA" || { echo "$SHA does not descend from $P6_SHA" >&2; exit 3; }
git -C "$SRC" merge-base --is-ancestor "$M4_SHA" "$SHA" || { echo "$SHA does not descend from $M4_SHA" >&2; exit 3; }
mkdir -p "$OUT"
git -C "$SRC" archive "$SHA" | tar -x -C "$OUT"
grep -q '^P6RecertEvolveTrialSHA:' "$OUT/DESCRIPTION" && { echo "DESCRIPTION already carries a P6 build field" >&2; exit 3; }
grep -q '^Version: 0.0.0.9001$' "$OUT/DESCRIPTION" || { echo "DESCRIPTION Version is not 0.0.0.9001" >&2; exit 3; }
[ -f "$OUT/src/simulate_chunk.cpp" ] && [ -f "$OUT/R/simulate_chunk_cpp.R" ] || { echo "no C++ chunk engine in $SHA" >&2; exit 3; }
printf 'P6RecertEvolveTrialSHA: %s\nEngineCppSHA: %s\n' "$P6_SHA" "$SHA" >> "$OUT/DESCRIPTION"
echo "built $OUT from $SHA (P6 lineage $P6_SHA, M4 $M4_SHA)"
