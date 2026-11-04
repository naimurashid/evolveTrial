#!/bin/bash
# p2_selftest.sh -- self-tests of p2_compare.R (step P2-3). Needs no phase 2 library.
#  1 reference vs itself             -> must PASS
#  2 phase 1 layer 2 dev outputs     -> must PASS   (evidence/layer2_l2a/compare.csv, the 75 dev rows)
#  3 power +0.02 on one design       -> must FAIL
#  4 type I +0.01 on one design      -> must FAIL
#  5 EN_null +2% on one design       -> must FAIL (extra)
# Per-trial EN variances for CSV input come from one tiny reference-engine run per design (R = 2000, first seed), built
# once into $T/var.csv through p2_cert.R mode ref (variance attribute x R).
set -uo pipefail
unset CDPATH; H=$(cd "$(dirname "$0")" >/dev/null && pwd); WT=$(cd "$H/../../.." >/dev/null && pwd); EV=$WT/dev/perf/evidence
REF=$EV/acceptance_set_reference_rows.csv; L2A=$EV/layer2_l2a/compare.csv
T=${1:-$EV/phase2_gate/selftest}; mkdir -p "$T/var_rds"
if [ ! -f "$T/var.csv" ]; then
  awk -F, 'NR>1 && !seen[$1]++ {print $1" "$7}' "$REF" | while read -r d s; do
    [ -f "$T/var_rds/$d.rds" ] || Rscript "$H/p2_cert.R" "$d" "$s" 2000 ref "$T/var_rds/$d.rds" > "$T/var_rds/$d.log" 2>&1
  done
  Rscript -e 'f <- list.files("'"$T"'/var_rds", "rds$", full.names = TRUE); d <- do.call(rbind, lapply(f, function(z) { x <- readRDS(z)
    data.frame(design = x$design, EN_null_var = x$EN_null_var * x$R, EN_alt_var = x$EN_alt_var * x$R) })); write.csv(d, "'"$T"'/var.csv", row.names = FALSE)'
fi
perturb() { # <in csv> <out csv> <column> <design> <delta in proportion of R, or relative for EN>
  Rscript -e 'a <- commandArgs(TRUE); x <- read.csv(a[1]); i <- x$design == a[4]
    if (a[3] %in% c("power", "type1")) { cc <- paste0(a[3], "_count"); x[[cc]][i] <- x[[cc]][i] + round(as.numeric(a[5]) * x$R[i]) } else x[[a[3]]][i] <- x[[a[3]]][i] * (1 + as.numeric(a[5]))
    write.csv(x, a[2], row.names = FALSE)' "$@"; }
rc=0; run() { # <name> <expected exit> <args...>
  n=$1; want=$2; shift 2; Rscript "$H/p2_compare.R" "$@" > "$T/$n.out" 2>&1; got=$?
  tag=PASS; [ "$want" -eq 1 ] && tag=FAIL
  if [ "$got" -eq "$want" ] && grep -q "^\[p2_compare\] $tag" "$T/$n.out"; then echo "selftest $n: exit $got (expected $want) OK -- $(grep -h '^\[p2_compare\]' "$T/$n.out")"; else echo "selftest $n: exit $got, expected $want  FAILED"; rc=1; fi; }
run ref_vs_ref 0 "$REF" "$REF" "$T/ref_vs_ref" --var-csv "$T/var.csv"
run l2a_dev_vs_ref 0 "$L2A" "$REF" "$T/l2a_dev_vs_ref" --var-csv "$T/var.csv"
perturb "$REF" "$T/pert_power.csv" power R23 0.02;  run perturb_power_R23 1 "$T/pert_power.csv" "$REF" "$T/perturb_power_R23" --var-csv "$T/var.csv"
perturb "$REF" "$T/pert_type1.csv" type1 R01 0.01;  run perturb_type1_R01 1 "$T/pert_type1.csv" "$REF" "$T/perturb_type1_R01" --var-csv "$T/var.csv"
perturb "$REF" "$T/pert_en.csv" EN_null R02 0.02;   run perturb_ENnull_R02 1 "$T/pert_en.csv" "$REF" "$T/perturb_ENnull_R02" --var-csv "$T/var.csv"
exit $rc
