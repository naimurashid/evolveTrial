# P3-2: branch counters for gate B, on top of the P3-1 instrumented R body (p31_instrument.R, its counters kept). Exit
# traces on calculate_current_probs_hc and calculate_current_probs_vs_ref read the two probabilities the R body just
# computed and the thresholds it is about to compare them with (IL:147-159, IL:240-253), and record per trial:
#   eff_crossed      pr_eff >= the efficacy threshold, with the R body's guard (rule TRUE, threshold non-NULL and finite)
#   fut_crossed      pr_fut >= the futility threshold (non-NULL and finite)
#   both_held        both at one look (the efficacy stop must win)
#   eff_rule_inert   rule TRUE but the threshold NULL or not finite (the R body skips the check)
# Efficacy and futility stops are counted from the R body's chunk result (sum_stop_eff, sum_stop_fut). The traces draw
# nothing; gate B checks the instrumented run against the plain R body on every chunk.
source(file.path(normalizePath("~/research/evolveTrial-merge"), "dev/perf/phase3/p31_instrument.R"))

.p32_look <- function(args, v, vsref) {
  H <- get(".p31B", envir = globalenv())
  rule <- if (vsref) args$efficacy_stopping_rule_vs_ref else args$efficacy_stopping_rule_hc
  te <- if (vsref) args$efficacy_threshold_vs_ref_prob else args$efficacy_threshold_hc_prob
  tf <- if (vsref) args$futility_threshold_vs_ref_prob else args$futility_threshold_hc_prob
  pre <- if (vsref) "vr_" else "hc_"
  e_on <- isTRUE(rule) && !is.null(te) && is.finite(te)
  if (isTRUE(rule) && !e_on) H$hit(paste0(pre, "eff_rule_inert"))
  e <- e_on && isTRUE(v$pr_eff >= te)
  f <- !is.null(tf) && is.finite(tf) && isTRUE(v$pr_fut >= tf)
  if (e) H$hit(paste0(pre, "eff_crossed"))
  if (f) H$hit(paste0(pre, "fut_crossed"))
  if (e && f) H$hit(paste0(pre, "both_held"))
}

# Run expr with the P3-1 entry trace and the P3-2 exit traces active.
p32_with_trace <- function(expr) {
  ns <- asNamespace("evolveTrial")
  assign(".p32_look", .p32_look, envir = globalenv())
  suppressMessages(trace("calculate_current_probs_hc", where = ns, print = FALSE,
                         tracer = quote(get(".p31B", envir = globalenv())$interim(args)),
                         exit = quote(get(".p32_look", envir = globalenv())(args, returnValue(), FALSE))))
  suppressMessages(trace("calculate_current_probs_vs_ref", where = ns, print = FALSE,
                         exit = quote(get(".p32_look", envir = globalenv())(args, returnValue(), TRUE))))
  on.exit({
    suppressMessages(untrace("calculate_current_probs_hc", where = ns))
    suppressMessages(untrace("calculate_current_probs_vs_ref", where = ns))
  }, add = TRUE)
  force(expr)
}
