# bucket.R -- attribute every Rprof sample to one bucket from its full call stack (innermost first).
#   Rscript bucket.R <rprof.out> <out.csv>
a <- commandArgs(TRUE)
L <- readLines(a[1])
L <- L[!grepl("^(#File|line profiling|sample.interval)", L)]
L <- L[nzchar(L)]
stk <- lapply(strsplit(L, " "), function(z) { z <- z[!grepl("^[0-9]+#[0-9]+$", z)]; gsub('"', "", z) })
has <- function(s, f) any(s %in% f)
POST_CPP <- c("draw_posterior_hazard_samples_cpp", "calculate_median_survival_matrix_cpp", "sample_vs_ref_medians_ph_cpp",
              "calculate_median_survival_piecewise_cpp")
POST_TOP <- c("draw_posterior_hazard_samples", "sample_vs_ref_medians", "sample_vs_ref_medians_ph",
              "calculate_current_probs_hc", "calculate_current_probs_vs_ref", "calculate_median_survival_matrix_cpp")
RNG <- c("rgamma", "rnorm", "rexp", "rweibull", "runif", "sample", "sample.int", "rbeta", "rbinom")
SLICE <- c("slice_arm_data_at_time", "get_active_registry", "calculate_interval_metrics_fast")
b <- vapply(stk, function(s) {
  inner <- s[1]
  if (has(s, POST_TOP) || inner %in% POST_CPP) {
    if (inner %in% POST_CPP) return("posterior: C++ call")
    if (inner %in% RNG) return("posterior: R-level RNG (rgamma/rnorm)")
    if (inner == "calculate_median_survival_piecewise" || has(s, "calculate_median_survival_piecewise")) return("posterior: R median_survival_piecewise loop")
    return("posterior: R-level sampler loop/other")
  }
  if (has(s, SLICE)) return("look: data slice (slice_arm_data_at_time / interval metrics)")
  if (has(s, "extract_events_fast")) return("look: extract_events_fast (rebalance)")
  if (has(s, c("interim_check", "interim_check_hc", "interim_check_vs_ref"))) return("look: gates/decision logic")
  if (has(s, "simulate_chunk")) {
    if (inner %in% RNG) return("trial loop: patient RNG (rexp/rweibull/runif/sample)")
    if (has(s, c("[<-.data.frame", "[<-", "[[<-.data.frame", "make_state", "data.frame"))) return("trial loop: registry/state data.frame bookkeeping")
    return("trial loop: other per-trial R code")
  }
  "outside simulate_chunk (setup/aggregation)"
}, "")
tab <- sort(table(b), decreasing = TRUE)
out <- data.frame(bucket = names(tab), samples = as.integer(tab), pct = round(100 * as.integer(tab) / length(b), 1))
# cross-cutting: self (innermost) time of generic R overhead functions, wherever they occur
inner <- vapply(stk, `[`, "", 1)
gen <- c("[<-.data.frame", "data.frame", "[.data.frame", "[[.data.frame", "as.data.frame", "deparse", ".deparseOpts", "make.names",
         "match.arg", "do.call", "rbind", "cbind", "vapply", "lapply", "%in%", "pmin", "pmax", "sys.call", "length", "structure")
cc <- data.frame(bucket = paste0("[self] ", gen), samples = vapply(gen, function(g) sum(inner == g), 0L))
cc$pct <- round(100 * cc$samples / length(b), 1)
# any frame that is rbind / do.call(rbind) inside the trial loop (O(n^2) growth check)
rb <- sum(vapply(stk, function(s) has(s, c("rbind", "rbind.data.frame")) && has(s, "simulate_chunk"), TRUE))
res <- rbind(out, cc[cc$samples > 0, ], data.frame(bucket = "[stack] rbind inside simulate_chunk", samples = rb, pct = round(100 * rb / length(b), 1)))
res$total_samples <- length(b)
write.csv(res, a[2], row.names = FALSE)
print(res, row.names = FALSE)
