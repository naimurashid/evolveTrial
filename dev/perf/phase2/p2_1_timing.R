# Step P2-1 timing: microseconds per call of the exact single-arm probability for K = 1..8 positive weights at realistic
# shapes (prior alpha = beta = 0.3, events 0 to 20 per interval), and of the between-arm hook, against the phase 1
# Monte Carlo it replaces (5,000 posterior draws + medians). Usage: Rscript dev/perf/phase2/p2_1_timing.R [lib]
a <- commandArgs(TRUE)
lib <- if (length(a) >= 1) a[1] else "dev/perf/lib/dev_p2"
library(evolveTrial, lib.loc = normalizePath(lib))
ns <- function(x) getFromNamespace(x, "evolveTrial")
pmb <- ns(".exact_pr_median_below")
set.seed(3)
reps <- 2000L
out <- NULL
for (K in 1:8) {
  Kt <- K + 1L
  lens <- rep(3, Kt)
  inputs <- replicate(50, list(E = as.integer(sample(0:20, Kt, TRUE)), PT = runif(Kt, 2, 80)), simplify = FALSE)
  thr <- 3 * (K - 1) + 1.7                          # inside interval K, so K positive weights
  al <- rep(0.3, Kt)
  t0 <- proc.time()[["elapsed"]]
  for (r in seq_len(reps / 50)) for (x in inputs) pmb(x$E, x$PT, al, al, lens, thr)
  us <- (proc.time()[["elapsed"]] - t0) / reps * 1e6
  ev <- mean(vapply(inputs, function(x) ns(".exact_pr_median")(x$E, x$PT, al, al, lens, thr)$n_eval, 0))
  t0 <- proc.time()[["elapsed"]]
  for (x in inputs[1:10]) {
    lam <- draw_posterior_hazard_samples_cpp(Kt, x$E, x$PT, al, al, 5000L)
    mean(calculate_median_survival_matrix_cpp(lam, lens) < thr)
  }
  us_mc <- (proc.time()[["elapsed"]] - t0) / 10 * 1e6
  out <- rbind(out, data.frame(K_pos = K, us_per_call_exact = round(us, 1), mean_integrand_evals = round(ev),
                               us_phase1_mc_5000 = round(us_mc)))
}
# R-call overhead of the hook (a trivial thr <= 0 call returns before any work)
t0 <- proc.time()[["elapsed"]]
for (r in 1:reps) pmb(c(1L, 1L), c(5, 5), c(0.3, 0.3), c(0.3, 0.3), c(3, 3), -1)
overhead <- (proc.time()[["elapsed"]] - t0) / reps * 1e6
print(out, row.names = FALSE)
cat(sprintf("R .Call overhead per hook call (thr <= 0, no work): %.1f us (included in the column above)\n", overhead))
lp <- ns(".exact_loghr_probs")
t0 <- proc.time()[["elapsed"]]
for (r in 1:reps) lp(0.1, 0.4, -0.3)
cat(sprintf("between-arm loghr_probs hook: %.1f us per call (incl. overhead)\n",
            (proc.time()[["elapsed"]] - t0) / reps * 1e6))
