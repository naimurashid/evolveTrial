# p2_2_sanity.R -- step P2-2 smoke test (not G5): phase 2 ("exact") against phase 1 ("cpp") on R01 and R23 at R = 2,000,
# one worker, sequential, the capture's cert seed. z = (exact - cpp) / sqrt(var_exact + var_cpp), variances of the mean
# from run_scenarios(return_variance = TRUE). Usage: Rscript dev/perf/phase2/p2_2_sanity.R [lib] [R] [out.csv]
a <- commandArgs(TRUE)
lib <- if (length(a) >= 1) a[1] else "dev/perf/lib/dev_p2"
R <- as.integer(if (length(a) >= 2) a[2] else 2000)
out <- if (length(a) >= 3) a[3] else "dev/perf/phase2/p2_2_sanity.csv"
library(evolveTrial, lib.loc = lib)
options(evolveTrial.looks_stop_at_final = TRUE)
rows <- list()
for (D in c("R01", "R23")) {
  d <- readRDS(sprintf("dev/perf/evidence/args/args_%s.rds", D))
  run <- function(eng) {
    options(evolveTrial.engine = eng)
    ba <- d$base_args; ba$num_simulations <- R; ba$parallel_replicates <- FALSE
    t0 <- proc.time()[["elapsed"]]
    tab <- run_scenarios(ba, d$scens, seed = d$cert_seed, return_variance = TRUE)
    list(tab = as.data.frame(tab), s = proc.time()[["elapsed"]] - t0)
  }
  p1 <- run("cpp"); p2 <- run("exact")
  for (s in seq_along(d$scens)) {
    t1 <- p1$tab[p1$tab$scenario == s, ]; t2 <- p2$tab[p2$tab$scenario == s, ]
    for (j in seq_len(nrow(t1))) {
      for (m in c("Type_I_Error_or_Power", "PET_Futility", "Exp_N")) {
        v1 <- t1[[paste0("Var_", m)]][j]; v2 <- t2[[paste0("Var_", m)]][j]
        rows[[length(rows) + 1L]] <- data.frame(design = D, scenario = s, arm = as.character(t1[[1]][j]), metric = m,
          cpp = t1[[m]][j], exact = t2[[m]][j], se = sqrt(v1 + v2),
          z = if (sqrt(v1 + v2) > 0) (t2[[m]][j] - t1[[m]][j]) / sqrt(v1 + v2) else NA_real_, stringsAsFactors = FALSE)
      }
    }
  }
  rows[[length(rows) + 1L]] <- data.frame(design = D, scenario = NA, arm = "all", metric = "seconds_per_1000_trials",
    cpp = p1$s / (length(d$scens) * R / 1000), exact = p2$s / (length(d$scens) * R / 1000), se = NA, z = NA)
}
res <- do.call(rbind, rows)
write.csv(res, out, row.names = FALSE)
print(res, digits = 4)
