# P3-2 gate S sizing, fixed before the full run (plan Section 5, handoff lesson 6). From the pilot runs
# (p32_gateS_pilot.sh) it estimates, per frame, each statistic's spread and each engine's cost per trial, and picks the
# number of R-body and exact trials so that every statistic has test power of at least 0.95 at a true difference of 0
# under the phase 2 margins (power 0.01, type I error 0.005, EN 1% of the R-body mean), with the test of
# p32_gateS_compare.R (two one-sided tests at 0.05 each).
#
# Test power at difference 0 is 1 - 2 * pnorm(z95 - m / se), so >= 0.95 needs m / se >= z95 + z975 = 3.605.
# With n1 R-body and n2 = k * n1 exact trials, se^2 = s2 * (1 / n1 + 1 / n2) = s2 * (1 + 1 / k) / n1, where s2 is
# p (1 - p) for a proportion and the per-trial variance of N for EN. Hence n1 >= s2 (1 + 1 / k) (3.605 / m)^2.
# k = sqrt(c_r / c_exact) (cost per trial from the pilot walls) minimizes the cost for a given se, bounded to [4, 16].
# Spread: the larger of the R-body and exact pilot estimates (proportions: the p closest to 0.5 of the two). Safety
# factor 1.2 on n1 for the uncertainty of the pilot estimates, then n1 rounded up to runs of 5,000 and n2 to runs of
# 50,000 (as in P3-1).
# Usage: Rscript p32_gateS_size.R <pilot dir> <out.tsv> [workers per run = 4] [parallel runs = 4]
a <- commandArgs(TRUE)
W <- if (length(a) >= 3) as.integer(a[3]) else 4L
P <- if (length(a) >= 4) as.integer(a[4]) else 4L
Z <- qnorm(0.95) + qnorm(0.975)
SAFETY <- 1.2; RUN_R <- 5000L; RUN_X <- 50000L
fs <- list.files(a[1], "\\.rds$", full.names = TRUE)
runs <- lapply(fs, readRDS)
exp_arm <- function(tab) { arms <- unique(tab$Arm_Name); if (length(arms) == 1L) arms else "Triplet" }
st <- do.call(rbind, lapply(runs, function(r) {
  t <- r$tab; t <- t[t$Arm_Name == exp_arm(t), ]
  t0 <- t[t$scenario == 1L, ]; t1 <- t[t$scenario == 2L, ]
  data.frame(frame = r$frame, engine = r$engine, n = r$n, wall = r$wall,
             type1 = t0$Type_I_Error_or_Power, power = t1$Type_I_Error_or_Power, EN_null = t0$Exp_N, EN_alt = t1$Exp_N,
             v_null = t0$Var_Exp_N * r$n, v_alt = t1$Var_Exp_N * r$n, pet_eff_alt = t1$PET_Efficacy,
             pet_eff_null = t0$PET_Efficacy, stringsAsFactors = FALSE)
}))
out <- list()
for (fr in sort(unique(st$frame))) {
  r <- st[st$frame == fr & st$engine == "r", ]; x <- st[st$frame == fr & st$engine == "exact", ]
  c_r <- sum(r$wall) / sum(r$n); c_x <- sum(x$wall) / sum(x$n)          # seconds per trial (both scenarios, W workers)
  k <- min(16, max(4, sqrt(c_r / c_x)))
  pp <- function(p) { p <- p[which.min(abs(p - 0.5))]; p * (1 - p) }
  s2 <- c(power = pp(c(r$power, x$power)), type1 = pp(c(r$type1, x$type1)),
          EN_null = max(r$v_null, x$v_null), EN_alt = max(r$v_alt, x$v_alt))
  m <- c(power = 0.01, type1 = 0.005, EN_null = 0.01 * r$EN_null, EN_alt = 0.01 * r$EN_alt)
  n1_need <- s2 * (1 + 1 / k) * (Z / m)^2
  n1 <- ceiling(SAFETY * max(n1_need) / RUN_R) * RUN_R
  n2 <- ceiling(k * n1 / RUN_X) * RUN_X
  se <- sqrt(s2 * (1 / n1 + 1 / n2))
  tp <- pmax(1 - 2 * pnorm(qnorm(0.95) - m / se), 0)
  wall_r <- c_r * n1; wall_x <- c_x * n2                                 # run-seconds (each run uses W workers)
  out[[fr]] <- data.frame(frame = fr, k = round(k, 2), c_r_ms = 1e3 * c_r, c_x_ms = 1e3 * c_x,
    p_power = r$power, p_type1 = r$type1, pet_eff_null_r = r$pet_eff_null, pet_eff_alt_r = r$pet_eff_alt,
    EN_null = r$EN_null, EN_alt = r$EN_alt, sd_EN_null = sqrt(s2[["EN_null"]]), sd_EN_alt = sqrt(s2[["EN_alt"]]),
    n1_need_power = ceiling(n1_need[["power"]]), n1_need_type1 = ceiling(n1_need[["type1"]]),
    n1_need_EN_null = ceiling(n1_need[["EN_null"]]), n1_need_EN_alt = ceiling(n1_need[["EN_alt"]]),
    binding = names(which.max(n1_need)), R_r = n1, R_exact = n2, runs_r = n1 / RUN_R, runs_exact = n2 / RUN_X,
    tp_power = tp[["power"]], tp_type1 = tp[["type1"]], tp_EN_null = tp[["EN_null"]], tp_EN_alt = tp[["EN_alt"]],
    wall_r_h = wall_r / 3600, wall_exact_h = wall_x / 3600, stringsAsFactors = FALSE)
}
tab <- do.call(rbind, out)
write.table(tab, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
print(tab, digits = 4, row.names = FALSE)
tot <- sum(tab$wall_r_h + tab$wall_exact_h)
cat(sprintf("\ntotal run-hours %.2f (each run %d workers); with %d runs in parallel the wall time is about %.1f h\n",
            tot, W, P, tot / P))
