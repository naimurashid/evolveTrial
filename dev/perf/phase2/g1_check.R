# Gate G1 (phase 2 plan dev/plans/2026-11-04-phase2-plan.md Section 3): exact posterior probabilities, step P2-1.
#
# Usage (from the worktree root):  Rscript dev/perf/phase2/g1_check.R [lib] [n_cases] [n_mc_cases] [cores]
#   lib defaults to dev/perf/lib/dev_p2. Writes dev/perf/phase2/g1_results.csv (one row per case and reference) and
#   prints the summary table. Exit status 0 when G1 passes, 1 when it fails.
#
# Single-arm. n_cases (>= 2000) random inputs to the package hook .exact_pr_median, covering 1 to 8 positive weights,
# events 0 to 300 per interval, weak and moderate priors (alpha = beta in 0.05 to 1, as in the acceptance designs, plus
# alpha = 0 cases), and thresholds inside an interval, exactly at a cutpoint, and beyond the last cutpoint. Each case is
# checked against:
#   * the gamma-sum tail recomputed here from the inputs (weights computed independently in R; checks the weight logic);
#   * "pgamma" when one term is left (R's pgamma), "integrate" when two are left (R integrate of the conditional pgamma
#     over the density of one term, split at quantiles), "moschopoulos" for any count where the long double series
#     converges within the term budget (dev/perf/phase2/ref_moschopoulos.cpp);
#   * "mc": for every case without a deterministic reference and n_mc_cases (>= 50) further cases, 10^7 draws of the phase 1 C++ sampler (draw_posterior_hazard_samples_cpp and
#     calculate_median_survival_matrix_cpp), z = (mc - exact) / sqrt(exact (1 - exact) / 10^7).
# Pass: |exact - reference| <= 1e-7 for every deterministic reference, |z| <= 4 for every MC case.
#
# Also: the median identity med < thr  <=>  sum_k w_k lambda_k > log 2 on 10^6 phase 1 draws (zero hazards included),
# and median_piecewise_r (the vs-reference median, via .chunk_ph_sample) against calculate_median_survival_matrix_cpp on
# replayed draws.
#
# Between-arm. (a) on random PH inputs, .exact_pr_loghr's mean and sd are identical to ph_beta_mode_var_cpp and its
# probabilities are identical to pnorm on them; (b) for 50 cases, 10^6 draws of the phase 1 sampler .chunk_ph_sample,
# mean(logHR >= m) and mean(logHR < m) within 4 standard errors; (c) the rnorm edge semantics.

a <- commandArgs(TRUE)
lib <- if (length(a) >= 1) a[1] else "dev/perf/lib/dev_p2"
n_cases <- if (length(a) >= 2) as.integer(a[2]) else 2400L
n_mc <- if (length(a) >= 3) as.integer(a[3]) else 64L
cores <- if (length(a) >= 4) as.integer(a[4]) else 16L
stopifnot(n_cases >= 2000L, n_mc >= 50L)
library(evolveTrial, lib.loc = normalizePath(lib))
stopifnot(dirname(getNamespaceInfo("evolveTrial", "path")) == normalizePath(lib))
ns <- function(x) getFromNamespace(x, "evolveTrial")
Rcpp::sourceCpp("dev/perf/phase2/ref_moschopoulos.cpp")
t_start <- Sys.time()
TOL <- 1e-7
rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- data.frame(..., stringsAsFactors = FALSE)

# ---- input generator -------------------------------------------------------------------------------------------------
gen_case <- function(i) {
  m <- ((i - 1L) %% 8L) + 1L                        # target number of positive weights
  mode <- c("inside", "at", "beyond")[((i - 1L) %/% 8L) %% 3L + 1L]
  Ktot <- if (mode == "beyond") m else m + sample(0:3, 1)
  if (mode == "at" && Ktot == m) Ktot <- m + 1L     # "at" needs the cutpoint that ends interval m
  lens <- if (runif(1) < 0.4) rep(sample(c(1.6, 2, 2.667, 3, 6), 1), Ktot) else runif(Ktot, 0.3, 6)
  start <- Reduce(`+`, lens, accumulate = TRUE)     # repeated double addition, as the engine
  start <- c(0, start[-Ktot])                       # start of interval k
  endm <- start[m] + lens[m]
  thr <- switch(mode,
    inside = if (m < Ktot) start[m] + runif(1, 0.001, 0.999) * lens[m] else start[m] + runif(1, 0.001, 0.999) * lens[m],
    at = endm,
    beyond = endm + runif(1, 0, 24))
  pr <- runif(1, 0.05, 1)
  alpha <- if (runif(1) < 0.15) runif(Ktot, 0.05, 1) else rep(pr, Ktot)
  if (runif(1) < 0.05) alpha[sample(Ktot, 1)] <- 0
  beta <- alpha
  beta[beta == 0] <- pr
  E <- ifelse(runif(Ktot) < 0.3, 0L, as.integer(pmin(300, round(exp(runif(Ktot, 0, log(300)))))))
  h <- exp(runif(Ktot, log(0.01), log(1)))
  PT <- E / h + runif(Ktot, 0, 60)
  PT[runif(Ktot) < 0.05] <- 0
  list(id = i, m = m, mode = mode, K = Ktot, E = E, PT = PT, alpha = alpha, beta = beta, lens = lens, thr = thr)
}

# independent weight computation (R): part of interval k below thr, last interval open-ended
gamma_terms <- function(cs) {
  K <- cs$K
  shape <- cs$alpha + cs$E
  rate <- cs$beta + cs$PT
  rate[rate <= 0] <- 1e-6
  start <- 0
  w <- numeric(K)
  for (k in seq_len(K)) {
    end <- start + cs$lens[k]
    w[k] <- if (cs$thr > start) (if (k < K) min(cs$thr, end) else cs$thr) - start else 0
    start <- end
  }
  keep <- w > 0 & shape > 0
  sc <- w[keep] / rate[keep]
  sh <- shape[keep]
  if (length(sc)) {                                 # merge equal scales (exact equality), as a single gamma
    u <- unique(sc)
    sh <- vapply(u, function(s) sum(sh[sc == s]), 0)
    sc <- u
  }
  list(shape = sh, scale = sc, n_pos = sum(w > 0))
}

# two-term reference: Pr(S > c) = Q_i(c / t_i) + int_0^{c / t_i} f_i(g) Q_j((c - t_i g) / t_j) dg, split at quantiles;
# for a_i < 1 every piece uses v = g^a_i (the density is singular at 0; in v the integrand is smooth).
ref_integrate <- function(sh, sc, c = log(2)) {
  i <- 1L; j <- 2L
  x <- c / sc[i]
  h <- function(g) pgamma((c - sc[i] * g) / sc[j], sh[j], lower.tail = FALSE)
  qs <- qgamma(c(1e-6, 1e-3, 0.05, 0.25, 0.5, 0.75, 0.95, 0.999, 1 - 1e-6), sh[i])
  qj <- (c - sc[j] * qgamma(c(1e-6, 1e-3, 0.05, 0.5, 0.95, 0.999, 1 - 1e-6), sh[j])) / sc[i]
  br <- sort(unique(c(0, qs[qs < x], qj[qj > 0 & qj < x], x)))
  tot <- pgamma(x, sh[i], lower.tail = FALSE)
  for (p in seq_len(length(br) - 1L)) {
    lo <- br[p]; hi <- br[p + 1L]
    if (sh[i] < 1) {
      ai <- sh[i]
      f <- function(v) { g <- v^(1 / ai); exp(-g - lgamma(ai + 1)) * h(g) }
      val <- integrate(f, lo^ai, hi^ai, rel.tol = 1e-12, abs.tol = 1e-15, subdivisions = 5000L, stop.on.error = FALSE)
    } else {
      f <- function(g) dgamma(g, sh[i]) * h(g)
      val <- integrate(f, lo, hi, rel.tol = 1e-12, abs.tol = 1e-15, subdivisions = 5000L, stop.on.error = FALSE)
    }
    if (val$message != "OK") return(NA_real_)
    tot <- tot + val$value
  }
  tot
}

# ---- single-arm: exact vs deterministic references -------------------------------------------------------------------
set.seed(20261104)
cases <- lapply(seq_len(n_cases), gen_case)
pm <- ns(".exact_pr_median")
gst <- ns(".exact_gamma_sum_tail")
sa <- vector("list", n_cases)
for (cs in cases) {
  ex <- pm(cs$E, cs$PT, cs$alpha, cs$beta, cs$lens, cs$thr)
  gt <- gamma_terms(cs)
  stopifnot(gt$n_pos == cs$m)                       # the generator hit its target count of positive weights
  n_t <- length(gt$shape)
  tail_r <- if (n_t == 0) 0 else gst(gt$shape, gt$scale, log(2))$tail
  base <- list(case = cs$id, n_pos = cs$m, n_terms = n_t, thr_mode = cs$mode, exact = ex$below)
  add(check = "sa_weights", case = cs$id, n_pos = cs$m, n_terms = n_t, thr_mode = cs$mode, reference = "R weights",
      exact = ex$below, ref = tail_r, abs_err = abs(ex$below - tail_r), z = NA_real_, info = "",
      pass = isTRUE(abs(ex$below - tail_r) <= 1e-13) && isTRUE(ex$above == 1 - ex$below))
  if (n_t == 1) {
    r <- pgamma(log(2), gt$shape, scale = gt$scale, lower.tail = FALSE)
    add(check = "sa_det", case = cs$id, n_pos = cs$m, n_terms = n_t, thr_mode = cs$mode, reference = "pgamma",
        exact = ex$below, ref = r, abs_err = abs(ex$below - r), z = NA_real_, info = "", pass = abs(ex$below - r) <= TOL)
  }
  if (n_t == 2) {
    r <- ref_integrate(gt$shape, gt$scale)
    if (!is.na(r)) add(check = "sa_det", case = cs$id, n_pos = cs$m, n_terms = n_t, thr_mode = cs$mode,
                       reference = "integrate", exact = ex$below, ref = r, abs_err = abs(ex$below - r), z = NA_real_,
                       info = "", pass = abs(ex$below - r) <= TOL)
  }
  if (n_t >= 2) {
    mo <- moschopoulos_tail(gt$shape, gt$scale, log(2), 4000L, 1e-13)
    if (mo$converged) add(check = "sa_det", case = cs$id, n_pos = cs$m, n_terms = n_t, thr_mode = cs$mode,
                          reference = "moschopoulos", exact = ex$below, ref = mo$tail, abs_err = abs(ex$below - mo$tail),
                          z = NA_real_, info = sprintf("terms=%d", mo$terms), pass = abs(ex$below - mo$tail) <= TOL)
  }
  sa[[cs$id]] <- c(base, list(method = ex$method, n_eval = ex$n_eval, err_est = ex$err_est, round_est = ex$round_est,
                              converged = ex$converged))
}
sa_df <- do.call(rbind, lapply(sa, as.data.frame))
det <- do.call(rbind, rows)
has_det <- unique(det$case[det$check == "sa_det"])
cat("single-arm cases:", n_cases, " with a deterministic reference:", length(has_det),
    " contour cases:", sum(sa_df$method == 2), " max contour err_est:", signif(max(sa_df$err_est), 3),
    " max round_est:", signif(max(sa_df$round_est), 3), " all converged:", all(sa_df$converged), "\n")

# ---- single-arm: Monte Carlo with the phase 1 sampler ----------------------------------------------------------------
# Every case without a deterministic reference, plus n_mc / 8 further cases per number of positive weights (preferring a
# non-degenerate probability).
draw <- evolveTrial::draw_posterior_hazard_samples_cpp
medm <- evolveTrial::calculate_median_survival_matrix_cpp
pick <- integer(0)
for (m in 1:8) {
  cand <- sa_df$case[sa_df$n_pos == m]
  no_det <- setdiff(cand, has_det)
  mid <- cand[sa_df$exact[cand] > 1e-3 & sa_df$exact[cand] < 1 - 1e-3]
  pool <- unique(c(sample(mid), sample(cand)))
  pick <- c(pick, no_det, head(setdiff(pool, no_det), ceiling(n_mc / 8)))   # every case without one, plus n_mc / 8
}
NMC <- 1e7
mc_one <- function(id) {
  cs <- cases[[id]]
  set.seed(1e5 + id)
  cnt <- 0
  for (b in 1:10) {
    lam <- draw(cs$K, as.integer(cs$E), cs$PT, cs$alpha, cs$beta, NMC / 10)
    med <- medm(lam, cs$lens)
    cnt <- cnt + sum(med < cs$thr)
  }
  cnt / NMC
}
mcp <- unlist(parallel::mclapply(pick, mc_one, mc.cores = cores, mc.preschedule = FALSE))
for (k in seq_along(pick)) {
  id <- pick[k]; e <- sa_df$exact[id]
  se <- sqrt(e * (1 - e) / NMC)
  z <- if (se > 0) (mcp[k] - e) / se else if (mcp[k] == e) 0 else Inf
  add(check = "sa_mc", case = id, n_pos = sa_df$n_pos[id], n_terms = sa_df$n_terms[id], thr_mode = sa_df$thr_mode[id],
      reference = "mc_1e7", exact = e, ref = mcp[k], abs_err = abs(mcp[k] - e), z = z, info = sprintf("se=%.3g", se),
      pass = abs(z) <= 4)
}

# ---- median identity on phase 1 draws ---------------------------------------------------------------------------------
# med < thr  <=>  H(thr) > log 2, with zero hazards (alpha = 0, E = 0) and thresholds inside, at and beyond cutpoints.
set.seed(77)
id_bad <- 0L; id_tot <- 0L; id_tie <- 0L
for (r in 1:20) {
  K <- sample(2:9, 1)
  lens <- runif(K, 0.5, 5)
  alpha <- runif(K, 0.05, 1); alpha[sample(K, 2, replace = TRUE)] <- 0
  E <- ifelse(runif(K) < 0.5, 0L, sample(1:50, K, replace = TRUE))
  lam <- draw(K, as.integer(E), runif(K, 0, 40), alpha, runif(K, 0.05, 1), 50000)
  med <- medm(lam, lens)
  start <- c(0, Reduce(`+`, lens, accumulate = TRUE)[-K])
  for (thr in c(runif(3, 0, sum(lens) + 10), start[-1], sum(lens) + 5)) {
    w <- ifelse(thr > start, ifelse(seq_len(K) < K, pmin(thr, start + lens), thr) - start, 0)
    H <- drop(lam %*% w)
    tie <- abs(H - log(2)) < 1e-12 * pmax(1, H) | abs(med - thr) < 1e-12 * max(1, thr)
    id_tie <- id_tie + sum(tie)
    id_bad <- id_bad + sum(((med < thr) != (H > log(2)))[!tie])
    id_tot <- id_tot + sum(!tie)
  }
}
add(check = "median_identity", case = NA, n_pos = NA, n_terms = NA, thr_mode = "", reference = "phase1 medians",
    exact = NA, ref = NA, abs_err = id_bad, z = NA, info = sprintf("rows=%d near_ties_excluded=%d", id_tot, id_tie),
    pass = id_bad == 0)

# median_piecewise_r (vs-reference path) against calculate_median_survival_matrix_cpp on replayed draws, including a
# zero last hazard (shape 0) so the open-ended tail and the Inf branch are exercised.
phs <- ns(".chunk_ph_sample")
mp_bad <- 0L; mp_tot <- 0L
for (r in 1:10) {
  K <- sample(3:12, 1)
  lens <- runif(K, 0.5, 4)
  alpha <- rep(runif(1, 0.05, 1), K)
  EC <- sample(0:20, K, TRUE); ET <- sample(0:20, K, TRUE)
  if (r %% 2 == 0) { alpha[K] <- 0; EC[K] <- 0L; ET[K] <- 0L }
  PTC <- runif(K, 1, 50); PTT <- runif(K, 1, 50); beta <- alpha; beta[beta == 0] <- 0.3
  n <- 4000
  set.seed(500 + r)
  o <- phs(as.integer(EC), PTC, as.integer(ET), PTT, alpha, beta, lens, 0, 1, 1e-6, 50L, n)
  set.seed(500 + r)
  bp <- evolveTrial::ph_beta_mode_var_cpp(as.integer(EC), PTC, as.integer(ET), PTT, alpha, beta, 0, 1, 1e-6, 50L)
  lhr <- rnorm(n, bp$mean, bp$sd)
  shapes <- alpha + EC + ET
  lam <- t(vapply(seq_len(n), function(i) rgamma(K, shapes, rate = (beta + PTC) + exp(lhr[i]) * PTT), numeric(K)))
  mc <- medm(lam, lens)
  mt <- medm(lam * exp(lhr), lens)
  mp_bad <- mp_bad + sum(!(identical(o$logHR, lhr))) + sum(o$medCtrl != mc) + sum(o$medTrt != mt)
  mp_tot <- mp_tot + 2 * n
}
add(check = "median_functions_agree", case = NA, n_pos = NA, n_terms = NA, thr_mode = "",
    reference = "median_piecewise_r vs calculate_median_survival_matrix_cpp", exact = NA, ref = NA, abs_err = mp_bad,
    z = NA, info = sprintf("medians=%d", mp_tot), pass = mp_bad == 0)

# ---- between-arm -------------------------------------------------------------------------------------------------------
plr <- ns(".exact_pr_loghr")
set.seed(4242)
gen_ph <- function() {
  K <- sample(c(8, 9, 12, 13, 14, 15), 1)
  alpha <- rep(runif(1, 0.05, 1), K)
  EC <- as.integer(ifelse(runif(K) < 0.3, 0, sample(0:40, K, TRUE)))
  ET <- as.integer(ifelse(runif(K) < 0.3, 0, sample(0:40, K, TRUE)))
  list(EC = EC, PTC = EC / runif(K, 0.02, 0.5) + runif(K, 0, 30), ET = ET, PTT = ET / runif(K, 0.02, 0.5) + runif(K, 0, 30),
       alpha = alpha, beta = alpha, lens = rep(24 / K, K), m = c(log(1 - runif(1, 0, 0.3)), log(runif(1, 0.6, 1.2)), 0))
}
ph_cases <- replicate(500, gen_ph(), simplify = FALSE)
ba_bad <- 0L
for (p in ph_cases) {
  ex <- plr(p$EC, p$PTC, p$ET, p$PTT, p$alpha, p$beta, p$lens, 0, 1, 1e-6, 50L, p$m)
  bp <- evolveTrial::ph_beta_mode_var_cpp(p$EC, p$PTC, p$ET, p$PTT, p$alpha, p$beta, 0, 1, 1e-6, 50L)
  ok <- identical(ex$mean, bp$mean) && identical(ex$sd, bp$sd) &&
    identical(ex$pr_ge, pnorm(p$m, bp$mean, bp$sd, lower.tail = FALSE)) &&
    identical(ex$pr_lt, pnorm(p$m, bp$mean, bp$sd))
  ba_bad <- ba_bad + !ok
}
add(check = "ba_closed_form", case = NA, n_pos = NA, n_terms = NA, thr_mode = "", reference = "pnorm on Laplace output",
    exact = NA, ref = NA, abs_err = ba_bad, z = NA, info = sprintf("cases=%d x 3 margins", length(ph_cases)),
    pass = ba_bad == 0)
NBA <- 1e6
ba_mc <- parallel::mclapply(1:50, function(i) {
  p <- ph_cases[[i]]
  set.seed(9000 + i)
  o <- phs(p$EC, p$PTC, p$ET, p$PTT, p$alpha, p$beta, p$lens, 0, 1, 1e-6, 50L, NBA)
  ex <- plr(p$EC, p$PTC, p$ET, p$PTT, p$alpha, p$beta, p$lens, 0, 1, 1e-6, 50L, p$m)
  data.frame(i = i, k = seq_along(p$m), ge = ex$pr_ge, lt = ex$pr_lt,
             mc_ge = vapply(p$m, function(m) mean(o$logHR >= m), 0), mc_lt = vapply(p$m, function(m) mean(o$logHR < m), 0))
}, mc.cores = cores, mc.preschedule = FALSE)
for (d in ba_mc) for (r in seq_len(nrow(d))) for (side in c("ge", "lt")) {
  e <- d[[side]][r]; mcv <- d[[paste0("mc_", side)]][r]
  se <- sqrt(e * (1 - e) / NBA)
  z <- if (se > 0) (mcv - e) / se else if (mcv == e) 0 else Inf
  add(check = "ba_mc", case = d$i[r], n_pos = NA, n_terms = NA, thr_mode = sprintf("margin%d_%s", d$k[r], side),
      reference = "mc_1e6", exact = e, ref = mcv, abs_err = abs(mcv - e), z = z, info = sprintf("se=%.3g", se),
      pass = abs(z) <= 4)
}
lp <- ns(".exact_loghr_probs")
edge <- rbind(lp(NaN, 1, 0), lp(0, NaN, 0), lp(0, Inf, 0), lp(0, -1, 0), lp(0.2, 0, c(0.1, 0.2, 0.3)),
              lp(Inf, 1, 5), lp(-Inf, 1, 5), lp(0, 1, NaN), lp(0, 1, c(-Inf, Inf)))
want <- rbind(c(NaN, NaN), c(NaN, NaN), c(NaN, NaN), c(NaN, NaN), c(1, 0), c(1, 0), c(0, 1), c(1, 0), c(0, 1),
              c(NaN, NaN), c(1, 0), c(0, 0))
want[12, ] <- c(0, 1)
edge_ok <- isTRUE(all.equal(unname(edge), want))
add(check = "ba_edges", case = NA, n_pos = NA, n_terms = NA, thr_mode = "", reference = "rnorm semantics",
    exact = NA, ref = NA, abs_err = NA, z = NA, info = "", pass = edge_ok)

# ---- output -------------------------------------------------------------------------------------------------------------
res <- do.call(rbind, rows)
write.csv(res, "dev/perf/phase2/g1_results.csv", row.names = FALSE)
summ <- do.call(rbind, lapply(split(res, paste(res$check, res$reference)), function(d) data.frame(
  check = d$check[1], reference = d$reference[1], rows = nrow(d),
  cases_n_pos = if (all(is.na(d$n_pos))) "" else paste(range(d$n_pos, na.rm = TRUE), collapse = "-"),
  max_abs_err = if (all(is.na(d$abs_err))) NA else signif(max(d$abs_err, na.rm = TRUE), 3),
  max_abs_z = if (all(is.na(d$z))) NA else round(max(abs(d$z), na.rm = TRUE), 2),
  failed = sum(!d$pass))))
print(summ, row.names = FALSE)
for (m in 1:8) {
  d <- res[res$n_pos %in% m & res$check %in% c("sa_det", "sa_mc"), ]
  cat(sprintf("n_pos=%d: det rows %d (max err %.2e), mc rows %d (max |z| %.2f)\n", m, sum(d$check == "sa_det"),
              max(c(0, d$abs_err[d$check == "sa_det"])), sum(d$check == "sa_mc"),
              max(c(0, abs(d$z[d$check == "sa_mc"])))))
}
cat("thr modes with a deterministic reference:", paste(names(table(det$thr_mode[det$check == "sa_det"])),
    table(det$thr_mode[det$check == "sa_det"]), collapse = ", "), "\n")
cat("cases without a deterministic reference:", n_cases - length(has_det), "(covered by MC:",
    length(intersect(setdiff(seq_len(n_cases), has_det), pick)), ")\n")
pass <- all(res$pass)
cat("G1", if (pass) "PASS" else "FAIL", " wall seconds", round(as.numeric(difftime(Sys.time(), t_start, units = "secs"))),
    "\n")
quit(status = if (pass) 0 else 1)
