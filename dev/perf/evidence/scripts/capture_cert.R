# capture_cert.R -- recapture one certified design's run_scenarios call exactly as the certify phase made it
# (plan dev/plans/2026-11-03-engine-cpp-plan.md, step 2a, task 2a.2; inventory items I1-I6).
#
#   Rscript capture_cert.R <design R01|R02|R03|R04|R23|R24|R25|R26|R27|R28> <R> <out.rds>
#
# The set-up is copied, not paraphrased, from the certify runners. Line numbers below cite
#   RT1 = ~/Downloads/atbo-longleaf/pilot/retune/R/retune_run.R        (md5 6e2a3905c9827523ace6dadb57467027; tier 1:
#         R01 R02 R03 R23 R24 R25 R27 R28)
#   RT2 = ~/Downloads/atbo-longleaf/pilot/retune2a/R/retune2a_run.R    (md5 b388c79f84c84d8c3fe8bc7f64bc89f4; tier 2a:
#         R04 R26). RT2 is RT1 plus the hunks marked "T2A:" (task table, held/protocol check, weights table, minimax ev
#         bound, staged n_init check, base_design id). R04 and R26 were certified by RT2, not RT1.
# Both md5 values equal the ones recorded in the Longleaf certify manifests (~/pilot_r1/retune_t1/c1a/certify_1 and
# ~/pilot_r1/retune_t2a/c2a_r/certify_51 manifest.txt); the script refuses to run if either runner file has changed.
# Local differences from the runner (each one is in capture_report.md): command-line parsing (design + R instead of
# --phase certify --design --seed --select), the frozen cert_design md5 is a literal (results/ holds no cert_design.md5),
# the runner's log files go to a temporary OUT, BATON is loaded from ~/R_libs_frozen/jasa_r1_v040, and a capture hook is
# installed innermost on run_scenarios before the two engine switches.
# Read-only on every repo; writes only <out.rds> and a temporary directory.
`%||%` <- function(x, y) if (is.null(x)) y else x                                                  # RT1:38
a <- commandArgs(TRUE)
if (length(a) != 3L) stop("usage: Rscript capture_cert.R <design> <R> <out.rds>")
DESIGN <- a[1]; RCAP <- as.integer(a[2]); OUTF <- a[3]
stopifnot(DESIGN %in% c("R01", "R02", "R03", "R04", "R23", "R24", "R25", "R26", "R27", "R28"), !is.na(RCAP), RCAP >= 1L)
if (file.exists(OUTF)) stop("output exists: ", OUTF)
PHASE <- "certify"; WORKERS <- 19L                                                                  # RT1:41 (--workers 19)
BATON_VERSION <- "0.4.0"                                                                            # RT1:45
COERCE <- TRUE                                                                                      # RT1:47 (default 1)
TEST_R <- if (RCAP == 10000L) 0L else RCAP                                                          # RT1:48 (--test-R)

# ---- paths (the runner's environment variables, RT1:12-19 and :50-53, resolved to the local copies) ---------------
WT <- normalizePath("~/research/evolveTrial-perf")
AL <- normalizePath("~/Downloads/atbo-longleaf")
TIER <- if (DESIGN %in% c("R04", "R26")) "t2a" else "t1"
RUNNER <- file.path(AL, if (TIER == "t2a") "pilot/retune2a/R/retune2a_run.R" else "pilot/retune/R/retune_run.R")
RUNNER_MD5 <- c(t1 = "6e2a3905c9827523ace6dadb57467027", t2a = "b388c79f84c84d8c3fe8bc7f64bc89f4")[[TIER]]
RR  <- file.path(AL, if (TIER == "t2a") "pilot/retune2a" else "pilot/retune")                      # RETUNE_ROOT
F0  <- normalizePath("~/Downloads/atbo-f0")                                                        # RETUNE_F0_DIR
LIB <- normalizePath(file.path(WT, "dev/perf/lib/ref"))                                            # RETUNE_EXPECTED_LIB
P6R <- normalizePath(file.path(AL, "pilot/p6recert/R/p6r_protocol.R"))                             # RETUNE_P6R_PROTOCOL
M4S <- normalizePath(file.path(AL, "pilot/verify/R/verify_m4_switch.R"))                           # RETUNE_M4_SWITCH
BL  <- normalizePath("~/R_libs_frozen/jasa_r1_v040")                                               # RETUNE_BATON_LIB
SEL <- file.path(AL, if (TIER == "t2a") "results/retune_t2a" else "results/retune_t1")             # --select
SEL_MD5 <- c(t1 = "1e1997bbfa44b232fced57f9d619d533", t2a = "ec1c987c6b7d80eb7ea25f071d0e1882")[[TIER]]
OUT <- tempfile(paste0("capture_", DESIGN, "_")); dir.create(OUT, recursive = TRUE); OUT <- normalizePath(OUT)

# ---- inputs: staged md5, F0 files md5, git check (RT1:59-67, RT2:63-71) -------------------------------------------
md5f <- function(f) unname(tools::md5sum(f))                                                        # RT1:59
if (md5f(RUNNER) != RUNNER_MD5) stop("runner ", RUNNER, " changed since the certify runs")
if (md5f(P6R) != "52d2b43cf69feccffeca53a9924d4b80" || md5f(M4S) != "405f7d60e1a4be532821e3ba603227bc")
  stop("engine switch files differ from the certify manifests")
chk_list <- function(listing, base) { x <- read.table(listing, stringsAsFactors = FALSE, col.names = c("md5", "path"))   # RT1:60-62
  bad <- x$path[md5f(file.path(base, x$path)) != x$md5 | !file.exists(file.path(base, x$path))]
  if (length(bad)) stop("md5 mismatch: ", paste(bad, collapse = ", ")) }
chk_list(file.path(RR, "staged", "STAGED.md5"), file.path(RR, "staged"))                           # RT1:63
chk_list(file.path(RR, "staged", "f0_files.md5"), F0)                                              # RT1:64
F0_SHA <- "1863cecb928eca56d4e52fe202542630568a5f35"                                               # RT1:65-67
git_sha <- tryCatch(suppressWarnings(system2("git", c("-C", F0, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)), error = function(e) character(0))
if (length(git_sha) && !identical(git_sha, F0_SHA)) stop("F0 worktree at ", git_sha, ", not 1863cecb")
TK <- utils::read.delim(file.path(RR, "staged", if (TIER == "t2a") "t2a_tasks.tsv" else "t1_tasks.tsv"),   # RT1:68, RT2:72
                        stringsAsFactors = FALSE, na.strings = "")

# ---- which task (RT1:71-81; certify branch). The cert seed is the design's first certification seed. ---------------
cd <- utils::read.csv(file.path(SEL, "cert_design.csv"), stringsAsFactors = FALSE)                 # RT1:73
if (md5f(file.path(SEL, "cert_design.csv")) != SEL_MD5) stop("cert_design.csv md5 != frozen")       # RT1:74 (literal md5)
CSEED <- min(cd$cert_seed[cd$design == DESIGN])
cr <- cd[cd$design == DESIGN & cd$cert_seed == CSEED, ]                                             # RT1:75-77
if (nrow(cr) != 1L) stop("design ", DESIGN, " seed ", CSEED, " is not exactly one row of the frozen cert_design.csv")
TASK <- cr$task
tr <- TK[TK$task == TASK, ]                                                                         # RT1:79-81
if (nrow(tr) != 1L) stop("unknown task ", TASK)
SID <- tr$scenario_id; PROT <- tr$protocol; SEED <- as.integer(tr$seed)
if (TIER == "t2a") {                                                                                # RT2:86-88
  if (identical(as.integer(tr$held), 1L)) stop("task ", TASK, " is held: ", tr$held_reason)
  if (!PROT %in% c("rws", "xphil")) stop("tier 2a protocol ", PROT, " not expected")
}

# ---- campaign settings (RT1:84-100), set before any package load. BO_SEED is the TASK seed, as in the runner. -----
stray <- grep("^BO_", names(Sys.getenv()), value = TRUE)
if (length(stray)) stop("BO_* variables already set (fail closed): ", paste(stray, collapse = ", "))
SLW <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (nzchar(SLW) && !identical(as.integer(SLW), WORKERS)) stop("SLURM_CPUS_PER_TASK ", SLW, " != --workers ", WORKERS)
ENV <- c(BO_LOW_REPS = "3000", BO_MED_REPS = "5000", BO_HIGH_REPS = "10000",
         BO_BUDGET_S1 = "700", BO_BUDGET_S2 = "300", BO_BUDGET_S3 = "200", BO_N_INIT = "90",
         BO_USE_BUFFER = "TRUE", BO_CONSTRAINT_BUFFER = "0.02",
         BO_VALIDATION_TIER_SIZES = "5,5", BO_VALIDATION_MIN_STOP = "10", BATON_DISABLE_KSNAP = "1",
         BO_SEED = as.character(SEED), SLURM_CPUS_PER_TASK = as.character(WORKERS))
do.call(Sys.setenv, as.list(ENV))                                                                   # RT1:100
if (length(setdiff(loadedNamespaces(), c("compiler", "graphics", "utils", "grDevices", "stats", "datasets", "methods", "base", "tools"))))
  stop("packages loaded before the environment block: ", paste(loadedNamespaces(), collapse = " "))

# ---- libraries (RT1:105-115) ---------------------------------------------------------------------------------------
.libPaths(unique(c(LIB, .libPaths())))
loadNamespace("evolveTrial", lib.loc = LIB)
if (dirname(normalizePath(getNamespaceInfo("evolveTrial", "path"))) != LIB) stop("evolveTrial not from ", LIB)
suppressPackageStartupMessages({ library(tidyverse); if (nzchar(BL)) library(BATON, lib.loc = BL) else library(BATON) })
if (!identical(as.character(utils::packageVersion("BATON")), BATON_VERSION)) stop("BATON is ", utils::packageVersion("BATON"), ", not ", BATON_VERSION)
if (nzchar(BL) && dirname(normalizePath(getNamespaceInfo("BATON", "path"))) != normalizePath(BL)) stop("BATON not from ", BL)
if (dirname(normalizePath(getNamespaceInfo("evolveTrial", "path"))) != LIB) stop("evolveTrial reloaded from another library")
ED <- utils::packageDescription("evolveTrial", lib.loc = LIB)
if (!identical(ED$VerifyM4PatchMD5, "83de2c8b2de1902630eddf00c334b508") || !identical(ED$P6RecertEvolveTrialSHA, "7565a17c625573d2a697d0d4a85efb2fa6c58ad9"))
  stop("library ", LIB, " is not the M4 build (P6 7565a17 + m4_nodropout.patch)")
# plan Section 0 assertions (in addition to the runner's)
if (dirname(getNamespaceInfo("evolveTrial", "path")) != normalizePath(LIB)) stop("plan assertion: evolveTrial path")
if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection"))) stop("RNGkind is ", paste(RNGkind(), collapse = "/"))

# ---- F0 code in the protocol's source order, and the row parsed as that runner did (RT1:118-130, RT2:125-137) ------
SRC <- list(rws = c("R/scenario_parser.R", "R/warmstart_core.R", "R/scenario_configs.R", "R/warmstart_wrappers.R", "scripts/sim_fun_evolveTrial_unified.R"),
            bws = c("R/scenario_parser.R", "R/warmstart_wrappers.R", "R/warmstart_core.R", "scripts/sim_fun_evolveTrial_unified.R"),
            xphil = c("R/warmstart_wrappers.R", "R/warmstart_core.R", "R/scenario_parser.R", "scripts/sim_fun_evolveTrial_unified.R"))
old <- setwd(F0)
suppressMessages(suppressWarnings(for (f in SRC[[PROT]]) source(f, local = FALSE)))
CSV <- file.path(F0, tr$csv)
spec <- switch(PROT,
  rws = { sc <- suppressMessages(parse_warmstart_scenarios(csv_path = CSV, filter_enabled = TRUE, verbose = FALSE)); sc[[SID]] },
  bws = { d <- readr::read_csv(CSV, show_col_types = FALSE); suppressMessages(parse_single_scenario(dplyr::filter(d, scenario_id == SID), verbose = FALSE)) },
  xphil = { d <- readr::read_csv(CSV, show_col_types = FALSE); r <- dplyr::filter(d, scenario_id == !!SID)
            if (nrow(r) != 1L) stop("row not unique"); suppressMessages(parse_single_scenario(r, verbose = FALSE)) })
setwd(old)
if (is.null(spec) || !identical(spec$scenario_id, SID)) stop("scenario ", SID, " not parsed")

# ---- assertions (RT1:134-165; tier 2a weights and minimax bound from RT2:144-162, staged n_init RT2:182) -----------
if (TIER == "t1") {                                                                                 # RT1:134-138
  WEIGHTS <- rbind(cohortA_H0_fut = c(0.01, 1, 0), cohortA_admissible_fut = c(0.5, 1, 0), cohortA_minimax_fut = c(0.99, 0.5, 0.5),
                   doublet_triplet_H0_fut = c(0.01, 1, 0), doublet_triplet_admissible_fut = c(0.5, 1, 0),
                   doublet_triplet_minimax_fut = c(0.99, 0.5, 0.5), doublet_triplet_admissible_N_fut = c(0.75, 1, 0),
                   doublet_triplet_balanced_fut = c(0.01, 0.5, 0.5))
  colnames(WEIGHTS) <- c("w_N", "w_null", "w_alt")
} else {                                                                                            # RT2:144-156
  WP <- rbind(H0 = c(0.01, 1, 0), minimax = c(0.99, 0.5, 0.5), balanced = c(0.01, 0.5, 0.5), H1 = c(0.01, 0, 1),
              admissible = c(0.5, 1, 0), admissible_EN = c(0.25, 1, 0), admissible_N = c(0.75, 1, 0))
  PH <- c(cohortA_admissible_EN_fut = "admissible_EN", cohortA_admissible_N_fut = "admissible_N", cohortA_H1_fut = "H1",
          cohortA_balanced_fut = "balanced", cohortB_H0_bidir = "H0", cohortB_minimax_bidir = "minimax", cohortB_balanced_bidir = "balanced",
          cohortB_H1_bidir = "H1", cohortB_H0_fut = "H0", cohortB_minimax_fut = "minimax", cohortB_balanced_fut = "balanced", cohortB_H1_fut = "H1",
          doublet_triplet_admissible_EN_fut = "admissible_EN", doublet_triplet_H1_fut = "H1", doublet_triplet_gen_optimal = "H0",
          doublet_triplet_gen_minimax = "minimax", doublet_triplet_gen_fleming = "balanced", doublet_triplet_gen_alt_optimal = "H1",
          doublet_triplet_gen_admissible = "admissible", doublet_triplet_gen_admissible_EN = "admissible_EN",
          doublet_triplet_gen_admissible_N = "admissible_N", doublet_triplet_H0_bidir = "H0", doublet_triplet_minimax_bidir = "minimax",
          doublet_triplet_balanced_bidir = "balanced", doublet_triplet_H1_bidir = "H1")
  WEIGHTS <- WP[PH, , drop = FALSE]; rownames(WEIGHTS) <- names(PH)
  if (!SID %in% rownames(WEIGHTS)) stop("no expected weights for ", SID)
  if (!isTRUE(all.equal(unname(WEIGHTS[SID, ]), c(tr$w_N_exp, tr$w_null_exp, tr$w_alt_exp), tolerance = 1e-12))) stop("staged weights differ from the literal table for ", SID)
  colnames(WEIGHTS) <- c("w_N", "w_null", "w_alt")   # RT2 keeps RT1:138 after the table
}
BOUNDS <- list(A = list(total_n = c(30, 120), eff = c(0.85, 0.99), fut = c(0.03, 0.65), ev = c(3, 50), beat = c(1, 16),   # RT1:139-142
                        margin = c(0, 1), prior_strength = c(0.15, 0.6), n_intervals = c(3, 10)),
               B = list(total_n = c(80, 170), eff = c(0.95, 0.99), fut = c(0.2, 0.8), ev = c(5, 25), beat = c(1, 16),
                        margin = c(0, 0.5), hr_margin = c(0, 0.25), prior_strength = c(0.05, 0.5), n_intervals = c(6, 16)))
B1 <- BOUNDS[[tr$group]]
if (TIER == "t1") { if (SID == "cohortA_minimax_fut") B1$ev <- c(8, 50) } else {                    # RT1:143
  if (tr$group == "A" && PH[[SID]] == "minimax") B1$ev <- c(8, 50) }                                # RT2:162
W0 <- spec$wrapper; we <- environment(W0)                                                           # RT1:144-147
got <- c(w_N = we$w_total_n, w_null = we$w_null, w_alt = we$w_alt)
if (!isTRUE(all.equal(got, WEIGHTS[SID, ], tolerance = 1e-12)) || !isTRUE(all.equal(we$w_EN, 1 - we$w_total_n)))
  stop(SID, ": wrapper weights ", paste(names(got), got, collapse = " "), " != expected ", paste(WEIGHTS[SID, ], collapse = " "))
if (!identical(spec$objective, "weighted_EN_N")) stop("objective ", spec$objective, " != weighted_EN_N")   # RT1:148
cs <- if (is.list(spec$constraints) && !is.null(spec$constraints$stage1)) spec$constraints$stage1 else spec$constraints   # RT1:149-156
BUF <- 0.02 / sqrt(3)
cv <- c(power = as.numeric(cs$power[2]), type1 = as.numeric(cs$type1[2]))
if (!setequal(names(cs), c("power", "type1")) || !identical(cs$power[1], "ge") || !identical(cs$type1[1], "le") ||
    !isTRUE(all.equal(unname(cv), c(0.80 + BUF, 0.10 - BUF), tolerance = 1e-6)))
  stop("stored constraints are not power >= 0.80 + buffer, type1 <= 0.10 - buffer: ", paste(names(cs), sapply(cs, paste, collapse = " "), collapse = "; "))
if (!setequal(names(spec$bounds_stage1), names(B1)) ||                                              # RT1:157-158
    !all(vapply(names(B1), function(p) isTRUE(all.equal(as.numeric(spec$bounds_stage1[[p]]), B1[[p]])), TRUE))) stop("stage-1 bounds differ from F0")
if (!identical(as.integer(we$cached_num_workers), WORKERS)) stop("wrapper cached ", we$cached_num_workers, " workers, expected ", WORKERS)   # RT1:159
FL <- c(low = 3000, med = 5000, high = 10000)                                                       # RT1:160-161
if (!isTRUE(all.equal(unlist(we$fidelity_levels)[c("low", "med", "high")], FL, check.attributes = FALSE))) stop("F0 wrapper fidelity is not 3000/5000/10000")
if (!identical(as.integer(spec$n_init), 60L) || !identical(as.integer(spec$batch_size), 3L)) stop("row n_init / batch size differ from F0 (60 / 3)")   # RT1:162
if (TIER == "t2a" && PROT == "rws" && !identical(as.integer(tr$n_init_expected), 60L)) stop("staged n_init is not 60")   # RT2:182
if (!identical(attr(W0, "trial_type"), tr$trial_type)) stop("wrapper trial_type ", attr(W0, "trial_type"), " != ", tr$trial_type)   # RT1:163
if (TEST_R > 0L) { e <- new.env(parent = we); fl <- we$fidelity_levels; fl[["high"]] <- TEST_R      # RT1:164-165
  assign("fidelity_levels", fl, envir = e); environment(W0) <- e }
# RT1:166-169 (search smoke) and RT1:175-186 (pareto_front_select_best coercion, BO selection only) do not apply to
# certification and are not copied.

# ---- capture hook, installed BEFORE the two engine switches so that it sits innermost (inventory I4) --------------
ns <- asNamespace("evolveTrial")
.cap <- new.env(parent = emptyenv()); .cap$calls <- list()
.cap$orig_rs <- get("run_scenarios", envir = ns)
utils::assignInNamespace("run_scenarios", function(base_args, scens, ...) {
  .cap$calls[[length(.cap$calls) + 1L]] <- list(base_args = base_args, scens = scens, dots = list(...),
                                               looks_stop_at_final_option = getOption("evolveTrial.looks_stop_at_final"),
                                               seed_before = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
  .cap$orig_rs(base_args, scens, ...)
}, ns = "evolveTrial")
HOOK <- get("run_scenarios", envir = ns)

# ---- engine: M4 switch off, then P6 (RT1:189-192) -------------------------------------------------------------------
source(M4S); source(P6R)
M4LOG <- file.path(OUT, "m4_switch_calls.tsv"); M4PROBE <- file.path(OUT, "m4_probe"); P6EV <- file.path(OUT, "p6_events")
verify_m4_install("off", M4LOG, M4PROBE)
p6r_install("p6", "probe", ev_dir = file.path(P6EV, "_init"), expected_lib = LIB)
if (!identical(.vm4$orig_rs, HOOK)) stop("M4 switch does not wrap the capture hook")
PN <- names(spec$bounds_stage1)                                                                     # RT1:193-201
rd_lines <- function(f) if (file.exists(f)) readLines(f, warn = FALSE) else character(0)
probe_lines <- function() unlist(lapply(sort(list.files(M4PROBE, full.names = TRUE)), rd_lines))
EVF <- file.path(OUT, "evals.csv"); .k <- 0L; .SITE <- NULL
site <- function() { cl <- vapply(sys.calls(), function(x) sub("^.*::", "", paste(deparse(x[[1]]), collapse = "")), "")
  if (any(cl %in% c("validate_feasible_designs", "validate_feasible_designs_adaptive"))) "gate"
  else if (any(cl == "run_verification_with_percentiles")) "verify" else if (any(cl == "bo_calibrate")) "bo" else "other" }
engine_fail <- function(msg) { writeLines(c(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), msg), file.path(OUT, "ENGINE-FAIL.txt"))
  cat("[retune] ENGINE FAIL:", msg, "\n"); quit(save = "no", status = 3) }
W <- function(theta, fidelity = "low", seed = NULL, ...) {                                           # RT1:202-246, verbatim
  .k <<- .k + 1L; st <- if (!is.null(.SITE)) .SITE else if (PHASE == "search") site() else PHASE; t0 <- Sys.time()
  n4 <- length(rd_lines(M4LOG)); np <- length(probe_lines()); n6 <- p6r_ncalls()
  p6r_set_ev_dir(file.path(P6EV, sprintf("c%05d", .k)))
  atomic <- is.atomic(theta); coerced <- atomic && COERCE
  if (coerced) { theta <- as.list(theta)
    cat(sprintf("%s call %d site %s seed %s: atomic theta passed as a list (%s)\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), .k, st,
                format(seed %||% NA), paste(names(theta), signif(unlist(theta), 10), sep = "=", collapse = " ")),
        file = file.path(OUT, "coerced_calls.log"), append = TRUE) }
  res <- tryCatch(W0(theta, fidelity = fidelity, seed = seed, ...), error = function(err) err)
  ok <- !inherits(res, "error")
  m4n <- rd_lines(M4LOG); m4n <- if (length(m4n) > n4) m4n[(n4 + 1):length(m4n)] else character(0)
  pr <- probe_lines(); pr <- if (length(pr) > np) pr[(np + 1):length(pr)] else character(0)
  nc <- p6r_ncalls() - n6; ck <- p6r_check_calls(n6, nc)
  ev <- tryCatch(p6r_read_events(.p6r$ev_dir)$W, error = function(e) NULL)
  m4sw <- unique(vapply(strsplit(m4n, "\t"), `[`, "", 3)); m4pass <- unique(vapply(strsplit(m4n, "\t"), `[`, "", 4))
  prv <- unique(vapply(strsplit(pr, "\t"), `[`, "", 3))
  p6F <- if (is.null(ev)) character(0) else unique(ev$frame_F); p6o <- if (is.null(ev)) character(0) else unique(ev$frame_opt)
  p6b <- if (is.null(ev)) character(0) else unique(ev$build)
  echo <- sprintf("p6[calls=%d F=%s opt=%s workers=%d frameF=%s frameOpt=%s build=%s] dropout[calls=%d switch=%s passed=%s workers=%d seen=%s]",
                  nc, paste(ck$F_values, collapse = "/"), paste(ck$option_values, collapse = "/"), if (is.null(ev)) 0L else length(unique(ev$pid)),
                  paste(p6F, collapse = "/"), paste(p6o, collapse = "/"), substr(paste(p6b, collapse = "/"), 1, 7),
                  length(m4n), paste(m4sw, collapse = "/"), paste(m4pass, collapse = "/"), length(pr), paste(prv, collapse = "/"))
  eng_ok <- !ok || (nc >= 1L && ck$ok && length(m4n) == nc && identical(m4sw, "off") && identical(m4pass, "FALSE") &&
                    length(pr) >= 1L && identical(prv, "FALSE") && identical(p6F, "6") && identical(p6o, "TRUE") &&
                    identical(p6b, P6R_CONST$et_sha))
  g <- function(n) { v <- if (ok) tryCatch(res[[n]], error = function(e) NULL) else NULL
    if (!is.null(v) && length(v) == 1L) as.numeric(v) else NA_real_ }
  th <- unlist(theta)[PN]; R <- tryCatch(unname(environment(W0)$fidelity_levels[[as.character(fidelity)]]), error = function(e) NA_real_)
  row <- data.frame(call = .k, phase = PHASE, site = st, fidelity = as.character(fidelity), R = R, seed = seed %||% NA,
                    theta_atomic = atomic, theta_coerced = coerced, as.list(th), status = if (ok) "ok" else paste("error:", conditionMessage(res)),
                    power = g("power"), type1 = g("type1"), EN_null = g("EN_null"), EN_alt = g("EN_alt"), EN = g("EN"),
                    weighted_EN_N = g("weighted_EN_N"), metric_total_n = g("metric_total_n"),
                    w_N = we$w_total_n, w_EN = we$w_EN, w_null = we$w_null, w_alt = we$w_alt,
                    engine_ok = eng_ok, engine_echo = echo, wall_s = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  utils::write.table(row, EVF, sep = ",", row.names = FALSE, col.names = !file.exists(EVF), append = file.exists(EVF), qmethod = "double")
  if (!eng_ok) engine_fail(sprintf("call %d: %s", .k, echo))
  if (!ok) stop(res)
  res
}
attributes(W) <- attributes(W0)[setdiff(names(attributes(W0)), "srcref")]                          # RT1:247
counts <- function(m, R) { f <- function(p) { x <- p * R; if (is.finite(x) && abs(x - round(x)) < 1e-6) round(x) else NA_real_ }   # RT1:261-262
  c(power_count = f(m[["power"]]), type1_count = f(m[["type1"]])) }
theta_of <- function(json) as.list(jsonlite::fromJSON(json))                                        # RT1:263

# ---- certify (RT1:371-405, certify branch; RT2:411-412 for the design id) -------------------------------------------
th <- theta_of(cr$theta_json)                                                                       # RT1:386
if (!setequal(names(th), PN)) stop("theta names differ from the stage-1 parameters")               # RT1:388-389
for (p in PN) if (th[[p]] < spec$bounds_stage1[[p]][1] - 1e-9 || th[[p]] > spec$bounds_stage1[[p]][2] + 1e-9) stop("theta ", p, " outside stage-1 bounds")
if (TIER == "t1") {                                                                                 # RT1:391-392
  if (!identical(sub("p$", "", tr$design), cr$design)) stop("task ", TASK, " is not a run of design ", cr$design)
} else if (!identical(tr$base_design, cr$design)) stop("task ", TASK, " is not a run of design ", cr$design)   # RT2:411-412
setwd(F0)                                                                                           # RT1:393
pt0 <- proc.time(); t0 <- Sys.time()                                                                # RT1:395-397
m <- W(th, fidelity = "high", seed = CSEED)
wall <- as.numeric(difftime(Sys.time(), t0, units = "secs")); cpu <- proc.time() - pt0
setwd(old)
RW <- if (TEST_R > 0L) TEST_R else 10000L                                                           # RT1:398-399
cc <- counts(m, RW)

# ---- capture checks (plan task 2a.2) ---------------------------------------------------------------------------------
if (length(.cap$calls) != 1L) stop("expected exactly one run_scenarios call, saw ", length(.cap$calls))
cap1 <- .cap$calls[[1]]; ba <- cap1$base_args
K <- as.integer(th$n_intervals)
chk <- c(cutpoints_eq_n_intervals = length(ba$interval_cutpoints_sim) - 1L == K,
         random_dropout_FALSE = identical(ba$random_dropout, FALSE),
         min_follow_up_at_final_6 = identical(ba$min_follow_up_at_final, 6),
         looks_stop_at_final_option = isTRUE(getOption("evolveTrial.looks_stop_at_final")) && isTRUE(cap1$looks_stop_at_final_option),
         num_workers_19 = identical(as.numeric(ba$num_workers), 19),
         num_simulations_R = identical(as.numeric(ba$num_simulations), as.numeric(RW)),
         n_rep_R = identical(as.numeric(attr(m, "n_rep") %||% NA), as.numeric(RW)),
         counts_integral = !anyNA(cc))
evrow <- utils::read.csv(EVF, stringsAsFactors = FALSE)
saveRDS(list(design = DESIGN, task = TASK, scenario_id = SID, trial_type = tr$trial_type, protocol = PROT, csv = tr$csv,
             theta = th, cert_seed = CSEED, base_args = ba, scens = cap1$scens, dots = cap1$dots,
             seed_before_run_scenarios = cap1$seed_before, wrapper_res = m, counts = cc, R = RW, tier = TIER,
             runner = RUNNER, runner_md5 = RUNNER_MD5, checks = chk, evals_row = evrow,
             wall_s = wall, cpu = cpu, env = ENV, sessionInfo = sessionInfo()), OUTF)
cat(sprintf("[capture] %s %s %s %s task %s cert_seed %d R %d K %d cutpoints %d power %.4f type1 %.4f E0N %.4f E1N %.4f wall %.1f s\n",
            DESIGN, SID, tr$trial_type, PROT, TASK, CSEED, RW, K, length(ba$interval_cutpoints_sim) - 1L,
            m[["power"]], m[["type1"]], m[["EN_null"]], m[["EN_alt"]], wall))
cat("[capture] engine echo:", evrow$engine_echo, "engine_ok", evrow$engine_ok, "\n")
cat("[capture] checks:", paste(names(chk), chk, sep = "=", collapse = " "), "\n")
if (!all(chk)) { cat("[capture] CHECK FAILED\n"); quit(save = "no", status = 4) }
