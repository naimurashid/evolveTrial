# P3-1: an instrumented copy of the R body (simulate_chunk_r) that counts, per trial, which P3-1 branches it took.
# The copy is built by deparsing the frame's closure and inserting counter calls at fixed statements; it draws exactly
# what the R body draws (gate B checks its result and RNG state against the real R body on every chunk). Counters for
# the interim threshold resolution come from a trace on calculate_current_probs_hc (namespace), active only while the
# instrumented copy runs.
.p31B <- new.env()
p31_counter_reset <- function() {
  .p31B$trials <- 0L
  .p31B$cur <- character(0)
  .p31B$trial_hits <- list()   # name -> number of trials with at least one hit
  .p31B$events <- list()       # name -> total number of hits
}
p31_counter_reset()
.p31_flush <- function() {
  if (.p31B$trials > 0L) for (k in unique(.p31B$cur)) .p31B$trial_hits[[k]] <- (.p31B$trial_hits[[k]] %||% 0L) + 1L
  .p31B$cur <- character(0)
}
.p31B$new_trial <- function() { .p31_flush(); .p31B$trials <- .p31B$trials + 1L }
.p31B$hit <- function(k) {
  .p31B$cur <- c(.p31B$cur, k)
  .p31B$events[[k]] <- (.p31B$events[[k]] %||% 0L) + 1L
}
.p31B$drop <- function(rc, te) if (is.finite(rc) && rc < te) .p31B$hit("dropout_censors_event")
.p31B$interim <- function(args) {
  .p31B$hit("hc_interim_draws")
  if (is.null(args$futility_threshold_hc_prob)) .p31B$hit("hc_interim_no_fut_rule")
  if (is.null(args$median_pfs_success_threshold_arms)) .p31B$hit("hc_interim_success_from_null_median")
  if (is.null(args$futility_median_arms)) .p31B$hit("hc_interim_futility_from_success")
}
p31_counts <- function() {
  .p31_flush()
  list(trials = .p31B$trials, trial_hits = unlist(.p31B$trial_hits), events = unlist(.p31B$events))
}

p31_instrument <- function(env) {
  x <- deparse(env$simulate_chunk_r, width.cutoff = 500L)
  tr <- trimws(x)
  out <- character(0)
  n_occ <- list()
  occ <- function(k) { n_occ[[k]] <<- (n_occ[[k]] %||% 0L) + 1L; n_occ[[k]] }
  i <- 1L
  while (i <= length(x)) {
    t <- tr[i]; line <- x[i]
    if (t == "for (sim_idx in sim_indices) {") {
      out <- c(out, line, ".p31B$new_trial()")
    } else if (t == "current_time <- max_trial_time") {
      out <- c(out, '{ .p31B$hit("cutoff_amd8"); current_time <- max_trial_time }')
    } else if (t == "if (looks_stop_at_final && current_time >= max_trial_time) {") {
      out <- c(out, line, '.p31B$hit("cutoff_p6")')
    } else if (t == "interim_ok_here <- FALSE") {
      k <- if (occ("iok") == 1L) "last_slot_no_look" else "vsref_phantom_no_look"
      out <- c(out, sprintf('{ .p31B$hit("%s"); interim_ok_here <- FALSE }', k))
    } else if (t == "if (current_time >= next_calendar_look) {") {
      out <- c(out, line, '.p31B$hit("accrual_look")',
               paste0('if (isTRUE(compare_arms_option) && !any(which(arm_names != reference_arm_name) %in% ',
                      'is_eligible(state))) .p31B$hit("vsref_phantom_look")'))
    } else if (t == "if (length(exp_active) == 0)" && trimws(x[i + 1L]) == "break") {
      out <- c(out, line, '{ .p31B$hit("vsref_exit_break"); break }')
      i <- i + 1L
    } else if (t == "t_random_censor <- Inf") {
      out <- c(out, line, ".p31B$drop(t_random_censor, t_event_true)")
    } else if (t == "if (current_time < t_final_planned) {") {
      out <- c(out, line, '.p31B$hit("j0_look")')
    } else if (startsWith(t, "state <- interim_check(state, post_enroll_time,")) {
      k <- if (occ("post") == 1L) "post_look_amd8" else "post_look_p6"
      out <- c(out, sprintf('.p31B$hit("%s")', k), line)
    } else if (t == "last_enroll_time + min_follow_up_at_final") {
      k <- c("inf_t_final_planned", "inf_t_star", "inf_final_time")[occ("inf")]
      out <- c(out, sprintf('{ .p31B$hit("%s"); last_enroll_time + min_follow_up_at_final }', k))
    } else if (t == "Inf") {
      out <- c(out, '{ .p31B$hit("final_success_fallback_Inf"); Inf }')
    } else if (t == "-Inf") {
      out <- c(out, '{ .p31B$hit("final_futility_fallback_negInf"); -Inf }')
    } else if (startsWith(t, "new_cuts <- rebalance_cutpoints_from_events(")) {
      out <- c(out, '.p31B$hit("rebalanced")', line)
    } else {
      out <- c(out, line)
    }
    i <- i + 1L
  }
  need <- c(occ = 2L)
  stopifnot(identical(n_occ$iok, 2L), identical(n_occ$post, 2L), identical(n_occ$inf, 3L))
  e2 <- new.env(parent = env)
  e2$.p31B <- .p31B
  f <- eval(parse(text = out), envir = e2)
  f
}

# Run fn with the interim trace active.
p31_with_trace <- function(expr) {
  ns <- asNamespace("evolveTrial")
  suppressMessages(trace("calculate_current_probs_hc", where = ns, print = FALSE,
                         tracer = quote(get(".p31B", envir = globalenv())$interim(args))))
  on.exit(suppressMessages(untrace("calculate_current_probs_hc", where = ns)), add = TRUE)
  force(expr)
}
