# modeling/ewmv/simulate_ewmv.R
# Generative simulator for the EWMV BART model. Mirrors the decision process in
# modeling/ewmv/ewmv_vec.stan exactly:
#   - within a trial, p_burst (belief) is fixed; at pump opportunity l (0-indexed
#     loss u_loss = l-1) the agent pumps with prob logistic(tau * u_pump),
#       u_pump = (1-pb) - lambda*pb*(l-1) + rho*pb*(1-pb)*(1 + lambda*(l-1))^2
#   - the balloon bursts if the agent pumps to its (fixed) breakpoint
#   - between trials, p_burst updates from cumulative pumps/successes:
#       p_burst = phi + (1 - exp(-n_pump*eta)) * ((n_pump - n_succ)/n_pump - phi)
#
# Environment: a vector of true breakpoints (one per trial), shared by all
# subjects (as in the real fixed-order task), plus the balloon max (color_max).

simulate_ewmv_subject <- function(phi, eta, rho, tau, lambda,
                                  breakpoints, maxpump = 128L) {
  ntrial <- length(breakpoints)
  npumps <- integer(ntrial)
  popped <- integer(ntrial)

  p_burst <- phi
  n_succ <- 0
  n_pump <- 0

  for (k in seq_len(ntrial)) {
    Bk <- breakpoints[k]
    l  <- 1L
    this_pumps <- 0L
    burst <- 0L
    repeat {
      u_loss <- l - 1
      u_pump <- (1 - p_burst) -
                lambda * p_burst * u_loss +
                rho * p_burst * (1 - p_burst) * (1 + lambda * u_loss)^2
      p_pump <- plogis(tau * u_pump)

      if (runif(1) < p_pump) {            # decide to pump
        if (l >= Bk) {                    # this pump hits the breakpoint -> burst
          this_pumps <- l; burst <- 1L; break
        }
        l <- l + 1L
        if (l > maxpump) {                # safety cap (breakpoints < maxpump)
          this_pumps <- maxpump; burst <- 0L; break
        }
      } else {                            # decide to stop / cash in
        this_pumps <- l - 1L; burst <- 0L; break
      }
    }

    npumps[k] <- this_pumps
    popped[k] <- burst
    n_succ <- n_succ + (this_pumps - burst)
    n_pump <- n_pump + this_pumps
    if (n_pump > 0) {
      p_burst <- phi + (1 - exp(-n_pump * eta)) *
                 ((n_pump - n_succ) / n_pump - phi)
    }
  }

  list(npumps = npumps, popped = popped)
}

# params_df: data.frame with columns participant_id, phi, eta, rho, tau, lambda
# Returns a trial-level data.frame ready for create_stan_params().
simulate_ewmv_dataset <- function(params_df, breakpoints, color_max = 128L,
                                  balloon_color = "b", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ntrial <- length(breakpoints)

  do.call(rbind, lapply(seq_len(nrow(params_df)), function(i) {
    r <- params_df[i, ]
    sim <- simulate_ewmv_subject(r$phi, r$eta, r$rho, r$tau, r$lambda,
                                 breakpoints, maxpump = color_max)
    data.frame(
      participant_id = r$participant_id,
      trial_number   = seq_len(ntrial),
      balloon_color  = balloon_color,
      inflations     = sim$npumps,
      popped         = sim$popped,
      color_max      = color_max,
      optimal_inflations = breakpoints,
      stringsAsFactors = FALSE
    )
  }))
}
