# modeling/fourpar/simulate_fourpar.R
# Generative simulator for the 4PAR BART model. Mirrors the decision process in
# modeling/fourpar/4par.stan exactly, including its numerical guards:
#   - belief about NOT bursting, from cumulative pumps/successes:
#       n_pump == 0 : belief = phi
#       otherwise   : belief = (phi + eta * n_succ) / (1 + eta * n_pump)
#     clamped to [0.001, 0.999]; p_burst = 1 - belief, clamped the same way
#   - target omega = -gam / log1m(p_burst), with log1m(p_burst) capped at
#     -0.001 and omega clamped to [-1000, 1000]
#   - at pump opportunity l the agent pumps with prob logistic(tau * (omega - l)),
#     the linear predictor clamped to [-100, 100]
#   - the balloon bursts if the agent pumps to its (fixed) breakpoint
#
# Environment: a vector of true breakpoints (one per trial), shared by all
# subjects (as in the real fixed-order task), plus the balloon max (color_max).

simulate_fourpar_subject <- function(phi, eta, gam, tau,
                                     breakpoints, maxpump = 128L) {
  ntrial <- length(breakpoints)
  npumps <- integer(ntrial)
  popped <- integer(ntrial)

  n_succ <- 0
  n_pump <- 0

  for (k in seq_len(ntrial)) {
    belief <- if (n_pump == 0) phi else (phi + eta * n_succ) / (1 + eta * n_pump)
    belief  <- min(0.999, max(0.001, belief))
    p_burst <- min(0.999, max(0.001, 1 - belief))

    denom <- min(-0.001, log1p(-p_burst))
    omega <- min(1000, max(-1000, -gam / denom))

    Bk <- breakpoints[k]
    l  <- 1L
    this_pumps <- 0L
    burst <- 0L
    repeat {
      linpred <- min(100, max(-100, tau * (omega - l)))
      if (runif(1) < plogis(linpred)) {   # decide to pump
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
  }

  list(npumps = npumps, popped = popped)
}

# params_df: data.frame with columns sub_id, phi, eta, gam, tau
# Returns a trial-level data.frame ready for create_stan_params().
simulate_fourpar_dataset <- function(params_df, breakpoints, color_max = 128L,
                                     balloon_color = "b", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ntrial <- length(breakpoints)

  do.call(rbind, lapply(seq_len(nrow(params_df)), function(i) {
    r <- params_df[i, ]
    sim <- simulate_fourpar_subject(r$phi, r$eta, r$gam, r$tau,
                                    breakpoints, maxpump = color_max)
    data.frame(
      sub_id         = r$sub_id,
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
