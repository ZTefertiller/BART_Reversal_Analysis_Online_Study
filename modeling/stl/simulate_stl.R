# modeling/stl/simulate_stl.R
# Generative simulator for the STL BART model. Mirrors the decision process in
# modeling/stl/stl_vec.stan exactly:
#   - each trial has a target omega; on trial 1, omega = nmax * omegaone
#   - at pump opportunity n the agent pumps with prob logistic(-beta * (n - omega))
#   - the balloon bursts if the agent pumps to its (fixed) breakpoint
#   - between trials omega is scaled by the previous trial's outcome:
#       burst : omega <- omega * (1 - vloss * (1 - npumps_prev / nmax))
#       cash  : omega <- omega * (1 + vwin  * (    npumps_prev / nmax))
#     (nmax is constant within a condition, so the nmax[k]/nmax[last] ratio in
#      the Stan code reduces to 1)
#
# Environment: a vector of true breakpoints (one per trial), shared by all
# subjects (as in the real fixed-order task), plus the balloon max (color_max).

simulate_stl_subject <- function(vwin, vloss, beta, omegaone,
                                 breakpoints, maxpump = 128L) {
  ntrial <- length(breakpoints)
  npumps <- integer(ntrial)
  popped <- integer(ntrial)

  omega <- maxpump * omegaone

  for (k in seq_len(ntrial)) {
    if (k > 1L) {
      prev_frac <- npumps[k - 1L] / maxpump
      omega <- if (popped[k - 1L] == 1L) {
        omega * (1 - vloss * (1 - prev_frac))
      } else {
        omega * (1 + vwin * prev_frac)
      }
    }

    Bk <- breakpoints[k]
    n <- 1L
    this_pumps <- 0L
    burst <- 0L
    repeat {
      p_pump <- plogis(-beta * (n - omega))
      if (runif(1) < p_pump) {            # decide to pump
        if (n >= Bk) {                    # this pump hits the breakpoint -> burst
          this_pumps <- n; burst <- 1L; break
        }
        n <- n + 1L
        if (n > maxpump) {                # safety cap (breakpoints < maxpump)
          this_pumps <- maxpump; burst <- 0L; break
        }
      } else {                            # decide to stop / cash in
        this_pumps <- n - 1L; burst <- 0L; break
      }
    }

    npumps[k] <- this_pumps
    popped[k] <- burst
  }

  list(npumps = npumps, popped = popped)
}

# params_df: data.frame with columns sub_id, vwin, vloss, beta, omegaone
# Returns a trial-level data.frame ready for create_stan_params().
simulate_stl_dataset <- function(params_df, breakpoints, color_max = 128L,
                                 balloon_color = "b", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ntrial <- length(breakpoints)

  do.call(rbind, lapply(seq_len(nrow(params_df)), function(i) {
    r <- params_df[i, ]
    sim <- simulate_stl_subject(r$vwin, r$vloss, r$beta, r$omegaone,
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
