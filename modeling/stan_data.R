# ===========================================================================
# modeling/stan_data.R
# ---------------------------------------------------------------------------
# create_stan_params(): trial-level data frame -> the Stan data list.
#
# Returns the SUPERSET of the fields the three models declare, so the same list
# can be handed to stl_vec.stan, 4par.stan and ewmv_vec.stan (CmdStan ignores
# data entries a model does not declare):
#
#   stl_vec.stan  nsub ntrial outcome npumps opportunity nmax maxpump d
#   4par.stan     nsub ntrial Tsubj maxpump npumps outcome   (builds its own d)
#   ewmv_vec.stan nsub subj_idx ntrial maxpump npumps outcome d
#
# Sourced by modeling/fit_models.R and modeling/run_recovery.R so the real fits
# and the recovery refits are built by one definition.
#
# Input needs: sub_id, trial_number, popped, inflations, color_max, balloon_color.
# Subjects must all have the same number of trials.
# ===========================================================================
suppressPackageStartupMessages({
  library(dplyr)
})

create_stan_params <- function(df) {
  need <- c("sub_id", "trial_number", "popped", "inflations", "color_max", "balloon_color")
  if (!all(need %in% names(df)))
    stop("Missing: ", paste(setdiff(need, names(df)), collapse = ", "))
  df <- df %>% dplyr::arrange(sub_id, trial_number)
  cnt <- df %>% dplyr::count(sub_id, name = "n")
  if (length(unique(cnt$n)) != 1L) stop("Unequal trials per subject.")
  nsub <- nrow(cnt); ntrial <- cnt$n[1]
  outcome <- matrix(df$popped,     nrow = nsub, ncol = ntrial, byrow = TRUE)
  npumps  <- matrix(df$inflations, nrow = nsub, ncol = ntrial, byrow = TRUE)
  nmax    <- matrix(df$color_max,  nrow = nsub, ncol = ntrial, byrow = TRUE)
  opportunity <- npumps + (1 - outcome); maxpump <- max(nmax)
  df$balloon_color_num <- match(df$balloon_color, c("b", "o", "y", "p"))
  balloon_color <- matrix(df$balloon_color_num, nrow = nsub, ncol = ntrial, byrow = TRUE)
  d <- array(NA_integer_, dim = c(nsub, ntrial, maxpump))
  for (i in 1:nsub) for (j in 1:ntrial) {
    pumps <- npumps[i, j]; out <- outcome[i, j]
    if (pumps > 0) d[i, j, 1:pumps] <- 1L
    if (pumps < maxpump && out == 0) d[i, j, pumps + 1] <- 0L
  }
  d[is.na(d)] <- 75L
  list(nsub = nsub, ntrial = ntrial, Tsubj = rep(ntrial, nsub), outcome = outcome,
       npumps = npumps, nmax = nmax, opportunity = opportunity, maxpump = maxpump,
       d = d, balloon_color = balloon_color,
       sub_ids = seq_len(nsub), subj_idx = seq_len(nsub))
}

# The three task phases, as filters on the published trial file.
CONDITION_FILTERS <- list(
  blue   = quote(balloon_color == "b" & trial_number < 91),   # pre-reversal
  orange = quote(balloon_color == "o" & trial_number > 90),   # post-reversal
  pink   = quote(balloon_color == "p")                        # control
)

# Trials of one condition, from the published trial file.
condition_trials <- function(df, cond) {
  f <- CONDITION_FILTERS[[cond]]
  if (is.null(f)) stop("Unknown condition: ", cond,
                       " (expected one of ", paste(names(CONDITION_FILTERS), collapse = ", "), ")")
  df %>% dplyr::filter(!!f)
}

# The condition's fixed breakpoint sequence (one per trial, shared by every
# subject, as in the real fixed-order task).
condition_breakpoints <- function(cdf) {
  cdf %>%
    dplyr::group_by(trial_number) %>%
    dplyr::summarise(brk = dplyr::first(optimal_inflations), .groups = "drop") %>%
    dplyr::arrange(trial_number) %>%
    dplyr::pull(brk)
}
