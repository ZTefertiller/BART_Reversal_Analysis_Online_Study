# modeling/ewmv/run_ewmv_recovery.R
# EWMV parameter recovery for one condition.
#   1. TRUE params are drawn from that condition's FITTED population model
#      (group means mu_pr + sd sigma from ewmv_<cond>_allrecruitments.rds, transformed
#      to the native scale). Drawing from the estimated population gives realistic
#      between-subject spread, so recovery is NOT deflated by the shrinkage you get
#      if you (mis)use the per-subject posterior means as "truth".
#   2. Simulate BART data from the EWMV generative model on that condition's
#      fixed breakpoints (simulate_ewmv.R).
#   3. Refit ewmv_vec.stan; compare TRUE vs RECOVERED per-subject means.
#
# Usage:  EMBED_RECOVERY_COND=blue|orange|yellow|pink Rscript run_ewmv_recovery.R
# Small core budget by default (4 chains x 1 thread).

suppressPackageStartupMessages({
  library(dplyr); library(cmdstanr); library(posterior); library(here)
})
source(here("modeling", "ewmv", "simulate_ewmv.R"))

param_names <- c("phi", "eta", "rho", "tau", "lambda")

create_stan_params <- function(df) {
  df <- df %>% arrange(participant_id, trial_number)
  cnt <- df %>% count(participant_id, name = "n")
  if (length(unique(cnt$n)) != 1L) stop("Unequal trials per subject.")
  nsub <- nrow(cnt); ntrial <- cnt$n[1]
  outcome <- matrix(df$popped,     nrow = nsub, ncol = ntrial, byrow = TRUE)
  npumps  <- matrix(df$inflations, nrow = nsub, ncol = ntrial, byrow = TRUE)
  nmax    <- matrix(df$color_max,  nrow = nsub, ncol = ntrial, byrow = TRUE)
  maxpump <- max(nmax)
  d <- array(NA_integer_, dim = c(nsub, ntrial, maxpump))
  for (i in 1:nsub) for (j in 1:ntrial) {
    pumps <- npumps[i, j]; out <- outcome[i, j]
    if (pumps > 0) d[i, j, 1:pumps] <- 1L
    if (pumps < maxpump && out == 0) d[i, j, pumps + 1] <- 0L
  }
  d[is.na(d)] <- 75L
  list(nsub = nsub, ntrial = ntrial, maxpump = maxpump,
       subj_idx = seq_len(nsub), npumps = npumps, outcome = outcome, d = d)
}

cond_spec <- list(
  blue   = list(rds = "ewmv_blue_allrecruitments.rds",   filt = quote(balloon_color == "b" & trial_number < 91)),
  orange = list(rds = "ewmv_orange_allrecruitments.rds", filt = quote(balloon_color == "o" & trial_number > 90)),
  yellow = list(rds = "ewmv_yellow_allrecruitments.rds", filt = quote(balloon_color == "y")),
  pink   = list(rds = "ewmv_pink_allrecruitments.rds",   filt = quote(balloon_color == "p"))
)

main <- function() {
  cond    <- Sys.getenv("EMBED_RECOVERY_COND", "blue")
  warmup  <- as.integer(Sys.getenv("EMBED_WARMUP", "1000"))
  sample  <- as.integer(Sys.getenv("EMBED_SAMPLE", "1000"))
  chains  <- as.integer(Sys.getenv("EMBED_CHAINS", "4"))
  threads <- as.integer(Sys.getenv("EMBED_THREADS", "1"))
  seed    <- as.integer(Sys.getenv("EMBED_SEED", "123"))
  spec    <- cond_spec[[cond]]
  if (is.null(spec)) stop("Unknown condition: ", cond)

  out_dir <- here("mcmc", "ewmv_recovery")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fit_path <- here("mcmc", "ewmv", spec$rds)
  if (!file.exists(fit_path)) stop("Missing fit: ", fit_path)
  fit <- readRDS(fit_path)

  # fitted population: group means + sds on the unconstrained scale
  mu  <- vapply(1:5, function(k) mean(fit$draws(sprintf("mu_pr[%d]", k), format = "df")[[sprintf("mu_pr[%d]", k)]]), numeric(1))
  sig <- vapply(1:5, function(k) mean(fit$draws(sprintf("sigma[%d]", k), format = "df")[[sprintf("sigma[%d]", k)]]), numeric(1))

  set.seed(seed)
  n <- 118
  zr <- function(k) mu[k] + sig[k] * rnorm(n)
  true_params <- tibble(
    participant_id = sprintf("s%03d", seq_len(n)),
    phi    = pnorm(zr(1)),
    eta    = pnorm(zr(2)),
    rho    = 0.5 - pnorm(zr(3)),
    tau    = exp(zr(4)),
    lambda = exp(zr(5))
  )

  data <- read.csv(here("data", "processed", "filtered_data", "full_all_recruitments.csv"))
  cdf  <- data %>% filter(!!spec$filt)
  bp   <- cdf %>% group_by(trial_number) %>%
    summarise(brk = first(optimal_inflations), .groups = "drop") %>%
    arrange(trial_number) %>% pull(brk)
  cmax <- max(cdf$color_max)
  cat(sprintf("[%s] ntrial=%d color_max=%d  mean true tau=%.2f\n",
              cond, length(bp), cmax, mean(true_params$tau)))

  sim <- simulate_ewmv_dataset(true_params, bp, color_max = cmax, seed = seed)
  stan_data <- create_stan_params(sim)

  mod <- cmdstan_model(here("modeling", "ewmv", "ewmv_vec.stan"),
                       cpp_options = list(stan_threads = TRUE))
  # cap treedepth (recovery only) so a single stuck chain can't crawl for hours
  treedepth <- as.integer(Sys.getenv("EMBED_TREEDEPTH", "10"))
  rfit <- mod$sample(
    data = stan_data, iter_warmup = warmup, iter_sampling = sample,
    chains = chains, parallel_chains = chains, threads_per_chain = threads,
    seed = 10191998, adapt_delta = 0.9, max_treedepth = treedepth, refresh = 200
  )

  recovered <- as.data.frame(sapply(param_names, function(p) {
    dr <- rfit$draws(variables = p, format = "df")
    cols <- grep(paste0("^", p, "\\["), names(dr), value = TRUE)
    idx <- as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", cols)); cols <- cols[order(idx)]
    vapply(cols, function(cc) mean(dr[[cc]]), numeric(1))
  }))
  names(recovered) <- paste0("rec_", param_names)

  cmp <- bind_cols(
    true_params %>% select(all_of(param_names)) %>%
      rename_with(~ paste0("true_", .x)),
    recovered
  ) %>% mutate(condition = cond)

  out_csv <- file.path(out_dir, paste0("recovery_", cond, ".csv"))
  write.csv(cmp, out_csv, row.names = FALSE)
  cat("Recovery r:",
      paste(sprintf("%s=%.2f", param_names,
                    vapply(param_names, function(p)
                      cor(cmp[[paste0("true_", p)]], cmp[[paste0("rec_", p)]]),
                      numeric(1))), collapse = "  "), "\n")
  cat("Wrote", out_csv, "\n")
}

main()
