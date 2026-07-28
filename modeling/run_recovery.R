#!/usr/bin/env Rscript
# ===========================================================================
# modeling/run_recovery.R
# ---------------------------------------------------------------------------
# Parameter recovery for any of the three BART models (STL, 4PAR, EWMV) in any
# of the three task phases (blue = pre-reversal, orange = post-reversal,
# pink = control).
#
#   1. TRUE parameters are drawn from that model x condition's POPULATION
#      distribution. Drawing from the population (rather than reusing the
#      per-subject posterior means as "truth") gives realistic between-subject
#      spread, so recovery is not deflated by shrinkage.
#   2. Data are simulated from the model's generative process on that
#      condition's fixed breakpoints, read from the published trial file.
#   3. The model is refit to the simulated data; TRUE vs RECOVERED per-subject
#      posterior means are written to
#      mcmc/recovery/recovery_<model>_<condition>.csv
#
# Usage
#   Rscript modeling/run_recovery.R                    # all 9 model x condition cells
#   Rscript modeling/run_recovery.R ewmv               # one model, all conditions
#   Rscript modeling/run_recovery.R ewmv blue          # one cell
#   Rscript modeling/run_recovery.R all orange         # all models, one condition
#
# Environment overrides
#   RECOVERY_TRUTH    auto (default) | fit | published
#                     Where the population parameters come from.
#                       fit       - the fitted group mean/sd from the real fit
#                                   objects (mcmc/model_fits_list.rds or
#                                   mcmc/fits/); the preferred source.
#                       published - moment-matched to the per-participant
#                                   posterior means in data_published/. Needs no
#                                   fit objects, so it runs on a clean clone.
#                                   Those means are shrunk toward the group, so
#                                   the implied population SD is a little narrow
#                                   and recovery r is, if anything, conservative.
#                       auto      - fit if the objects are present, else published.
#   RECOVERY_N        subjects to simulate (default 118, the study N)
#   RECOVERY_WARMUP / RECOVERY_SAMPLE / RECOVERY_CHAINS / RECOVERY_THREADS
#   RECOVERY_SEED     default 123
#   RECOVERY_TREEDEPTH  default 10 (capped for recovery so one stuck chain
#                       cannot crawl for hours)
# ===========================================================================
suppressPackageStartupMessages({
  library(dplyr); library(cmdstanr); library(posterior); library(here)
})

source(here("modeling", "stan_data.R"))
source(here("modeling", "stl", "simulate_stl.R"))
source(here("modeling", "fourpar", "simulate_fourpar.R"))
source(here("modeling", "ewmv", "simulate_ewmv.R"))

log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

# --- per-model spec ---------------------------------------------------------
# link: native <-> unconstrained, matching each .stan file's transform.
#   ewmv     phi, eta ~ Phi(z); rho = 0.5 - Phi(z); tau, lambda = exp(z)
#   fourpar  phi = 0.01 + 0.98*Phi(z); eta, gam, tau = 0.01 + exp(z)
#   stl      subject parameters are sampled on the native scale, truncated to
#            their declared bounds, so the link is the identity.
.link_probit  <- list(to_nat = function(z) pnorm(z),        to_unc = function(x) qnorm(x))
.link_rho     <- list(to_nat = function(z) 0.5 - pnorm(z),  to_unc = function(x) qnorm(0.5 - x))
.link_log     <- list(to_nat = function(z) exp(z),          to_unc = function(x) log(x))
.link_4par_phi<- list(to_nat = function(z) 0.01 + 0.98 * pnorm(z),
                      to_unc = function(x) qnorm((x - 0.01) / 0.98))
.link_4par_pos<- list(to_nat = function(z) 0.01 + exp(z),
                      to_unc = function(x) log(pmax(x - 0.01, .Machine$double.eps)))
.link_id      <- list(to_nat = function(z) z,               to_unc = function(x) x)

MODELS <- list(
  stl = list(
    stan   = c("modeling", "stl", "stl_vec.stan"),
    params = c("vwin", "vloss", "beta", "omegaone"),
    links  = list(vwin = .link_id, vloss = .link_id, beta = .link_id, omegaone = .link_id),
    # native-scale truncated normal, bounds as declared in stl_vec.stan
    bounds = list(vwin = c(0, 1), vloss = c(0, 1), beta = c(0, 3), omegaone = c(0, 1)),
    mu_vars = c("mu_vwin", "mu_vloss", "mu_beta", "mu_omegaone"),
    simulate = simulate_stl_dataset
  ),
  fourpar = list(
    stan   = c("modeling", "fourpar", "4par.stan"),
    params = c("phi", "eta", "gam", "tau"),
    links  = list(phi = .link_4par_phi, eta = .link_4par_pos,
                  gam = .link_4par_pos, tau = .link_4par_pos),
    bounds = NULL,
    mu_vars = sprintf("mu_pr[%d]", 1:4),
    simulate = simulate_fourpar_dataset
  ),
  ewmv = list(
    stan   = c("modeling", "ewmv", "ewmv_vec.stan"),
    params = c("phi", "eta", "rho", "tau", "lambda"),
    links  = list(phi = .link_probit, eta = .link_probit, rho = .link_rho,
                  tau = .link_log, lambda = .link_log),
    bounds = NULL,
    mu_vars = sprintf("mu_pr[%d]", 1:5),
    simulate = simulate_ewmv_dataset
  )
)

# --- population parameters --------------------------------------------------
.draw_mean <- function(fit, v) mean(as_draws_matrix(fit$draws(variables = v))[, 1])

# Locate a saved fit for this cell: the assembled list first, then mcmc/fits/.
.find_fit <- function(model, cond) {
  bundle <- here("mcmc", "model_fits_list.rds")
  if (file.exists(bundle)) {
    fl <- readRDS(bundle)
    nm <- paste0(model, "_", cond)
    if (!is.null(fl[[nm]])) return(fl[[nm]])
  }
  hits <- list.files(here("mcmc", "fits", model), full.names = TRUE,
                     pattern = sprintf("^%s_%s(_n\\d+)?\\.rds$", model, cond))
  if (length(hits)) return(readRDS(hits[which.max(file.mtime(hits))]))
  NULL
}

# mu / sigma on the model's own scale, from the real fit.
.pop_from_fit <- function(spec, fit) {
  list(mu    = vapply(spec$mu_vars, function(v) .draw_mean(fit, v), numeric(1)),
       sigma = vapply(seq_along(spec$params),
                      function(k) .draw_mean(fit, sprintf("sigma[%d]", k)), numeric(1)))
}

# mu / sigma moment-matched to the published per-participant posterior means.
.pop_from_published <- function(spec, model, cond) {
  f <- here("data_published", "bart_rl_online_summary.csv")
  if (!file.exists(f)) stop("Missing published summary: ", f)
  d <- read.csv(f, stringsAsFactors = FALSE)
  cols <- sprintf("%s_%s_%s", model, spec$params, cond)
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop("Published summary lacks: ", paste(miss, collapse = ", "))
  unc <- lapply(seq_along(spec$params), function(k) {
    spec$links[[spec$params[k]]]$to_unc(d[[cols[k]]])
  })
  list(mu    = vapply(unc, function(z) mean(z[is.finite(z)]), numeric(1)),
       sigma = vapply(unc, function(z) sd(z[is.finite(z)]),   numeric(1)))
}

# Truncated normal via inverse CDF (for STL's bounded native-scale parameters).
.rtnorm <- function(n, mu, sd, lo, hi) {
  qnorm(runif(n, pnorm(lo, mu, sd), pnorm(hi, mu, sd)), mu, sd)
}

.draw_true <- function(spec, pop, n) {
  out <- lapply(seq_along(spec$params), function(k) {
    p <- spec$params[k]
    if (!is.null(spec$bounds)) {
      b <- spec$bounds[[p]]
      .rtnorm(n, pop$mu[k], pop$sigma[k], b[1], b[2])
    } else {
      spec$links[[p]]$to_nat(pop$mu[k] + pop$sigma[k] * rnorm(n))
    }
  })
  names(out) <- spec$params
  cbind(data.frame(sub_id = seq_len(n)), as.data.frame(out))
}

# --- one recovery cell ------------------------------------------------------
run_cell <- function(model, cond, trials, opts) {
  spec <- MODELS[[model]]
  if (is.null(spec)) stop("Unknown model: ", model)

  truth <- opts$truth
  fit <- if (truth %in% c("auto", "fit")) .find_fit(model, cond) else NULL
  if (truth == "fit" && is.null(fit))
    stop("RECOVERY_TRUTH=fit but no saved fit found for ", model, "_", cond)
  pop <- if (!is.null(fit)) .pop_from_fit(spec, fit) else .pop_from_published(spec, model, cond)
  src <- if (!is.null(fit)) "fit" else "published"

  cdf  <- condition_trials(trials, cond)
  bp   <- condition_breakpoints(cdf)
  cmax <- max(cdf$color_max)

  set.seed(opts$seed)
  true_params <- .draw_true(spec, pop, opts$n)

  log_msg(sprintf("[%s/%s] truth=%s  ntrial=%d  color_max=%d  n=%d",
                  model, cond, src, length(bp), cmax, opts$n))

  sim <- spec$simulate(true_params, bp, color_max = cmax,
                       balloon_color = substr(cond, 1, 1), seed = opts$seed)
  stan_data <- create_stan_params(sim)

  mod <- cmdstan_model(do.call(here, as.list(spec$stan)),
                       cpp_options = list(stan_threads = TRUE))
  rfit <- mod$sample(
    data = stan_data, iter_warmup = opts$warmup, iter_sampling = opts$sample,
    chains = opts$chains, parallel_chains = opts$chains,
    threads_per_chain = opts$threads,
    seed = opts$seed, adapt_delta = 0.9, max_treedepth = opts$treedepth,
    refresh = 200
  )

  recovered <- as.data.frame(lapply(spec$params, function(p) {
    dr   <- rfit$draws(variables = p, format = "df")
    cols <- grep(paste0("^", p, "\\["), names(dr), value = TRUE)
    idx  <- as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", cols))
    cols <- cols[order(idx)]
    vapply(cols, function(cc) mean(dr[[cc]]), numeric(1))
  }))
  names(recovered) <- paste0("rec_", spec$params)

  cmp <- bind_cols(
    true_params %>% select(all_of(spec$params)) %>% rename_with(~ paste0("true_", .x)),
    recovered
  ) %>% mutate(model = model, condition = cond, truth_source = src)

  out_dir <- here("mcmc", "recovery")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_csv <- file.path(out_dir, sprintf("recovery_%s_%s.csv", model, cond))
  write.csv(cmp, out_csv, row.names = FALSE)

  rs <- vapply(spec$params, function(p)
    cor(cmp[[paste0("true_", p)]], cmp[[paste0("rec_", p)]]), numeric(1))
  log_msg(sprintf("[%s/%s] recovery r: %s", model, cond,
                  paste(sprintf("%s=%.2f", spec$params, rs), collapse = "  ")))
  log_msg("wrote", out_csv)
  invisible(cmp)
}

# --- entry point ------------------------------------------------------------
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  arg_model <- if (length(args) >= 1) args[1] else Sys.getenv("RECOVERY_MODEL", "all")
  arg_cond  <- if (length(args) >= 2) args[2] else Sys.getenv("RECOVERY_COND", "all")

  models <- if (identical(arg_model, "all")) names(MODELS) else strsplit(arg_model, ",")[[1]]
  conds  <- if (identical(arg_cond, "all")) names(CONDITION_FILTERS) else strsplit(arg_cond, ",")[[1]]
  bad <- setdiff(models, names(MODELS))
  if (length(bad)) stop("Unknown model(s): ", paste(bad, collapse = ", "))
  bad <- setdiff(conds, names(CONDITION_FILTERS))
  if (length(bad)) stop("Unknown condition(s): ", paste(bad, collapse = ", "))

  opts <- list(
    truth     = Sys.getenv("RECOVERY_TRUTH", "auto"),
    n         = as.integer(Sys.getenv("RECOVERY_N", "118")),
    warmup    = as.integer(Sys.getenv("RECOVERY_WARMUP", "1000")),
    sample    = as.integer(Sys.getenv("RECOVERY_SAMPLE", "1000")),
    chains    = as.integer(Sys.getenv("RECOVERY_CHAINS", "4")),
    threads   = as.integer(Sys.getenv("RECOVERY_THREADS", "1")),
    seed      = as.integer(Sys.getenv("RECOVERY_SEED", "123")),
    treedepth = as.integer(Sys.getenv("RECOVERY_TREEDEPTH", "10"))
  )
  if (!opts$truth %in% c("auto", "fit", "published"))
    stop("RECOVERY_TRUTH must be auto, fit or published")

  trials <- read.csv(here("data_published", "bart_rl_online_trials.csv"),
                     stringsAsFactors = FALSE)

  log_msg(sprintf("recovery: %d cell(s) = {%s} x {%s}",
                  length(models) * length(conds),
                  paste(models, collapse = ","), paste(conds, collapse = ",")))
  for (m in models) for (cc in conds) run_cell(m, cc, trials, opts)
  log_msg("ALL DONE ->", here("mcmc", "recovery"))
}

if (sys.nframe() == 0L || identical(environment(), globalenv())) main()
