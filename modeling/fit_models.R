#!/usr/bin/env Rscript
# ===========================================================================
# Fit the three hierarchical BART models (STL, 4PAR, EWMV) separately for each
# task phase (blue = pre-reversal, orange = post-reversal, pink = control) to
# all 118 participants. Subjects are indexed by sub_id (1..118).
#
# Sampler: 4 chains x 1000 warmup / 2000 sampling.
# Writes per-fit .rds objects, per-fit summaries, and a diagnostics table under
# mcmc/fits/. Run modeling/build_model_fits_list.R afterwards to assemble
# mcmc/model_fits_list.rds, which the analysis report reads.
# ===========================================================================
suppressPackageStartupMessages({
  library(dplyr); library(cmdstanr); library(posterior)
  library(here); library(readr); library(tibble)
})

OUT <- here("mcmc", "fits")
for (m in c("", "stl", "ewmv", "fourpar"))
  dir.create(file.path(OUT, m), recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")

# create_stan_params() and the condition filters are shared with
# modeling/run_recovery.R so the fits and the recovery refits are built from one
# definition.
source(here("modeling", "stan_data.R"))

df <- read.csv(here("data_published", "bart_rl_online_trials.csv"),
               stringsAsFactors = FALSE)

blue   <- create_stan_params(condition_trials(df, "blue"))
orange <- create_stan_params(condition_trials(df, "orange"))
pink   <- create_stan_params(condition_trials(df, "pink"))
log_msg("stan params built: nsub per color =",
        blue$nsub, orange$nsub, pink$nsub, "| ntrial =", blue$ntrial)

cmd_stan_fit <- function(data, fit_name, model_path, out_sub) {
  fit_name <- paste0(fit_name, "_n", data$nsub)
  mod <- cmdstan_model(model_path, cpp_options = list(stan_threads = TRUE))
  t0 <- Sys.time()
  fit <- mod$sample(
    data = data, iter_warmup = 1000, iter_sampling = 2000,
    chains = 4, parallel_chains = 4, threads_per_chain = 2,
    seed = 10191998, adapt_delta = 0.99, step_size = 1, max_treedepth = 12,
    refresh = 250, save_warmup = FALSE,
    output_dir = file.path(OUT, out_sub))
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  fit$save_object(file.path(OUT, out_sub, paste0(fit_name, ".rds")))
  ds  <- fit$diagnostic_summary()
  sm  <- fit$summary()
  diag <- tibble(
    fit = fit_name, minutes = round(mins, 2),
    num_divergent = sum(ds$num_divergent),
    num_max_treedepth = sum(ds$num_max_treedepth),
    ebfmi_min = round(min(ds$ebfmi), 3),
    rhat_max = round(max(sm$rhat, na.rm = TRUE), 4),
    ess_bulk_min = round(min(sm$ess_bulk, na.rm = TRUE), 0),
    ess_tail_min = round(min(sm$ess_tail, na.rm = TRUE), 0))
  log_msg(sprintf("DONE %-14s %.1f min | div=%d | rhat_max=%.3f | ess_bulk_min=%.0f",
                  fit_name, mins, diag$num_divergent, diag$rhat_max, diag$ess_bulk_min))
  saveRDS(sm, file.path(OUT, out_sub, paste0(fit_name, "_summary.rds")))
  diag
}

specs <- tribble(
  ~name,            ~data,   ~model,                                  ~dir,
  "stl_blue",       "blue",  here("modeling","stl","stl_vec.stan"),   "stl",
  "stl_orange",     "orange",here("modeling","stl","stl_vec.stan"),   "stl",
  "stl_pink",       "pink",  here("modeling","stl","stl_vec.stan"),   "stl",
  "fourpar_blue",   "blue",  here("modeling","fourpar","4par.stan"),  "fourpar",
  "fourpar_orange", "orange",here("modeling","fourpar","4par.stan"),  "fourpar",
  "fourpar_pink",   "pink",  here("modeling","fourpar","4par.stan"),  "fourpar",
  "ewmv_blue",      "blue",  here("modeling","ewmv","ewmv_vec.stan"), "ewmv",
  "ewmv_orange",    "orange",here("modeling","ewmv","ewmv_vec.stan"), "ewmv",
  "ewmv_pink",      "pink",  here("modeling","ewmv","ewmv_vec.stan"), "ewmv")

data_env <- list(blue = blue, orange = orange, pink = pink)
log_msg("STARTING 9 fits (4 chains x 1000/2000)")
diags <- bind_rows(lapply(seq_len(nrow(specs)), function(i) {
  s <- specs[i, ]
  log_msg("fitting", s$name, "...")
  cmd_stan_fit(data_env[[s$data]], s$name, s$model, s$dir)
}))
write_csv(diags, file.path(OUT, "diagnostics_summary.csv"))
log_msg("ALL DONE. diagnostics ->", file.path(OUT, "diagnostics_summary.csv"))
print(as.data.frame(diags))
