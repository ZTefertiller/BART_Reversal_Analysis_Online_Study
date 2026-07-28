#!/usr/bin/env Rscript
# Sort MCMC draws into a clear, self-describing layout:
#
#   mcmc/<model>/<color>/<date>_<model>_<color>_n<P>_warmup<W>_draws<D>.rds
#   mcmc/<model>/<color>/<...>_summary.rds
#
# Each .rds is a self-contained cmdstanr fit (save_object materialised the draws),
# so the raw cmdstan chain CSVs are redundant backups and are left in
# mcmc/fits/ (or pruned later). Idempotent: re-run as more fits finish.
#
# Naming carries: date, model/fit, n participants, warmup, and sampling draws.

SRC  <- "mcmc/fits"
DATE <- format(Sys.Date())
N    <- 118L           # participants
W    <- 1000L          # warmup iterations
D    <- 2000L          # sampling draws (per chain; 4 chains)

fits <- list.files(SRC, pattern = "^(stl|ewmv|fourpar)_(blue|orange|pink)\\.rds$",
                   recursive = TRUE, full.names = TRUE)
if (!length(fits)) { cat("No finished fits to sort yet.\n"); quit(save = "no") }

for (f in fits) {
  nm    <- sub("\\.rds$", "", basename(f))          # e.g. stl_blue
  parts <- strsplit(nm, "_")[[1]]
  model <- parts[1]; color <- parts[2]
  dest  <- file.path("mcmc", model, color)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  base  <- sprintf("%s_%s_%s_n%d_warmup%d_draws%d", DATE, model, color, N, W, D)

  file.rename(f, file.path(dest, paste0(base, ".rds")))
  sm <- file.path(dirname(f), paste0(nm, "_summary.rds"))
  if (file.exists(sm)) file.rename(sm, file.path(dest, paste0(base, "_summary.rds")))
  cat("sorted", nm, "->", file.path(dest, paste0(base, ".rds")), "\n")
}
