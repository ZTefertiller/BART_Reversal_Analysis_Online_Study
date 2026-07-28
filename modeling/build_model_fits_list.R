#!/usr/bin/env Rscript
# ===========================================================================
# build_model_fits_list.R
# ---------------------------------------------------------------------------
# Assemble mcmc/model_fits_list.rds -- the named list of the 9 fits
# (<model>_<color>) that the analysis report reads.
#
# Reads the per-fit .rds files written by fit_models.R and writes the combined
# named list. If several files match a slot the NEWEST is used.
#
#   Rscript modeling/build_model_fits_list.R [FITS_DIR]
#   (default FITS_DIR = mcmc/fits)
# ===========================================================================
suppressPackageStartupMessages(library(here))
args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else here("mcmc", "fits")
OUT  <- here("mcmc", "model_fits_list.rds")

models <- c("stl", "fourpar", "ewmv")
colors <- c("blue", "orange", "pink")

find_one <- function(model, color) {
  dir  <- file.path(ROOT, model)
  hits <- list.files(dir, pattern = sprintf("_%s_%s_.*\\.rds$", model, color),
                     full.names = TRUE)
  if (length(hits) == 0L) stop("no fit for ", model, "_", color, " under ", dir)
  hits[order(file.mtime(hits), decreasing = TRUE)][1L]   # newest if several
}

fits <- list()
for (m in models) for (col in colors) {
  nm <- paste0(m, "_", col)
  f  <- find_one(m, col)
  cat("loading", nm, "<-", basename(f), "\n")
  fits[[nm]] <- readRDS(f)
}
stopifnot(length(fits) == 9L)
saveRDS(fits, OUT)
cat("\nwrote", OUT, "\n  ", length(fits), "fits:", paste(names(fits), collapse = ", "), "\n")
