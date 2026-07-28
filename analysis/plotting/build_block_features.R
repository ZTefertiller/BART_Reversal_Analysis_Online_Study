#!/usr/bin/env Rscript
# analysis/build_block_features.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
})

# -------------------------------
# Core: append block-wise averages
# -------------------------------
append_block_averages <- function(df,
                                  block_size = 10L,
                                  blocks_per_color = c(b = 6L, o = 6L, y = 6L, p = 3L),
                                  verbose = TRUE,
                                  sanity_out = "output/block_sanity_checks.csv") {
  req <- c("participant_id","balloon_color","trial_number","inflations","popped")
  stopifnot(all(req %in% names(df)))
  
  # Robust popped flag (logical)
  popped_flag <- function(x) tolower(as.character(x)) %in% c("1","true","t","yes","y")
  
  # 1) Order and assign **within-color** blocks of size 10 by appearance
  df_blk <- df %>%
    arrange(participant_id, balloon_color, trial_number) %>%
    group_by(participant_id, balloon_color) %>%
    mutate(.row_within_color = row_number(),
           block = ceiling(.row_within_color / block_size)) %>%
    ungroup() %>%
    mutate(popped_logical = popped_flag(popped))
  
  # 2) Keep only requested number of blocks per color
  allowed <- tibble(
    balloon_color = names(blocks_per_color),
    max_block     = as.integer(blocks_per_color)
  )
  
  df_blk <- df_blk %>%
    left_join(allowed, by = "balloon_color") %>%
    filter(!is.na(max_block), block <= max_block) %>%
    select(-.row_within_color)
  
  # 3) Build block membership tables (trial numbers used)
  #    - trials_all: the 10 trials in the block
  #    - trials_adj: the subset of those trials where popped == FALSE
  membership <- df_blk %>%
    group_by(participant_id, balloon_color, block) %>%
    summarise(
      trials_all = list(trial_number),
      trials_adj = list(trial_number[!popped_logical]),
      n_all      = length(trial_number),
      n_adj      = sum(!popped_logical),
      .groups = "drop"
    ) %>%
    mutate(
      block_label = paste0(balloon_color, block),
      trials_all_str = vapply(trials_all, function(v) paste(v, collapse = ","), character(1)),
      trials_adj_str = vapply(trials_adj, function(v) paste(v, collapse = ","), character(1))
    ) %>%
    select(participant_id, balloon_color, block, block_label, n_all, n_adj, trials_all_str, trials_adj_str)
  
  # 4) Console sanity checks
  if (verbose) {
    cat("\n--- Block sanity checks (per participant / color / block) ---\n")
    membership %>%
      arrange(participant_id, balloon_color, block) %>%
      group_split(participant_id) %>%
      walk(function(tt) {
        pid <- tt$participant_id[1]
        cat(sprintf("\n[participant %s]\n", pid))
        for (i in seq_len(nrow(tt))) {
          r <- tt[i, ]
          cat(sprintf("  %s: n_all=%d  n_adj=%d\n    trials_all: %s\n    trials_adj: %s\n",
                      r$block_label, r$n_all, r$n_adj, r$trials_all_str, r$trials_adj_str))
        }
      })
  }
  
  # 5) Write sanity checks to CSV (one row per participant/color/block)
  if (!is.null(sanity_out)) {
    dir.create(dirname(sanity_out), recursive = TRUE, showWarnings = FALSE)
    write_csv(membership, sanity_out)
  }
  
  # 6) Compute participant-level means per (color, block)
  #    avg_all  = mean inflations on **all** trials in that block
  #    avg_adj  = mean inflations on the **same block trials** but with pops excluded
  means_all <- df_blk %>%
    group_by(participant_id, balloon_color, block) %>%
    summarise(mean_all = mean(inflations, na.rm = TRUE), .groups = "drop")
  
  means_adj <- df_blk %>%
    group_by(participant_id, balloon_color, block) %>%
    summarise(mean_adj = mean(inflations[!popped_logical], na.rm = TRUE), .groups = "drop")
  
  # 7) Wide feature names: e.g., b1_avg, b1_avg_adj, o6_avg, o6_avg_adj
  wide_all <- means_all %>%
    mutate(var = paste0(balloon_color, block, "_avg")) %>%
    select(participant_id, var, mean_all) %>%
    pivot_wider(names_from = var, values_from = mean_all)
  
  wide_adj <- means_adj %>%
    mutate(var = paste0(balloon_color, block, "_avg_adj")) %>%
    select(participant_id, var, mean_adj) %>%
    pivot_wider(names_from = var, values_from = mean_adj)
  
  # 8) Join wide features back to every row of the participant (by design)
  out <- df %>%
    left_join(wide_all, by = "participant_id") %>%
    left_join(wide_adj, by = "participant_id")
  
  # 9) Return both the augmented data and the membership table (invisibly)
  attr(out, "block_membership") <- membership
  out
}

# ---------------------------------
# CLI usage (optional, can be sourced)
# ---------------------------------
# Run from command line:
#   Rscript analysis/build_block_features.R /path/to/input.csv /path/to/output.csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  in_csv  <- args[1]
  out_csv <- if (length(args) >= 2) args[2] else "output/data_with_block_avgs.csv"
  
  stopifnot(file.exists(in_csv))
  dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
  
  data <- readr::read_csv(in_csv, show_col_types = FALSE)
  
  data2 <- append_block_averages(
    data,
    block_size = 10L,
    blocks_per_color = c(b = 6L, o = 6L, y = 6L, p = 3L),
    verbose = TRUE,
    sanity_out = "output/block_sanity_checks.csv"
  )
  
  readr::write_csv(data2, out_csv)
  cat(sprintf("\nWrote augmented data to: %s\nSanity checks: output/block_sanity_checks.csv\n", out_csv))
}

