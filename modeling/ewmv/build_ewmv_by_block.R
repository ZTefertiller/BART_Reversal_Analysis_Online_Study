# Rebuild modeling/ewmv/ewmv_by_block.csv from the CANONICAL model fits.
#
# This intermediate is .gitignored (*.csv) and does not travel with the repo.
# reporting/BART_presentation.Rmd reads it for the "Schizotypy and Inverse
# Temperature" section.
#
# Source of truth = mcmc/model_fits_list.rds (the SAME fits used by
# reporting/BART_analysis_report.Rmd and its cut-down param x questionnaire
# heatmaps).
#
# ALIGNMENT (2026-07): parameter column [k] belongs to the k-th subject in the
# fit's index order == GLOBAL ALPHABETICAL participant_id (create_stan_params:
# arrange(participant_id, trial_number) -> distinct(participant_id)). The
# combined data file is recruitment-blocked, so params must be attached in the fit's
# order, NOT the file order. .extract_params_from_fit() now sorts qdf by
# participant_id, so the params land on the right person; we also arrange the
# final table by participant_id so its row order == the fit draw column order
# (required by the posterior-propagation step in ewmv_trait_robustness.R).
#
# NOTE: an earlier version of this file assumed the aligned tau x SPQ ~ -.12 was
# a "DIFFERENT EWMV run" and rebuilt to force -.23. That -.23 was the SAME fits
# joined by row position to the wrong participants; -.12 is the correct value.
#
#   pre_*  <- ewmv_blue   (blue  = pre-reversal high-value balloon)
#   post_* <- ewmv_orange (orange = post-reversal high-value balloon)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

source(here("analysis", "plotting", "param_q_heatmap.R"))  # .extract_params_from_fit, .cor_test
source(here("analysis", "plotting", "spq_factors.R"))       # add_spq_factors

param_names <- c("phi", "eta", "rho", "tau", "lambda")

data <- read.csv(here("data", "processed", "filtered_data", "full_all_recruitments.csv"))
data <- add_spq_factors(data)

fits <- readRDS(here("mcmc", "model_fits_list.rds"))
get_fit <- function(pattern) {
  idx <- grep(pattern, names(fits), ignore.case = TRUE)
  if (length(idx) != 1) stop("Pattern '", pattern, "' matched ", length(idx), " fits")
  fits[[idx]]
}
fit_blue   <- get_fit("^ewmv.*blue")
fit_orange <- get_fit("^ewmv.*orange")
fit_pink   <- get_fit("^ewmv.*pink")

# Per-participant rows in the SAME order the report's heatmap uses.
qdf_blue   <- data %>% filter(balloon_color == "b") %>% distinct(participant_id, .keep_all = TRUE)
qdf_orange <- data %>% filter(balloon_color == "o") %>% distinct(participant_id, .keep_all = TRUE)
qdf_pink   <- data %>% filter(balloon_color == "p") %>% distinct(participant_id, .keep_all = TRUE)

pre <- .extract_params_from_fit(fit_blue, "ewmv_blue", param_names,
                                qdf_blue %>% select(participant_id),
                                prefix = "pre_")
post <- .extract_params_from_fit(fit_orange, "ewmv_orange", param_names,
                                 qdf_orange %>% select(participant_id),
                                 prefix = "post_")
con <- .extract_params_from_fit(fit_pink, "ewmv_pink", param_names,
                                qdf_pink %>% select(participant_id),
                                prefix = "con_")

q_cols <- c("spq_total", "spq_b_total", "spq_cogper", "spq_interp", "spq_disorg",
            "pdi_total", "pdi_distress", "pdi_frequency", "pdi_conviction",
            "caps_total", "caps_distress", "caps_intrusiveness", "caps_frequency",
            "phq_total", "mdq_total", "ipip_total", "ppgm_total")
q_cols <- intersect(q_cols, names(data))

questionnaires <- data %>%
  distinct(participant_id, .keep_all = TRUE) %>%
  select(participant_id, all_of(q_cols))

out <- pre %>%
  select(participant_id, starts_with("pre_")) %>%
  full_join(post %>% select(participant_id, starts_with("post_")),
            by = "participant_id") %>%
  full_join(con %>% select(participant_id, starts_with("con_")),
            by = "participant_id") %>%
  left_join(questionnaires, by = "participant_id") %>%
  arrange(participant_id)   # row order == fit draw column order (alphabetical)

out_path <- here("modeling", "ewmv", "ewmv_by_block.csv")
write_csv(out, out_path)
cat("Wrote", out_path, "-", nrow(out), "participants,", ncol(out), "columns\n")

# Sanity check on a headline correlation.
ct <- .cor_test(out$spq_total, out$post_tau)
cat(sprintf("CHECK  SPQ Total x post_tau:  r = %.3f  (p = %.3f, n = %d)\n",
            ct$r, ct$p, ct$n))
cti <- .cor_test(out$spq_interp, out$post_tau)
cat(sprintf("CHECK  SPQ Interp x post_tau: r = %.3f  (p = %.3f)\n", cti$r, cti$p))
ctp <- .cor_test(out$pdi_total, out$post_tau)
cat(sprintf("CHECK  PDI Total x post_tau:  r = %.3f  (p = %.3f)\n", ctp$r, ctp$p))
