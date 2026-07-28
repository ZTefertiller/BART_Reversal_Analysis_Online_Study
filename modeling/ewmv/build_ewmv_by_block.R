# Rebuild modeling/ewmv/ewmv_by_block.csv from the model fits.
#
# This intermediate is .gitignored (*.csv) and does not travel with the repo.
# The published equivalent is data_published/bart_rl_online_ewmv_params.csv,
# which is what the analysis report reads; this script is how that file is made.
#
# Source of truth = mcmc/model_fits_list.rds (the same fits the report uses).
#
# ALIGNMENT: parameter column [k] belongs to the k-th subject in the fit's index
# order, which create_stan_params() sets as arrange(sub_id, trial_number) -> the
# k-th distinct sub_id. .extract_params_from_fit() sorts qdf by sub_id so the
# params land on the right person, and we arrange the final table by sub_id so
# its row order == the fit draw column order (required by the
# posterior-propagation step in ewmv_trait_robustness.R).
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

data <- read.csv(here("data_published", "bart_rl_online_trials.csv"))
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
qdf_blue   <- data %>% filter(balloon_color == "b") %>% distinct(sub_id, .keep_all = TRUE)
qdf_orange <- data %>% filter(balloon_color == "o") %>% distinct(sub_id, .keep_all = TRUE)
qdf_pink   <- data %>% filter(balloon_color == "p") %>% distinct(sub_id, .keep_all = TRUE)

pre <- .extract_params_from_fit(fit_blue, "ewmv_blue", param_names,
                                qdf_blue %>% select(sub_id),
                                prefix = "pre_")
post <- .extract_params_from_fit(fit_orange, "ewmv_orange", param_names,
                                 qdf_orange %>% select(sub_id),
                                 prefix = "post_")
con <- .extract_params_from_fit(fit_pink, "ewmv_pink", param_names,
                                qdf_pink %>% select(sub_id),
                                prefix = "con_")

q_cols <- c("spq_total", "spq_b_total", "spq_cogper", "spq_interp", "spq_disorg",
            "pdi_total", "pdi_distress", "pdi_frequency", "pdi_conviction",
            "caps_total", "caps_distress", "caps_intrusiveness", "caps_frequency",
            "phq_total", "mdq_total", "ipip_total", "ppgm_total")
q_cols <- intersect(q_cols, names(data))

questionnaires <- data %>%
  distinct(sub_id, .keep_all = TRUE) %>%
  select(sub_id, all_of(q_cols))

out <- pre %>%
  select(sub_id, starts_with("pre_")) %>%
  full_join(post %>% select(sub_id, starts_with("post_")),
            by = "sub_id") %>%
  full_join(con %>% select(sub_id, starts_with("con_")),
            by = "sub_id") %>%
  left_join(questionnaires, by = "sub_id") %>%
  arrange(sub_id)   # row order == fit draw column order (alphabetical)

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
