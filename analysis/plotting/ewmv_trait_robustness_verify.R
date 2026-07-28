# analysis/plotting/ewmv_trait_robustness_verify.R
# STEP 0 of the EWMV trait-parameter robustness checks: REPRODUCE the published
# Pearson correlation matrix exactly, using the same pipeline as
# analysis/plotting/ewmv_cutdown_heatmap.R::ewmv_cutdown_heatmap():
#   - same source CSV (modeling/ewmv/ewmv_by_block.csv, built from the full fits)
#   - same parameter set (.ECH_PARAMS) and condition prefixes (pre_/post_/con_)
#   - same questionnaire rows (.ECH_Q)
#   - same Pearson .cor_test (complete.cases, exact = FALSE)
#   - same Holm correction: per (param, q) across the three conditions (m = 3)
#   - same participant set (whatever rows the CSV carries; .cor_test drops NA pairwise)
#
# Values are the per-participant EWMV parameter means joined to trait scores.
#
# This script ONLY reads + reproduces; it modifies nothing. Source it and call
# ewmv_trait_repro_matrix(); it returns the long df and prints a wide table.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr)
})

# Reuse the EXACT helpers / definitions from the published pipeline so there is
# zero chance of a divergent parameter list, prefix, or correlation routine.
source(here::here("analysis", "plotting", "param_q_heatmap.R"))     # .cor_test
source(here::here("analysis", "plotting", "ewmv_cutdown_heatmap.R")) # .ECH_PARAMS/.ECH_Q/.ech_cor_long

# Build the full long correlation table (param x q x condition) + Holm p_adj,
# exactly mirroring ewmv_cutdown_heatmap() + .ech_plot()'s holm step.
ewmv_trait_long <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    method = "pearson") {
  d <- utils::read.csv(by_block_csv)
  long <- bind_rows(
    .ech_cor_long(d, "pre_",  "pre reversal",  method),
    .ech_cor_long(d, "post_", "post reversal", method),
    .ech_cor_long(d, "con_",  "control",       method))
  # Holm per (param, q) across the conditions shown — identical to .ech_plot().
  long %>%
    group_by(param, q) %>%
    mutate(p_adj = p.adjust(p, method = "holm")) %>%
    ungroup()
}

# Print the reproduced matrix as r (p_holm) wide tables, one per condition, and
# run explicit assertions against the reported headline numbers.
ewmv_trait_repro_matrix <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    method = "pearson", tol = 0.005, verbose = TRUE) {

  long <- ewmv_trait_long(by_block_csv, method)

  if (verbose) {
    for (cd in c("pre reversal", "post reversal", "control")) {
      cat(sprintf("\n=== %s — r (p_holm) ===\n", cd))
      wide <- long %>%
        filter(cond == cd) %>%
        mutate(cell = sprintf("%+.2f (%.3f)", r, p_adj),
               q_label = .ECH_Q[q], param_label = .ECH_PARAMS[param]) %>%
        select(q_label, param_label, cell) %>%
        pivot_wider(names_from = param_label, values_from = cell)
      print(as.data.frame(wide), row.names = FALSE)
    }
  }

  # --- assertions against the reported headline values ---------------------
  get_r <- function(p, q, cd) {
    v <- long$r[long$param == p & long$q == q & long$cond == cd]
    if (length(v) != 1) NA_real_ else v
  }
  checks <- tibble::tribble(
    ~label,                              ~param, ~q,           ~cond,           ~expected,
    "SPQ Interp x tau (pre)",            "tau",  "spq_interp", "pre reversal",  -0.13,
    "SPQ Interp x tau (post)",           "tau",  "spq_interp", "post reversal", -0.08,
    "SPQ Interp x tau (control)",        "tau",  "spq_interp", "control",       -0.08
  ) %>%
    mutate(actual = purrr::pmap_dbl(list(param, q, cond),
                                    ~ get_r(..1, ..2, ..3)),
           diff = abs(actual - expected),
           ok = diff <= tol)

  cat("\n=== REPRODUCTION CHECK vs reported values ===\n")
  print(as.data.frame(checks %>% select(label, expected, actual, diff, ok)),
        row.names = FALSE)

  all_ok <- all(checks$ok, na.rm = TRUE) && !any(is.na(checks$actual))
  if (all_ok) {
    cat("\n*** MATCH: reproduced r values agree with reported values (tol =",
        tol, "). Safe to proceed with robustness checks. ***\n")
  } else {
    cat("\n*** MISMATCH: reproduced r values DO NOT match reported values.",
        "STOP and investigate before running robustness checks. ***\n")
  }

  invisible(list(long = long, checks = checks, all_ok = all_ok))
}

if (sys.nframe() == 0) {
  ewmv_trait_repro_matrix()
}
