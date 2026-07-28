# analysis/plotting/ewmv_cutdown_heatmap.R
# Self-contained cut-down EWMV parameter x questionnaire correlation heatmaps for
# the conference presentation. Cut-down rows = the schizotypy/PLE measures only:
#   CAPS Total, PDI Total, SPQ Total, SPQ Cognitive-Perceptual, SPQ Interpersonal,
#   SPQ Disorganized.
# Columns = the five EWMV parameters; facets = condition (pre / post / control).
# Larger, readable fonts (Arial); base_size controls overall scale.
#
# Two builders:
#   ewmv_cutdown_heatmap()          main fits, conditions pre / post / control,
#                                   read from modeling/ewmv/ewmv_by_block.csv.
#   ewmv_cutdown_heatmap_nofirst()  no-first-trial fits (blue=pre, orange=post),
#                                   read from the 1 MB *_no_first.rds + task.csv.
#
# Significance: Holm-corrected per (parameter x questionnaire) across the
# conditions shown; outlined cells + stars mark p < .05.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})
source(here::here("analysis", "plotting", "param_q_heatmap.R"))  # .extract_params_from_fit, .cor_test
source(here::here("analysis", "plotting", "spq_factors.R"))       # add_spq_factors

.ECH_PARAMS <- c(phi = "Prior weight (ψ)", eta = "Updating exp. (ξ)",
                 rho = "Risk pref. (ρ)", tau = "Inv. temp. (τ)",
                 lambda = "Loss aversion (λ)")
# cut-down questionnaire rows (top -> bottom)
.ECH_Q <- c(caps_total = "CAPS Total", pdi_total = "PDI Total",
            spq_total  = "SPQ Total",  spq_cogper = "SPQ Cog-Per",
            spq_interp = "SPQ Interp", spq_disorg = "SPQ Disorg")
.ECH_COND_COL <- c("pre reversal" = "#2630F5", "post reversal" = "#E68D33",
                   "control" = "#CF3160")

# Build the long correlation df from a per-participant data.frame that has
# columns <prefix><param> and the questionnaire columns, for one condition.
.ech_cor_long <- function(d, prefix, cond_label, method = "pearson") {
  out <- lapply(names(.ECH_PARAMS), function(p) {
    pc <- paste0(prefix, p)
    lapply(names(.ECH_Q), function(q) {
      cr <- .cor_test(d[[q]], d[[pc]], method = method)
      data.frame(param = p, q = q, cond = cond_label,
                 r = cr$r, p = cr$p, n = cr$n, stringsAsFactors = FALSE)
    }) %>% bind_rows()
  }) %>% bind_rows()
  out
}

# Assemble + style the faceted heatmap from a stacked long df (param,q,cond,r,p).
# cell_mult scales the in-cell r-value text relative to base_size (larger =
# bigger numbers in the tiles). All text is black and uses one family.
.ech_plot <- function(long, cond_levels, base_size = 13, family = "Arial",
                      cell_mult = 0.5) {
  long <- long %>%
    group_by(param, q) %>%
    mutate(p_adj = p.adjust(p, method = "holm")) %>%
    ungroup() %>%
    mutate(
      stars = case_when(is.na(p_adj) ~ "", p_adj < .001 ~ "***",
                        p_adj < .01 ~ "**", p_adj < .05 ~ "*", TRUE ~ ""),
      r_lab = ifelse(is.na(r), "NA", sprintf("%.2f%s", r, stars)),
      sig   = !is.na(p_adj) & p_adj < .05,
      param_label = factor(.ECH_PARAMS[param], levels = unname(.ECH_PARAMS)),
      q_label     = factor(.ECH_Q[q], levels = rev(unname(.ECH_Q))),
      cond        = factor(cond, levels = cond_levels)
    )

  ggplot(long, aes(x = param_label, y = q_label, fill = r)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_tile(data = ~ subset(.x, sig), fill = NA, color = "black",
              linewidth = 1.1) +
    geom_text(aes(label = r_lab), size = base_size * cell_mult, color = "black",
              fontface = "bold", family = family) +
    scale_fill_gradient2(low = "#075AFF", mid = "white", high = "#FF0000",
                         midpoint = 0, limits = c(-1, 1), name = "r",
                         na.value = "grey85") +
    facet_wrap(~ cond, nrow = 1) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(
      text          = element_text(family = family, colour = "black"),
      panel.grid    = element_blank(),
      # condition strips (pre/post/control) sized to MATCH the x/y axis labels.
      strip.text    = element_text(face = "bold", size = base_size + 3,
                                   family = family, colour = "black"),
      axis.text.x   = element_text(size = base_size + 3, angle = 40, hjust = 1,
                                   family = family, colour = "black"),
      axis.text.y   = element_text(size = base_size + 3, family = family,
                                   colour = "black"),
      legend.title  = element_text(size = base_size + 1, family = family,
                                   colour = "black"),
      legend.text   = element_text(size = base_size, family = family,
                                   colour = "black"),
      panel.spacing = unit(1, "lines")
    )
}

# ---- main fits: pre / post / control from ewmv_by_block.csv -----------------
ewmv_cutdown_heatmap <- function(
    by_block_csv = here::here("data_published", "bart_rl_online_ewmv_params.csv"),
    base_size = 13, family = "Arial", method = "pearson", cell_mult = 0.5) {
  d <- utils::read.csv(by_block_csv)
  long <- bind_rows(
    .ech_cor_long(d, "pre_",  "pre reversal",  method),
    .ech_cor_long(d, "post_", "post reversal", method),
    .ech_cor_long(d, "con_",  "control",       method))
  .ech_plot(long, names(.ECH_COND_COL), base_size, family, cell_mult)
}

# ---- no-first-trial fits: blue (pre) + orange (post) ------------------------
# Reads the small *_no_first.rds fits, rebuilds per-participant params in the
# exact Stan order, joins task.csv questionnaires (+ SPQ factors), correlates.
ewmv_cutdown_heatmap_nofirst <- function(
    task_csv   = here::here("data", "processed", "filtered_data", "task.csv"),
    blue_rds   = here::here("mcmc", "ewmv", "ewmv_blue_no_first.rds"),
    orange_rds = here::here("mcmc", "ewmv", "ewmv_orange_post_no_first.rds"),
    base_size = 13, family = "Arial", method = "pearson", cell_mult = 0.5) {
  df <- utils::read.csv(task_csv)
  pn <- names(.ECH_PARAMS)

  # per-participant rows in the SAME order passed to Stan (first trial removed)
  qorder <- function(color, pre) {
    df %>%
      filter(balloon_color == color,
             if (pre) trial_number < 91 else trial_number > 90) %>%
      arrange(sub_id, trial_number) %>%
      group_by(sub_id) %>% slice(-1) %>% ungroup() %>%
      distinct(sub_id, .keep_all = TRUE) %>%
      arrange(sub_id)
  }
  q_per_id <- df %>% distinct(sub_id, .keep_all = TRUE) %>%
    add_spq_factors() %>%
    select(sub_id, all_of(names(.ECH_Q)))

  one <- function(rds, color, pre, cond_label) {
    fit <- readRDS(rds)
    qdf <- qorder(color, pre)
    pp  <- .extract_params_from_fit(fit, "ewmv_nofirst", pn, qdf, prefix = "p_") %>%
      select(sub_id, paste0("p_", pn)) %>%
      left_join(q_per_id, by = "sub_id")
    .ech_cor_long(pp, "p_", cond_label, method)
  }

  long <- bind_rows(
    one(blue_rds,   "b", TRUE,  "pre reversal"),
    one(orange_rds, "o", FALSE, "post reversal"))
  .ech_plot(long, c("pre reversal", "post reversal"), base_size, family,
            cell_mult)
}
