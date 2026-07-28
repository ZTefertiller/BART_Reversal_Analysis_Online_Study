# analysis/plotting/ewmv_recovery_plot.R
# True vs recovered EWMV parameters from modeling/ewmv/run_ewmv_recovery.R.
# ewmv_recovery_grid(): 4 conditions (rows) x 5 parameters (cols) scatter grid
# with identity + fit lines and a recovery-correlation table.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

.REC_PARAMS <- c(
  phi    = "Prior Weight (ψ)",
  eta    = "Updating Exp. (ξ)",
  rho    = "Risk Pref. (ρ)",
  tau    = "Inv. Temp. (τ)",
  lambda = "Loss Aversion (λ)"
)

.REC_CONDS <- c(
  blue   = "pre reversal",
  orange = "post reversal",
  yellow = "constant",
  pink   = "control"
)

# Reads whatever recovery_<cond>.csv files exist; returns NULL if none.
ewmv_recovery_grid <- function(
    reg_dir   = here::here("mcmc", "ewmv_recovery"),
    conds     = names(.REC_CONDS),
    family    = "Arial",
    base_size = 12
) {
  pn <- names(.REC_PARAMS)

  dfs <- lapply(conds, function(cc) {
    f <- file.path(reg_dir, paste0("recovery_", cc, ".csv"))
    if (!file.exists(f)) return(NULL)
    d <- read.csv(f); d$condition <- cc; d
  })
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) return(NULL)

  long <- bind_rows(lapply(dfs, function(d) {
    bind_rows(lapply(pn, function(p) {
      tibble(condition = d$condition[1], param = p,
             true = d[[paste0("true_", p)]], rec = d[[paste0("rec_", p)]])
    }))
  })) %>%
    filter(is.finite(true), is.finite(rec)) %>%
    # Standardize each parameter (by the TRUE mean/SD) so every panel shares one
    # scale and the identity line sits at 45 degrees. r is scale-invariant.
    group_by(condition, param) %>%
    mutate(tz = as.numeric(scale(true)),
           rz = (rec - mean(true, na.rm = TRUE)) / sd(true, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      cond_lab  = factor(.REC_CONDS[condition],
                         levels = .REC_CONDS[intersect(conds, unique(condition))]),
      param_lab = factor(.REC_PARAMS[param], levels = .REC_PARAMS[pn])
    )

  # Per-condition point colours (match the rest of the deck).
  cond_cols <- c(blue = "#2630F5", orange = "#E68D33",
                 yellow = "#C9A800", pink = "#CF3160")

  cors <- long %>%
    group_by(condition, cond_lab, param, param_lab) %>%
    summarise(r = cor(true, rec), n = dplyr::n(), .groups = "drop")

  # r annotation: plain black text, top-left of each panel (no fill chip).
  ann <- cors %>%
    mutate(r_lab = sprintf("r = %.2f", r), x = -2.85, y = 2.85)

  LIM <- 3
  grid <- ggplot(long, aes(tz, rz)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey45", linewidth = 0.7) +
    geom_point(aes(colour = condition), alpha = 0.40, size = 1.9,
               show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE, colour = "grey15",
                linewidth = 1.4, formula = y ~ x) +
    geom_text(data = ann, aes(x = x, y = y, label = r_lab),
              hjust = 0, vjust = 1, size = base_size / 2.3, fontface = "bold",
              colour = "black", family = family) +
    scale_colour_manual(values = cond_cols) +
    # condition strips on the RIGHT (vertical), so the left y-axis title sits
    # close to the panels; parameter strips on top.
    facet_grid(cond_lab ~ param_lab) +
    coord_fixed(ratio = 1, xlim = c(-LIM, LIM), ylim = c(-LIM, LIM)) +
    scale_x_continuous(breaks = c(-2, 0, 2)) +
    scale_y_continuous(breaks = c(-2, 0, 2)) +
    labs(x = "True parameter value", y = "Recovered parameter value") +
    theme_minimal(base_size = base_size, base_family = family) +
    theme(panel.grid.minor = element_blank(),
          panel.spacing = unit(0.8, "lines"),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = 0.4),
          strip.text.y = element_text(angle = -90, face = "bold",
                                      size = base_size, family = family),
          strip.text.x = element_text(face = "bold", size = base_size, family = family),
          strip.placement = "outside",
          axis.title.y = element_text(size = base_size, family = family,
                                      margin = margin(r = 2)),
          axis.title = element_text(size = base_size, family = family),
          axis.text = element_text(size = base_size - 1, colour = "grey45",
                                   family = family))

  # wide table: rows = condition, cols = parameters (in canonical order), cells = r
  tbl <- cors %>%
    select(Condition = cond_lab, param_lab, r) %>%
    pivot_wider(names_from = param_lab, values_from = r) %>%
    select(Condition, all_of(unname(.REC_PARAMS[pn]))) %>%
    mutate(across(where(is.numeric), ~ sprintf("%.2f", .x)))

  list(grid = grid, table = tbl, cors = cors,
       conditions = unique(long$condition))
}
