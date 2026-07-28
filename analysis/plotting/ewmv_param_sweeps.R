# analysis/plotting/ewmv_param_sweeps.R
# R port of the computational-model parameter-sweep figures (originally
# pug_presenation_scripts/bart_param_figures.py). Same Park et al. (2021)
# equations, same sweep values / ranges, but native ggplot so they always
# regenerate here and use "pump"/"pumps"/"explosion" wording (never "inflation").
#
# Public: ewmv_param_sweeps(family, base_size) -> named list of ggplots:
#   ewmv_xi, ewmv_psi, ewmv_lambda, ewmv_rho, ewmv_tau,
#   par4_eta, par4_gamma, compare_prior, compare_risk

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(ggtext)
})

# ---- model equations (verbatim from Park et al. 2021) ----------------------
.eps_belief  <- function(N, psi, xi, P) { w <- exp(-xi * N); w * psi + (1 - w) * P }
.eps_utility <- function(l, p, lam, rho, r = 1) {
  (1 - p) * r - p * lam * (l - 1) * r + rho * p * (1 - p) * (r + lam * (l - 1) * r)^2
}
.eps_ppump   <- function(l, p, lam, rho, tau, r = 1)
  1 / (1 + exp(-tau * .eps_utility(l, p, lam, rho, r)))
.eps_par4_weight  <- function(N, eta) 1 / (1 + eta * N)
.eps_par4_belief  <- function(N, phi, eta, P) {
  w <- .eps_par4_weight(N, eta); w * (1 - phi) + (1 - w) * P
}
.eps_par4_utility <- function(l, p, gamma, r = 1) (1 - p)^l * (l * r)^gamma

# evenly spaced colours along a viridis-family ramp (high->low like the original)
# Constrain the ramp away from the near-white yellow end (end = 0.82) so the
# lightest line in plasma/magma/viridis stays a readable gold rather than glaring.
.eps_cols <- function(n, option = "viridis")
  viridisLite::viridis(n, option = option, begin = 0.05, end = 0.82)

.eps_theme <- function(base_size, family) {
  theme_minimal(base_size = base_size, base_family = family) +
    theme(plot.title = element_text(face = "bold", family = family,
                                    size = base_size + 1),
          legend.title = element_blank(),
          legend.text = element_text(size = base_size - 2, family = family),
          legend.position = "inside",
          legend.position.inside = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                           colour = NA),
          axis.title = element_text(size = base_size, family = family),
          axis.text = element_text(size = base_size - 2, family = family),
          panel.grid.minor = element_blank())
}

# generic single-sweep line plot from a long df (x, y, lvl ordered hi->lo)
.eps_sweep <- function(d, title, xlab, ylab, cols, base_size, family,
                       hline = NULL, hline_lty = "dotted", ylim = NULL) {
  p <- ggplot(d, aes(x, y, colour = lvl, group = lvl)) +
    { if (!is.null(hline))
        geom_hline(yintercept = hline, linetype = hline_lty, colour = "grey50",
                   linewidth = 0.6) } +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = cols) +
    labs(title = title, x = xlab, y = ylab) +
    .eps_theme(base_size, family)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

# helper: build long df sweeping one parameter
.eps_long <- function(values, xs, fun, fmt) {
  bind_rows(lapply(values, function(v)
    data.frame(x = xs, y = fun(v), lvl = sprintf(fmt, v)))) %>%
    mutate(lvl = factor(lvl, levels = sprintf(fmt, values)))
}

# Strip a sweep plot down to a clean bordered box to pair with the contrast
# individual boxes: full black border, small inside legend, no panel grid.
# `title` (HTML markdown, e.g. "Prior Weight (&psi;)") becomes a bold strip-style
# heading on top — this is where the PARAMETER NAME lives in the combined grid,
# so the paired contrast box below can drop its own strip. Axis titles/text are
# scaled up (title_size / axis_size) to read at the same size as that heading.
sweep_to_box <- function(p, base_size = 13, family = "Arial",
                         border_lwd = 1.4, title = NULL,
                         title_size = base_size + 4,
                         axis_size  = base_size,
                         tick_size  = base_size - 1,
                         legend_size = base_size - 2) {
  p <- p +
    labs(title = title) +
    theme_classic(base_size = base_size, base_family = family) +
    theme(
      text         = element_text(family = family),
      axis.title   = element_text(size = axis_size, family = family,
                                  face = "bold", colour = "black"),
      axis.text    = element_text(size = tick_size, family = family,
                                  face = "bold", colour = "black"),
      axis.line    = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA,
                                  linewidth = border_lwd),
      legend.title = element_blank(),
      legend.text  = element_text(size = legend_size, family = family,
                                  colour = "black"),
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.key.size = unit(1.25, "lines"),
      legend.background = element_rect(fill = scales::alpha("white", 0.7),
                                       colour = NA),
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 8))
  if (is.null(title)) {
    p + theme(plot.title = element_blank())
  } else {
    p + theme(plot.title = ggtext::element_markdown(
      size = title_size, face = "bold", family = family, colour = "black",
      hjust = 0.5, margin = margin(b = 4)))
  }
}

ewmv_param_sweeps <- function(family = "Arial", base_size = 13,
                              by_block_csv = here::here("modeling", "ewmv",
                                                        "ewmv_by_block.csv"),
                              mono_option = NULL) {
  out <- list()

  # Colour helper: by default each sweep keeps its own viridis-family palette
  # (the `option` argument). When `mono_option` is set, EVERY sweep uses that one
  # ramp so higher vs lower parameter values are coloured consistently across all
  # panels (the "shared colouring" variant). The line order in each sweep is
  # high->low, so the ramp direction is identical everywhere.
  cf <- function(n, option) .eps_cols(n, if (is.null(mono_option)) option
                                          else mono_option)

  # --- EWMV xi : belief vs cumulative pumps -------------------------------
  xi_v <- c(0.05, 0.02, 0.01, 0.005, 0.002); N <- seq(0, 600, length.out = 400)
  out$ewmv_xi <- .eps_sweep(
    .eps_long(xi_v, N, function(v) .eps_belief(N, psi = 0.30, xi = v, P = 0.10),
              "ξ = %g"),
    "ξ  —  updating exponent", "Cumulative pumps",
    "Believed explosion probability", cf(length(xi_v), "viridis"),
    base_size, family, hline = 0.10, ylim = c(0, max(0.30, 0.10) * 1.08))

  # --- EWMV psi : belief vs cumulative pumps ------------------------------
  psi_v <- c(0.70, 0.50, 0.30, 0.15, 0.05)
  out$ewmv_psi <- .eps_sweep(
    .eps_long(psi_v, N, function(v) .eps_belief(N, psi = v, xi = 0.01, P = 0.10),
              "ψ = %g"),
    "ψ  —  prior belief of explosion", "Cumulative pumps",
    "Believed explosion probability", cf(length(psi_v), "cividis"),
    base_size, family, hline = 0.10, ylim = c(0, max(psi_v) * 1.08))

  # --- EWMV lambda : utility vs pump number -------------------------------
  lam_v <- c(8, 4, 2, 1); l <- 1:25
  out$ewmv_lambda <- .eps_sweep(
    .eps_long(lam_v, l, function(v) .eps_utility(l, p = 0.05, lam = v, rho = 0.0),
              "λ = %g"),
    "λ  —  loss aversion", "Pump number", "Utility of pumping",
    cf(length(lam_v), "plasma"), base_size, family,
    hline = 0, hline_lty = "solid", ylim = c(-2, 1.15))

  # --- EWMV rho : utility vs pump number ----------------------------------
  rho_v <- c(0.03, 0.01, 0.001, -0.005, -0.02)
  out$ewmv_rho <- .eps_sweep(
    .eps_long(rho_v, l, function(v) .eps_utility(l, p = 0.05, lam = 4, rho = v),
              "ρ = %g"),
    "ρ  —  risk preference (variance weight)", "Pump number",
    "Utility of pumping", cf(length(rho_v), "viridis"),
    base_size, family, hline = 0, hline_lty = "solid", ylim = c(-4, 3))

  # --- EWMV tau : probability of pumping vs pump number -------------------
  tau_v <- c(8, 3, 1, 0.5)
  out$ewmv_tau <- .eps_sweep(
    .eps_long(tau_v, l, function(v)
      .eps_ppump(l, p = 0.05, lam = 4, rho = 0.0, tau = v), "τ = %g"),
    "τ  —  inverse temperature", "Pump number", "Probability of pumping",
    cf(length(tau_v), "magma"), base_size, family,
    hline = 0.5, ylim = c(0, 1.02))

  # --- Four-parameter eta : belief vs cumulative pumps --------------------
  eta_v <- c(0.02, 0.01, 0.005, 0.002, 0.001)
  out$par4_eta <- .eps_sweep(
    .eps_long(eta_v, N, function(v) .eps_par4_belief(N, phi = 0.70, eta = v, P = 0.10),
              "η = %g"),
    "η  —  updating coefficient (Four-Parameter)", "Cumulative pumps",
    "Believed explosion probability", .eps_cols(length(eta_v), "viridis"),
    base_size, family, hline = 0.10,
    ylim = c(0, if ((1 - 0.70) > 0.10) (1 - 0.70) * 1.1 else 0.10 * 1.3))

  # --- Four-parameter gamma : normalised EU vs pump number ----------------
  gam_v <- c(1.5, 1.0, 0.6, 0.3); lg <- 1:40
  out$par4_gamma <- .eps_sweep(
    .eps_long(gam_v, lg, function(v) {
      u <- .eps_par4_utility(lg, p = 0.05, gamma = v); u / max(u) }, "γ = %g"),
    "γ  —  risk-taking propensity (Four-Parameter)", "Pump number",
    "Expected utility (normalised)", .eps_cols(length(gam_v), "plasma"),
    base_size, family)

  # --- comparison: weight on prior (hyperbolic vs exponential) -----------
  x <- seq(0, 6, length.out = 400)
  cmp <- bind_rows(
    data.frame(x = x, y = 1 / (1 + x),
               lvl = "Four-Parameter:  ω = 1 / (1 + x)  (hyperbolic)"),
    data.frame(x = x, y = exp(-x),
               lvl = "EWMV:  ω = e^(−x)  (exponential)"))
  out$compare_prior <- ggplot(cmp, aes(x, y, colour = lvl)) +
    geom_line(linewidth = 1.3) +
    scale_colour_manual(values = c("#c1666b", "#3b6ea5")) +
    coord_cartesian(ylim = c(0, 1.02)) +
    labs(title = "Weight on the prior belief",
         x = "x  =  (updating rate) × cumulative pumps",
         y = "Weight on the prior  (ω)") +
    .eps_theme(base_size, family)

  # --- comparison: risk representation (gamma vs lambda/rho) --------------
  gl <- 1:40
  cmpL <- .eps_long(c(1.2, 0.6, 0.3), gl, function(v) {
    u <- .eps_par4_utility(gl, p = 0.05, gamma = v); u / max(u) }, "γ = %g")
  pL <- .eps_sweep(cmpL, "Four-Parameter:  risk = curvature (γ)", "Pump number",
                   "Expected utility (normalised)", .eps_cols(3, "plasma"),
                   base_size, family)
  l2 <- 1:25
  cmpR <- bind_rows(
    data.frame(x = l2, y = .eps_utility(l2, 0.05, 2, 0.0),  lvl = "baseline  (λ=2, ρ=0)"),
    data.frame(x = l2, y = .eps_utility(l2, 0.05, 6, 0.0),  lvl = "higher loss aversion  (λ=6)"),
    data.frame(x = l2, y = .eps_utility(l2, 0.05, 2, 0.02), lvl = "higher variance weight  (ρ=0.02)")) %>%
    mutate(lvl = factor(lvl, levels = c("baseline  (λ=2, ρ=0)",
                                        "higher loss aversion  (λ=6)",
                                        "higher variance weight  (ρ=0.02)")))
  pR <- ggplot(cmpR, aes(x, y, colour = lvl)) +
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.6) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = c(
      "baseline  (λ=2, ρ=0)"             = "#404040",
      "higher loss aversion  (λ=6)"      = "#c1666b",
      "higher variance weight  (ρ=0.02)" = "#3b6ea5")) +
    coord_cartesian(ylim = c(-2.5, 1.6)) +
    labs(title = "EWMV:  risk = loss aversion (λ) + variance (ρ)",
         x = "Pump number", y = "Utility of pumping") +
    .eps_theme(base_size, family)
  out$compare_risk <- pL + pR

  out
}
