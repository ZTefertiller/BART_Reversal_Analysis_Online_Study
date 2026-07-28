# analysis/05_looic_comparison.R
# LOOIC comparisons using in-memory fits ONLY.
# Assumes the following objects ALREADY EXIST in the environment:
# stl_blue, stl_orange, stl_pink
# fourpar_blue, fourpar_orange, fourpar_pink
# ewmv_blue, ewmv_orange, ewmv_pink

suppressPackageStartupMessages({
  library(cmdstanr)
  library(loo)
  library(posterior)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(officer)
  library(writexl)
  library(purrr)
})

RESULTS_DIR <- here::here("output", "model_comparison")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

# ===================================================================
# LOG-LIK EXTRACTION
# ===================================================================
.get_loglik <- function(fit, name) {
  if (is.null(fit)) {
    warning("Fit ", name, " is NULL")
    return(NULL)
  }
  
  ll <- try(fit$draws("log_lik", format = "matrix"), silent = TRUE)
  
  if (inherits(ll, "try-error") || is.null(ll) || ncol(ll) == 0) {
    warning("No log_lik found for: ", name)
    return(NULL)
  }
  
  r_eff <- try(
    relative_eff(
      exp(ll),
      chain_id = rep(1:fit$num_chains(), each = fit$num_iters())
    ),
    silent = TRUE
  )
  
  if (inherits(r_eff, "try-error")) r_eff <- NA
  
  list(loglik = ll, r_eff = r_eff)
}

# ===================================================================
# CALCULATE LOO FOR A LIST OF FITS
# ===================================================================
calculate_loo_list <- function(fit_list) {
  out <- list()
  for (nm in names(fit_list)) {
    z <- .get_loglik(fit_list[[nm]], nm)
    if (is.null(z)) next
    
    out[[nm]] <- tryCatch(
      loo(z$loglik, r_eff = z$r_eff),
      error = function(e) {
        warning("LOO failed for ", nm, ": ", e$message)
        NULL
      }
    )
  }
  out
}

# ===================================================================
# FINAL TABLE CREATION
# ===================================================================
create_final_loo_table <- function(loo_list) {
  if (!length(loo_list)) return(tibble())
  
  model_names <- names(loo_list)
  
  looic_vals <- map_dbl(
    loo_list,
    ~ .x$estimates["looic", "Estimate"],
    .default = NA_real_
  )
  
  df <- tibble(Model = model_names, LOOIC = looic_vals)
  
  cmp <- NULL
  if (length(loo_list) > 1) {
    cmp <- tryCatch(loo_compare(loo_list), error = function(e) NULL)
  }
  
  if (!is.null(cmp)) {
    cmp_df <- as.data.frame(cmp)
    cmp_df$Model <- rownames(cmp_df)
    df <- df |>
      left_join(cmp_df |> select(Model, elpd_diff), by = "Model") |>
      mutate(dLOOIC = -2 * elpd_diff) |>
      select(-elpd_diff)
  } else {
    df$dLOOIC <- 0
  }
  
  # ---- ADD LOO MODEL WEIGHTS ----
  if (length(loo_list) > 1) {
    w <- tryCatch(
      loo_model_weights(loo_list, method = "pseudobma", BB = TRUE),
      error = function(e) NULL
    )
    if (!is.null(w)) {
      df <- df |> mutate(`LOO weight` = unname(w[Model]))
    } else {
      df <- df |> mutate(`LOO weight` = NA_real_)
    }
  } else {
    df <- df |> mutate(`LOO weight` = 1)
  }
  # ------------------------------
  
  pk <- tibble(
    Model = model_names,
    `Max Pareto k` = map_dbl(
      loo_list,
      ~ max(.x$diagnostics$pareto_k, na.rm = TRUE),
      .default = NA_real_
    ),
    `% k > 0.7` = map_dbl(
      loo_list,
      ~ mean(.x$diagnostics$pareto_k > 0.7, na.rm = TRUE) * 100,
      .default = NA_real_
    )
  )
  
  df <- df |>
    left_join(pk, by = "Model") |>
    arrange(LOOIC) |>
    mutate(across(where(is.numeric), round, 2))
  
  df
}

# ===================================================================
# WORD EXPORT
# ===================================================================
save_table_to_word <- function(df, title, fname) {
  if (!nrow(df)) return(invisible(NULL))
  
  doc <- read_docx()
  doc <- body_add_par(doc, title, style = "heading 2")
  doc <- body_add_table(doc, value = df)
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Pareto-k interpretation:", style = "heading 3")
  doc <- body_add_par(doc, "k > 0.7 indicates unreliable LOO estimates; compare models cautiously.")
  
  print(doc, target = file.path(RESULTS_DIR, fname))
}

# ===================================================================
# PLOTTING
# ===================================================================
plot_loo_loss <- function(cmp_table, title) {
  if (is.null(cmp_table) || !"elpd_diff" %in% names(cmp_table)) return(NULL)
  
  df <- cmp_table |>
    mutate(
      elpd_loss = -elpd_diff,
      se_diff   = ifelse(is.na(se_diff), 0, se_diff)
    )
  
  ggplot(df, aes(x = reorder(model, elpd_loss), y = elpd_loss)) +
    geom_col(fill = "lightblue") +
    geom_errorbar(
      aes(
        ymin = pmax(0, elpd_loss - se_diff),
        ymax = elpd_loss + se_diff
      ),
      width = 0.2
    ) +
    coord_flip() +
    labs(
      title = title,
      x     = "Model",
      y     = "ELPD Loss vs Best Model"
    ) +
    theme_minimal(base_family = "Arial")
}

# ===================================================================
# CLEAN REPORTING HELPERS (combined table + dot/CI plot)
# ===================================================================

# Tidy table across the three conditions, adding the SE of ΔLOOIC (= 2*se_diff
# from loo_compare) so differences are interpretable.
# Older cached loo_res objects carry the delta column as a literal "ΔLOOIC".
# Whether that name matches the UTF-8 literal in this file depends on the locale
# the cache was written under, so match it structurally instead: any column that
# is not one of the known ASCII ones is the delta column. Everything downstream
# uses the ASCII name dLOOIC.
.normalise_looic_tbl <- function(tbl) {
  if (is.null(tbl) || !nrow(tbl)) return(tbl)
  if ("dLOOIC" %in% names(tbl)) return(tbl)
  known <- c("Model", "LOOIC", "LOO weight", "Max Pareto k", "% k > 0.7")
  cand <- setdiff(names(tbl), known)
  if (length(cand) == 1L) names(tbl)[names(tbl) == cand] <- "dLOOIC"
  else stop("build_looic_combined: cannot identify the delta-LOOIC column; ",
            "got columns: ", paste(names(tbl), collapse = ", "))
  tbl
}

build_looic_combined <- function(loo_res) {
  conds <- c(Blue = "blue", Orange = "orange", Pink = "pink")
  purrr::map_dfr(names(conds), function(cond) {
    tbl <- .normalise_looic_tbl(loo_res[[conds[cond]]])
    if (is.null(tbl) || !nrow(tbl)) return(NULL)
    cmp <- loo_res$cmp[[conds[cond]]]
    se <- if (!is.null(cmp) && "se_diff" %in% names(cmp)) {
      tibble(Model = cmp$model,
             dLOOIC_se = 2 * ifelse(is.na(cmp$se_diff), 0, cmp$se_diff))
    } else {
      tibble(Model = tbl$Model, dLOOIC_se = NA_real_)
    }
    tbl |>
      left_join(se, by = "Model") |>
      mutate(Condition = factor(cond, levels = c("Blue", "Orange", "Pink")))
  })
}

# Clean ΔLOOIC dot + CI plot (0 = best model; higher = worse), faceted by block.
plot_looic_combined <- function(loo_res) {
  d <- build_looic_combined(loo_res)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d |>
    mutate(Model = factor(Model, levels = c("4PAR", "STL", "EWMV")))
  cond_cols <- c(Blue = "#2630F5", Orange = "#E68D33", Pink = "#CF3160")
  ggplot(d, aes(x = dLOOIC, y = Model, colour = Condition)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = pmax(0, dLOOIC - dLOOIC_se),
                       xmax = dLOOIC + dLOOIC_se),
                   height = 0, linewidth = 0.7, na.rm = TRUE) +
    geom_point(size = 2.8) +
    facet_wrap(~ Condition, nrow = 1) +
    scale_colour_manual(values = cond_cols, guide = "none") +
    labs(x = "ΔLOOIC vs. best model  (0 = best; lower is better)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold"))
}

# ===================================================================
# MAIN WRAPPER (in-memory fits ONLY)
# ===================================================================
looic_compare_all <- function() {
  
  blue_fits <- list(
    STL     = stl_blue,
    `4PAR` = fourpar_blue,
    EWMV    = ewmv_blue
  ) |> purrr::discard(is.null)
  
  orange_fits <- list(
    STL     = stl_orange,
    `4PAR` = fourpar_orange,
    EWMV    = ewmv_orange
  ) |> purrr::discard(is.null)
  
  pink_fits <- list(
    STL     = stl_pink,
    `4PAR` = fourpar_pink,
    EWMV    = ewmv_pink
  ) |> purrr::discard(is.null)
  
  # Compute LOOs
  blue_loos   <- calculate_loo_list(blue_fits)
  orange_loos <- calculate_loo_list(orange_fits)
  pink_loos   <- calculate_loo_list(pink_fits)
  
  # Create tables
  blue_tbl   <- create_final_loo_table(blue_loos)
  orange_tbl <- create_final_loo_table(orange_loos)
  pink_tbl   <- create_final_loo_table(pink_loos)
  
  # Export tables
  save_table_to_word(blue_tbl,   "Blue LOOIC",   "looic_blue.docx")
  save_table_to_word(orange_tbl, "Orange LOOIC", "looic_orange.docx")
  save_table_to_word(pink_tbl,   "Pink LOOIC",   "looic_pink.docx")
  
  if (requireNamespace("writexl", quietly = TRUE)) {
    write_xlsx(
      list(Blue = blue_tbl, Orange = orange_tbl, Pink = pink_tbl),
      file.path(RESULTS_DIR, "looic_tables.xlsx")
    )
  }
  
  # Comparison tables (always compute)
  blue_cmp_raw   <- tryCatch(loo_compare(blue_loos),   error = function(e) NULL)
  orange_cmp_raw <- tryCatch(loo_compare(orange_loos), error = function(e) NULL)
  pink_cmp_raw   <- tryCatch(loo_compare(pink_loos),   error = function(e) NULL)
  
  blue_cmp_table <- if (!is.null(blue_cmp_raw)) {
    df <- as.data.frame(blue_cmp_raw)
    df$model <- rownames(df)
    df
  } else NULL
  
  orange_cmp_table <- if (!is.null(orange_cmp_raw)) {
    df <- as.data.frame(orange_cmp_raw)
    df$model <- rownames(df)
    df
  } else NULL
  
  pink_cmp_table <- if (!is.null(pink_cmp_raw)) {
    df <- as.data.frame(pink_cmp_raw)
    df$model <- rownames(df)
    df
  } else NULL
  
  # PRINT THEM
  cat("\n\n=== BLUE MODEL COMPARISON TABLE ===\n");   print(blue_cmp_table)
  cat("\n\n=== ORANGE MODEL COMPARISON TABLE ===\n"); print(orange_cmp_table)
  cat("\n\n=== PINK MODEL COMPARISON TABLE ===\n");   print(pink_cmp_table)
  
  # Build plots from the comparison tables
  plots <- list(
    blue   = if (!is.null(blue_cmp_table))   plot_loo_loss(blue_cmp_table,   "Blue LOOIC")   else NULL,
    orange = if (!is.null(orange_cmp_table)) plot_loo_loss(orange_cmp_table, "Orange LOOIC") else NULL,
    pink   = if (!is.null(pink_cmp_table))   plot_loo_loss(pink_cmp_table,   "Pink LOOIC")   else NULL
  )
  
  # RETURN EVERYTHING
  list(
    blue   = blue_tbl,
    orange = orange_tbl,
    pink   = pink_tbl,
    cmp    = list(
      blue   = blue_cmp_table,
      orange = orange_cmp_table,
      pink   = pink_cmp_table
    ),
    plots  = plots
  )
}