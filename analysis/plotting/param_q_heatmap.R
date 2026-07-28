# analysis/plotting/param_q_heatmap.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(grid)
})

# ---------- shared helpers ----------
.param_cols <- function(draws_df, p) {
  cols <- grep(paste0("^", p, "\\["), names(draws_df), value = TRUE)
  if (!length(cols)) return(character(0))
  idx <- as.integer(sub("^.*\\[(\\d+)\\]$", "\\1", cols))
  cols[order(idx)]
}

# Extract per-participant parameter means from a fit object or draws df
.extract_params_from_fit <- function(fit, fit_name, param_names, qdf, prefix = "ewmv_") {
  # Figure out what 'fit' actually is and get a draws data.frame
  if (inherits(fit, "CmdStanMCMC") || inherits(fit, "CmdStanGQ")) {
    # Proper CmdStanR object
    draws <- fit$draws(variables = param_names, format = "df")
  } else if (is.data.frame(fit)) {
    # Already a draws data.frame
    draws <- fit
  } else if (is.list(fit) && !is.null(fit$draws_df) && is.data.frame(fit$draws_df)) {
    # Custom wrapper: list(draws_df = ...)
    draws <- fit$draws_df
  } else {
    stop(
      sprintf(
        "In .extract_params_from_fit: fit '%s' is of class [%s] and has no usable $draws() or draws_df.",
        fit_name, paste(class(fit), collapse = ", ")
      )
    )
  }
  
  # ---- ALIGNMENT FIX (2026-07) ----------------------------------------------
  # Parameter column [k] belongs to the k-th subject IN THE FIT'S INDEX ORDER,
  # which create_stan_params() set as arrange(participant_id, trial_number) ->
  # distinct(participant_id): i.e. GLOBAL ALPHABETICAL participant_id. Callers
  # build qdf from the combined data file, whose rows are recruitment-blocked
  # (arrange(recruitment, participant_id)), so qdf can arrive in a DIFFERENT order.
  # Assigning means[k] to qdf row k by position (below) is only correct once qdf
  # is sorted into the fit's index order; otherwise param[k] is attached to the
  # wrong participant (was true for 117/118 rows since the 2025-11-25 recruitment
  # merge). Sorting here also guarantees the returned row order == fit draw
  # column order, which downstream posterior-propagation code relies on.
  if (!"participant_id" %in% names(qdf))
    stop(".extract_params_from_fit: qdf must contain participant_id to align to ",
         "the fit's index order")
  qdf <- dplyr::arrange(qdf, participant_id)

  n_q <- nrow(qdf)
  out <- qdf %>% dplyr::mutate(.fit_name = fit_name)
  
  for (p in param_names) {
    cols    <- .param_cols(draws, p)
    n_fit   <- length(cols)
    col_out <- paste0(prefix, p)
    
    if (n_fit == 0) {
      warning(sprintf("[%s] No columns for '%s' — skipping", fit_name, p))
      out[[col_out]] <- NA_real_
      next
    }
    
    n_use <- min(n_q, n_fit)
    
    means <- vapply(
      seq_len(n_use),
      function(j) {
        col_data <- draws[[cols[j]]]
        if (length(col_data) > 0) mean(col_data, na.rm = TRUE) else NA_real_
      },
      numeric(1)
    )
    
    out[[col_out]] <- NA_real_
    out[[col_out]][seq_len(n_use)] <- means
  }
  
  out
}

# --- Correlation test helper (Pearson) ---
.cor_test <- function(x, y, method = "pearson") {
  complete_idx <- complete.cases(x, y)
  x_complete <- x[complete_idx]; y_complete <- y[complete_idx]
  n_complete <- length(x_complete)
  if (n_complete < 3 ||
      length(unique(x_complete)) < 2 ||
      length(unique(y_complete)) < 2 ||
      sd(x_complete, na.rm = TRUE) == 0 ||
      sd(y_complete, na.rm = TRUE) == 0) {
    return(tibble(r = NA_real_, p = NA_real_, n = n_complete))
  }
  test <- tryCatch(
    stats::cor.test(x_complete, y_complete, method = method, exact = FALSE),
    error = function(e) NA_real_
  )
  if (length(test) == 1 && is.na(test)) {
    return(tibble(r = NA_real_, p = NA_real_, n = n_complete))
  }
  tibble(r = as.numeric(test$estimate), p = as.numeric(test$p.value), n = n_complete)
}

# --- plotting helper: PARAM strip over three columns; Blue/Orange/Pink as top ticks ---
# expects df_subset to already contain: r, r_label, Block, q_label, param_label, sig (logical)
.plot_heatmap_subset <- function(df_subset, model_name, param_order, out_dir,
                                 show_title = TRUE,
                                 reverse_y = TRUE,     # reverse alphabetical for y
                                 y_levels = NULL,      # explicit y order (bottom->top); overrides alphabetical
                                 suffix = NULL,
                                 width_per_subcol = 1.0,
                                 height_per_row = 0.46,
                                 base_margin = 3.2,
                                 text_size = 4.6,
                                 dpi_png = 600) {

  # explicit y order if supplied, otherwise alphabetize (ascending or descending)
  y_levels_alpha <- if (!is.null(y_levels)) {
    y_levels
  } else if (reverse_y) {
    rev(sort(unique(df_subset$q_label)))
  } else {
    sort(unique(df_subset$q_label))
  }

  # Condition top-labels: coloured pre / post / con (high-value balloon colours)
  block_labs <- c(
    Blue   = "<span style='color:#2630F5'>**pre**</span>",
    Orange = "<span style='color:#E68D33'>**post**</span>",
    Pink   = "<span style='color:#CF3160'>**con**</span>"
  )
  df_subset <- df_subset %>%
    mutate(
      param_label = factor(param_label, levels = param_order),
      Block       = factor(Block, levels = c("Blue","Orange","Pink"),
                           labels = block_labs[c("Blue","Orange","Pink")]),
      q_label     = factor(q_label, levels = y_levels_alpha)
    )
  
  n_params <- length(param_order)
  plot_width  <- n_params * 3 * width_per_subcol + base_margin
  plot_height <- length(y_levels_alpha) * height_per_row + 3.2
  
  p_heatmap <- ggplot(df_subset, aes(x = Block, y = q_label, fill = r)) +
    geom_tile(color = "white", linewidth = 0.5) +
    # overlay border for significant cells
    geom_tile(data = ~subset(.x, sig %in% TRUE),
              fill = NA, color = "black", linewidth = 0.75) +
    geom_text(aes(label = r_label), size = text_size, color = "black") +
    scale_fill_gradient2(
      low = "#075AFF", mid = "white", high = "#FF0000",
      midpoint = 0, limits = c(-1, 1), name = "r", na.value = "grey80"
    ) +
    scale_x_discrete(position = "top", drop = FALSE) +
    labs(x = NULL, y = NULL) +  # remove y-axis title
    facet_grid(
      rows = vars(), cols = vars(param_label),
      scales = "free_x", space = "free_x", switch = "x"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.placement     = "outside",
      strip.text.x        = ggtext::element_markdown(face = "bold", size = 12.5, margin = margin(b = 6)),
      strip.background.x  = element_rect(fill = NA, colour = NA),
      axis.text.x.top     = ggtext::element_markdown(size = 11.5),
      axis.text.y         = ggtext::element_markdown(size = 10.5, hjust = 1),
      panel.grid          = element_blank(),
      legend.position     = "right",
      panel.spacing.x     = unit(0, "pt")
    )
  
  if (show_title) {
    p_heatmap <- p_heatmap +
      ggtitle(paste(model_name, "Model")) +
      theme(
        plot.title = element_text(face = "bold", size = 15.5,
                                  hjust = 0.5, margin = margin(b = 10))
      )
  }
  
  file_suffix <- if (!is.null(suffix)) paste0("_", suffix) else ""
  plot_fname_base <- file.path(out_dir, paste0("param_q_heatmap_", model_name, file_suffix))
  
  # save SVG + high-DPI PNG
  ggsave(paste0(plot_fname_base, ".svg"), p_heatmap,
         width = plot_width, height = plot_height, bg = "white", limitsize = FALSE)
  ggsave(paste0(plot_fname_base, ".png"), p_heatmap,
         width = plot_width, height = plot_height, dpi = dpi_png, units = "in",
         bg = "white", limitsize = FALSE)
  
  print(p_heatmap)
}

# ---------- Combined Correlation Heatmap ----------
# Holm correction matches your scatter wrappers:
# per (Model, base_param, questionnaire) across the 3 colors only (m = 3).

run_param_q_heatmap <- function(
    df_blue_pre, df_orange_post, df_pink,
    ewmv_blue, ewmv_orange, ewmv_pink,
    stl_blue, stl_orange, stl_pink,
    fourpar_blue, fourpar_orange, fourpar_pink,
    q_cols,
    out_dir = here::here("reporting","figures","param_q","heatmap")
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # --- parameter maps (HTML for headers) ---
  param_map_ewmv <- c(
    phi    = "Prior Weight (&psi;)",
    rho    = "Risk Preference (&rho;)",
    lambda = "Loss Aversion (&lambda;)",
    tau    = "Inverse Temperature (&tau;)",
    eta    = "Updating Exponent (&xi;)"
  ); param_names_ewmv <- names(param_map_ewmv)
  
  param_map_stl <- c(
    vwin     = "Reward Learning Rate (v<sub>win</sub>)",
    vloss    = "Loss Learning Rate (v<sub>loss</sub>)",
    beta     = "Behavioural Consistency (&beta;)",
    omegaone = "Initial Pump Target (&omega;<sub>1</sub>)"
  ); param_names_stl <- names(param_map_stl)
  
  param_map_fourpar <- c(
    phi = "Prior Weight (&phi;)",
    eta = "Updating Exponent (&eta;)",
    gam = "Risk/Curvature (&gamma;)",
    tau = "Inverse Temperature (&tau;)"
  ); param_names_fourpar <- names(param_map_fourpar)
  
  # --- questionnaire bases (one row per participant per block) ---
  qdf_blue   <- df_blue_pre    %>% distinct(participant_id, .keep_all = TRUE)
  qdf_orange <- df_orange_post %>% distinct(participant_id, .keep_all = TRUE)
  qdf_pink   <- df_pink        %>% distinct(participant_id, .keep_all = TRUE)
  
  # --- extract parameters from the fits you pass in ---
  dfb_ewmv <- .extract_params_from_fit(ewmv_blue,   "ewmv_blue",   param_names_ewmv,  qdf_blue,   prefix = "ewmv_")
  dfo_ewmv <- .extract_params_from_fit(ewmv_orange, "ewmv_orange", param_names_ewmv,  qdf_orange, prefix = "ewmv_")
  dfp_ewmv <- .extract_params_from_fit(ewmv_pink,   "ewmv_pink",   param_names_ewmv,  qdf_pink,   prefix = "ewmv_")
  
  dfb_stl  <- .extract_params_from_fit(stl_blue,    "stl_blue",    param_names_stl,   qdf_blue,   prefix = "stl_")
  dfo_stl  <- .extract_params_from_fit(stl_orange,  "stl_orange",  param_names_stl,   qdf_orange, prefix = "stl_")
  dfp_stl  <- .extract_params_from_fit(stl_pink,    "stl_pink",    param_names_stl,   qdf_pink,   prefix = "stl_")
  
  dfb_fourpar <- .extract_params_from_fit(fourpar_blue,   "fourpar_blue",   param_names_fourpar, qdf_blue,   prefix = "fourpar_")
  dfo_fourpar <- .extract_params_from_fit(fourpar_orange, "fourpar_orange", param_names_fourpar, qdf_orange, prefix = "fourpar_")
  dfp_fourpar <- .extract_params_from_fit(fourpar_pink,   "fourpar_pink",   param_names_fourpar, qdf_pink,   prefix = "fourpar_")
  
  # --- combine & long ---
  df_all <- bind_rows(
    dfb_ewmv, dfo_ewmv, dfp_ewmv,
    dfb_stl,  dfo_stl,  dfp_stl,
    dfb_fourpar, dfo_fourpar, dfp_fourpar
  )
  
  param_cols_all <- c(
    paste0("ewmv_",    param_names_ewmv),
    paste0("stl_",     param_names_stl),
    paste0("fourpar_", param_names_fourpar)
  )
  
  df_long <- df_all %>%
    pivot_longer(cols = all_of(param_cols_all),
                 names_to = "param_name", values_to = "param_value") %>%
    filter(!is.na(param_value)) %>%
    pivot_longer(cols = all_of(q_cols),
                 names_to = "q_name", values_to = "q_value")
  
  # --- correlations (raw r, raw p) ---
  df_cor <- df_long %>%
    group_by(.fit_name, param_name, q_name) %>%
    reframe(.cor_test(param_value, q_value)) %>%
    ungroup()
  
  # --- add Model/Block + base_param labels ---
  df_labeled <- df_cor %>%
    mutate(
      Model = case_when(
        grepl("^ewmv_", .fit_name)    ~ "EWMV",
        grepl("^stl_", .fit_name)     ~ "STL",
        grepl("^fourpar_", .fit_name) ~ "4PAR",
        TRUE ~ NA_character_
      ),
      Block = case_when(
        grepl("blue", .fit_name)   ~ "Blue",
        grepl("orange", .fit_name) ~ "Orange",
        grepl("pink", .fit_name)   ~ "Pink",
        TRUE ~ NA_character_
      ),
      base_param = sub("^(ewmv|stl|fourpar)_", "", param_name)
    )
  
  # --- HOLM ONLY (per Model × base_param × q over 3 colors) ---
  df_labeled_holm <- df_labeled %>%
    group_by(Model, base_param, q_name) %>%
    mutate(p_adj = p.adjust(p, method = "holm")) %>%
    ungroup()
  
  build_df_plot <- function(src_df, mode = c("raw","holm")) {
    mode <- match.arg(mode)
    p_use_sym <- if (mode == "holm") quote(p_adj) else quote(p)
    
    src_df %>%
      mutate(
        param_label = case_when(
          Model == "EWMV"    & base_param == "phi"    ~ param_map_ewmv["phi"],
          Model == "EWMV"    & base_param == "rho"    ~ param_map_ewmv["rho"],
          Model == "EWMV"    & base_param == "lambda" ~ param_map_ewmv["lambda"],
          Model == "EWMV"    & base_param == "tau"    ~ param_map_ewmv["tau"],
          Model == "EWMV"    & base_param == "eta"    ~ param_map_ewmv["eta"],
          Model == "STL"     & base_param == "vwin"   ~ param_map_stl["vwin"],
          Model == "STL"     & base_param == "vloss"  ~ param_map_stl["vloss"],
          Model == "STL"     & base_param == "beta"   ~ param_map_stl["beta"],
          Model == "STL"     & base_param == "omegaone" ~ param_map_stl["omegaone"],
          Model == "4PAR" & base_param == "phi"    ~ param_map_fourpar["phi"],
          Model == "4PAR" & base_param == "eta"    ~ param_map_fourpar["eta"],
          Model == "4PAR" & base_param == "gam"    ~ param_map_fourpar["gam"],
          Model == "4PAR" & base_param == "tau"    ~ param_map_fourpar["tau"],
          TRUE ~ NA_character_
        ),
        param_order_within_model = case_when(
          Model == "EWMV"    ~ match(base_param, names(param_map_ewmv)),
          Model == "STL"     ~ match(base_param, names(param_map_stl)),
          Model == "4PAR" ~ match(base_param, names(param_map_fourpar)),
          TRUE ~ NA_integer_
        ),
        q_label   = toupper(gsub("_"," ", q_name)),
        p_use_val = !!p_use_sym,
        p_stars   = dplyr::case_when(
          is.na(p_use_val)  ~ "",
          p_use_val < 0.001 ~ "***",
          p_use_val < 0.01  ~ "**",
          p_use_val < 0.05  ~ "*",
          TRUE              ~ ""
        ),
        r_label = ifelse(is.na(r), "NA", sprintf("%.2f%s", r, p_stars)),
        sig     = !is.na(p_use_val) & (p_use_val < 0.05),
        Model   = factor(Model, levels = c("EWMV","STL","4PAR")),
        Block   = factor(Block, levels = c("Blue","Orange","Pink"))
      )
  }
  
  # ONLY HOLM
  df_plot_holm <- build_df_plot(df_labeled_holm,  "holm")
  
  order_for <- function(dfm) {
    dfm %>% arrange(param_order_within_model) %>% pull(param_label) %>% unique()
  }
  param_order_ewmv_holm    <- order_for(filter(df_plot_holm, Model == "EWMV"))
  param_order_stl_holm     <- order_for(filter(df_plot_holm, Model == "STL"))
  param_order_fourpar_holm <- order_for(filter(df_plot_holm, Model == "4PAR"))
  
  # --- export heatmaps (HOLM ONLY; one per model; labeled exactly like you had it) ---
  if (nrow(filter(df_plot_holm, Model == "EWMV")) > 0) {
    .plot_heatmap_subset(
      filter(df_plot_holm, Model == "EWMV"),
      "EWMV", param_order_ewmv_holm, out_dir,
      show_title = TRUE, reverse_y = TRUE, suffix = "holm"
    )
  }
  
  if (nrow(filter(df_plot_holm, Model == "STL")) > 0) {
    .plot_heatmap_subset(
      filter(df_plot_holm, Model == "STL"),
      "STL", param_order_stl_holm, out_dir,
      show_title = TRUE, reverse_y = TRUE, suffix = "holm"
    )
  }
  
  if (nrow(filter(df_plot_holm, Model == "4PAR")) > 0) {
    .plot_heatmap_subset(
      filter(df_plot_holm, Model == "4PAR"),
      "4PAR", param_order_fourpar_holm, out_dir,
      show_title = TRUE, reverse_y = TRUE, suffix = "holm"
    )
  }
  # --- PRINT SIGNIFICANT HOLM-CORRECTED CORRELATIONS TO CONSOLE ---
  sig_tbl <- df_plot_holm %>%
    filter(sig %in% TRUE) %>%
    transmute(
      Model,
      Block,
      Parameter = base_param,
      Questionnaire = q_name,
      r = round(r, 3),
      p_holm = signif(p_use_val, 3),
      n
    ) %>%
    arrange(Model, Parameter, Questionnaire, Block)
  
  cat("\n\n=== SIGNIFICANT HOLM-CORRECTED PARAM × QUESTIONNAIRE CORRELATIONS (p < 0.05) ===\n")
  if (nrow(sig_tbl) == 0) {
    cat("(none)\n")
  } else {
    print(sig_tbl, n = nrow(sig_tbl))
  }
  # ----------------------------------------------------------------
  invisible(list(holm = df_plot_holm))
}