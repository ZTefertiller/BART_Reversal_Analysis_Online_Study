# analysis/tables/compact_summaries.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
})

# ── Public API ────────────────────────────────────────────────────────────────
# build_compact_tables(df, out_prefix = "bart_summary_compact")
#   - Writes:
#       "<out_prefix>_bygender_ttests.xlsx"
#       "<out_prefix>_collapsed.xlsx"
#   - Returns: list(by_gender = <df>, collapsed = <df>)
# Requirements in df:
#   sub_id, gender (male/female/diverse), session (1/2),
#   balloon_color (b/o/y/p), trial_number, inflations, popped (1/0)
#
# The by-gender table and its t-tests compare Men vs Women only; participants
# who reported another gender are kept in the collapsed table but are not one
# of the two compared groups.

build_compact_tables <- function(df, out_prefix = "bart_summary_compact") {
  # Guardrails
  need <- c("sub_id","gender","session","balloon_color",
            "trial_number","inflations","popped")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

  # ── Normalize + compute earnings (do NOT overwrite df) ─────────────────────
  color_map <- c(b = "Blue", o = "Orange", y = "Yellow", p = "Pink")

  dat <- df %>%
    transmute(
      id     = sub_id,
      gender = dplyr::case_when(
        tolower(trimws(as.character(gender))) == "female" ~ "Women",
        tolower(trimws(as.character(gender))) == "male"   ~ "Men",
        TRUE                                              ~ "Other"
      ),
      sess   = as.integer(session),
      color  = recode(balloon_color, !!!color_map),
      trial  = as.integer(trial_number),
      infl   = inflations,
      exp    = popped
    ) %>%
    mutate(
      earn = ifelse(sess == 1, infl * 0.003, infl * 0.01)  # GBP rules
    )
  
  # ── Helpers ────────────────────────────────────────────────────────────────
  m   <- function(x) mean(x, na.rm = TRUE)
  sd_ <- function(x) sd(x, na.rm = TRUE)
  
  pair_cell <- function(M, SD, money = FALSE) {
    if (money) paste0("£", sprintf("%.2f", M), "(", "£", sprintf("%.2f", SD), ")")
    else       paste0(sprintf("%.2f", M),  "(", sprintf("%.2f", SD),  ")")
  }
  
  tt_compact <- function(x_m, x_w) {
    x_m <- x_m[is.finite(x_m)]
    x_w <- x_w[is.finite(x_w)]
    if (length(x_m) < 2 || length(x_w) < 2) return("")
    tt <- try(stats::t.test(x_m, x_w, var.equal = FALSE), silent = TRUE)
    if (inherits(tt, "try-error")) return("")
    tval <- sprintf("%.2f", unname(tt$statistic))
    pval <- ifelse(tt$p.value < .001, "<.001", sprintf("%.3f", tt$p.value))
    paste0("t=", tval, ", p=", pval)
  }
  
  # per-id stats inside a section already filtered to a single (sess,color slice)
  per_id_stats <- function(d) {
    if (nrow(d) == 0) return(tibble())
    core <- d %>%
      group_by(id, gender, color) %>%
      summarise(
        earn_sum  = sum(earn, na.rm = TRUE),
        exp_sum   = sum(exp == 1, na.rm = TRUE),
        adj_total = m(infl[exp == 0]),
        .groups = "drop"
      )
    
    buckets <- d %>%
      arrange(trial) %>%
      group_by(id, gender, color) %>%
      mutate(idx = row_number(),
             bucket = case_when(
               idx <= 10 ~ "First10",
               idx <= 20 ~ "Middle10",
               idx <= 30 ~ "Last10",
               TRUE ~ NA_character_
             )) %>%
      filter(!is.na(bucket), exp == 0) %>%
      group_by(id, gender, color, bucket) %>%
      summarise(val = m(infl), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = bucket, values_from = val)
    
    core %>% left_join(buckets, by = c("id","gender","color"))
  }
  
  # Men/Women rows + per-color t-test row
  summarise_section_with_t <- function(d) {
    if (nrow(d) == 0) return(tibble())
    pid <- per_id_stats(d)
    
    grp <- pid %>%
      group_by(color, gender) %>%
      summarise(
        Earnings   = pair_cell(m(earn_sum),  sd_(earn_sum),  money = TRUE),
        Explosions = pair_cell(m(exp_sum),   sd_(exp_sum),   money = FALSE),
        Total      = pair_cell(m(adj_total), sd_(adj_total), money = FALSE),
        `First 10` = pair_cell(m(First10),   sd_(First10),   money = FALSE),
        `Middle 10`= pair_cell(m(Middle10),  sd_(Middle10),  money = FALSE),
        `Last 10`  = pair_cell(m(Last10),    sd_(Last10),    money = FALSE),
        .groups = "drop"
      ) %>%
      arrange(factor(color, levels = c("Blue","Orange","Yellow","Pink")),
              factor(gender, levels = c("Men","Women")))
    
    tt_rows <- pid %>%
      group_by(color) %>%
      group_modify(~{
        msk_m <- .x$gender == "Men"
        msk_w <- .x$gender == "Women"
        tibble(
          gender = "t-test (M vs W)",
          Earnings   = tt_compact(.x$earn_sum[msk_m],  .x$earn_sum[msk_w]),
          Explosions = tt_compact(.x$exp_sum[msk_m],   .x$exp_sum[msk_w]),
          Total      = tt_compact(.x$adj_total[msk_m], .x$adj_total[msk_w]),
          `First 10` = tt_compact(.x$First10[msk_m],   .x$First10[msk_w]),
          `Middle 10`= tt_compact(.x$Middle10[msk_m],  .x$Middle10[msk_w]),
          `Last 10`  = tt_compact(.x$Last10[msk_m],    .x$Last10[msk_w])
        )
      }) %>% ungroup()
    
    bind_rows(
      grp %>% filter(gender == "Men"),
      grp %>% filter(gender == "Women"),
      tt_rows
    ) %>%
      arrange(factor(color, levels = c("Blue","Orange","Yellow","Pink")),
              match(gender, c("Men","Women","t-test (M vs W)"))) %>%
      mutate(`Dependent measure` = paste(color, gender)) %>%
      select(`Dependent measure`, Earnings, Explosions, Total, `First 10`, `Middle 10`, `Last 10`)
  }
  
  # Collapsed over gender
  summarise_section_collapsed <- function(d) {
    if (nrow(d) == 0) return(tibble())
    pid <- per_id_stats(d)
    pid %>%
      group_by(color) %>%
      summarise(
        Earnings   = pair_cell(m(earn_sum),  sd_(earn_sum),  money = TRUE),
        Explosions = pair_cell(m(exp_sum),   sd_(exp_sum),   money = FALSE),
        Total      = pair_cell(m(adj_total), sd_(adj_total), money = FALSE),
        `First 10` = pair_cell(m(First10),   sd_(First10),   money = FALSE),
        `Middle 10`= pair_cell(m(Middle10),  sd_(Middle10),  money = FALSE),
        `Last 10`  = pair_cell(m(Last10),    sd_(Last10),    money = FALSE),
        .groups = "drop"
      ) %>%
      arrange(factor(color, levels = c("Blue","Orange","Yellow","Pink"))) %>%
      mutate(`Dependent measure` = color) %>%
      select(`Dependent measure`, Earnings, Explosions, Total, `First 10`, `Middle 10`, `Last 10`)
  }
  
  # ── Sections ───────────────────────────────────────────────────────────────
  sec1_pre  <- dat %>% filter(sess == 1, trial >= 1,  trial <= 90,  color %in% c("Blue","Orange","Yellow"))
  sec1_post <- dat %>% filter(sess == 1, trial >= 91, trial <= 180, color %in% c("Blue","Orange","Yellow"))
  sec2      <- dat %>% filter(sess == 2, color == "Pink")
  
  # By-gender + t-tests
  tab_pre   <- summarise_section_with_t(sec1_pre)  %>% mutate(Section = "Session 1 — Pre-reversal")
  tab_post  <- summarise_section_with_t(sec1_post) %>% mutate(Section = "Session 1 — Post-reversal")
  tab_s2    <- summarise_section_with_t(sec2)      %>% mutate(Section = "Session 2")
  final_bygender <- bind_rows(tab_pre, tab_post, tab_s2) %>%
    relocate(Section, .before = `Dependent measure`)
  
  # Collapsed
  c_pre   <- summarise_section_collapsed(sec1_pre)  %>% mutate(Section = "Session 1 — Pre-reversal")
  c_post  <- summarise_section_collapsed(sec1_post) %>% mutate(Section = "Session 1 — Post-reversal")
  c_sess2 <- summarise_section_collapsed(sec2)      %>% mutate(Section = "Session 2")
  final_collapsed <- bind_rows(c_pre, c_post, c_sess2) %>%
    relocate(Section, .before = `Dependent measure`)
  
  # ── Export ─────────────────────────────────────────────────────────────────
  # By-gender file
  f1 <- paste0(out_prefix, "_bygender_ttests.xlsx")
  wb1 <- createWorkbook()
  addWorksheet(wb1, "By gender + t-tests")
  writeData(wb1, "By gender + t-tests", final_bygender)
  hdr <- createStyle(textDecoration = "bold", halign = "center", border = "Bottom")
  addStyle(wb1, "By gender + t-tests", hdr, rows = 1, cols = 1:ncol(final_bygender), gridExpand = TRUE)
  left   <- createStyle(halign = "left")
  center <- createStyle(halign = "center")
  addStyle(wb1, "By gender + t-tests", left,   rows = 2:(nrow(final_bygender)+1), cols = 1:2, gridExpand = TRUE)
  addStyle(wb1, "By gender + t-tests", center, rows = 2:(nrow(final_bygender)+1), cols = 3:ncol(final_bygender), gridExpand = TRUE)
  setColWidths(wb1, "By gender + t-tests", cols = 1, widths = 28)
  setColWidths(wb1, "By gender + t-tests", cols = 2, widths = 24)
  setColWidths(wb1, "By gender + t-tests", cols = 3:ncol(final_bygender), widths = 18)
  saveWorkbook(wb1, f1, overwrite = TRUE)
  
  # Collapsed file
  f2 <- paste0(out_prefix, "_collapsed.xlsx")
  wb2 <- createWorkbook()
  addWorksheet(wb2, "Collapsed (no gender)")
  writeData(wb2, "Collapsed (no gender)", final_collapsed)
  addStyle(wb2, "Collapsed (no gender)", hdr, rows = 1, cols = 1:ncol(final_collapsed), gridExpand = TRUE)
  addStyle(wb2, "Collapsed (no gender)", left,   rows = 2:(nrow(final_collapsed)+1), cols = 1:2, gridExpand = TRUE)
  addStyle(wb2, "Collapsed (no gender)", center, rows = 2:(nrow(final_collapsed)+1), cols = 3:ncol(final_collapsed), gridExpand = TRUE)
  setColWidths(wb2, "Collapsed (no gender)", cols = 1, widths = 28)
  setColWidths(wb2, "Collapsed (no gender)", cols = 2, widths = 22)
  setColWidths(wb2, "Collapsed (no gender)", cols = 3:ncol(final_collapsed), widths = 18)
  saveWorkbook(wb2, f2, overwrite = TRUE)
  
  message("Wrote: ", f1)
  message("Wrote: ", f2)
  
  invisible(list(by_gender = final_bygender, collapsed = final_collapsed))
}



