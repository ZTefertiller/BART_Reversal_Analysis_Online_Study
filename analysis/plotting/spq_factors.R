# analysis/plotting/spq_factors.R
# Add SPQ three-factor scores (Raine 1994 classic model) as columns.
#
# Mapping (Suspiciousness loads on BOTH Cognitive-Perceptual and Interpersonal,
# per Raine et al. 1994):
#   Cognitive-Perceptual = ideas of reference + magical thinking +
#                          unusual perceptual experiences + suspiciousness
#   Interpersonal        = social anxiety + no close friends +
#                          constricted affect + suspiciousness
#   Disorganized         = odd/eccentric behaviour + odd speech
#
# Subscale columns expected (Likert-summed): spq_ideas, spq_magic,
# spq_perception, spq_suspicion, spq_social, spq_friends, spq_affect,
# spq_behavior, spq_speech.

# Display labels for the three factors (HTML/markdown-friendly).
SPQ_FACTOR_LABELS <- c(
  spq_cogper  = "SPQ Cognitive-Perceptual",
  spq_interp  = "SPQ Interpersonal",
  spq_disorg  = "SPQ Disorganized"
)

add_spq_factors <- function(df) {
  needed <- c("spq_ideas", "spq_magic", "spq_perception", "spq_suspicion",
              "spq_social", "spq_friends", "spq_affect",
              "spq_behavior", "spq_speech")
  miss <- setdiff(needed, names(df))
  if (length(miss)) {
    stop("add_spq_factors: missing SPQ subscale columns: ",
         paste(miss, collapse = ", "))
  }

  # na.rm = FALSE: if any component is missing the factor is NA, so that
  # participant is dropped from downstream correlations (consistent with totals).
  df$spq_cogper <- rowSums(
    df[, c("spq_ideas", "spq_magic", "spq_perception", "spq_suspicion")],
    na.rm = FALSE
  )
  df$spq_interp <- rowSums(
    df[, c("spq_social", "spq_friends", "spq_affect", "spq_suspicion")],
    na.rm = FALSE
  )
  df$spq_disorg <- rowSums(
    df[, c("spq_behavior", "spq_speech")],
    na.rm = FALSE
  )
  df
}
