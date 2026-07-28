# analysis/plotting/questionnaire_meta.R
# Single source of truth for questionnaire display: label, colour family, scale
# type, and (where known) theoretical range. Used by the density plots and the
# summary-stat cards so colours/labels never overlap or drift between figures.
#
# Colour families (no cross-family overlap):
#   SPQ + SPQ-B + 3 factors + 9 subscales  -> cool blues / indigos / purples / cyans
#   CAPS total + 3 subscales               -> warm corals / oranges
#   PDI  total + 3 subscales               -> magentas / pinks
#   MDQ green, PHQ-9 gold, IPIP slate, PPGM brown
#
# tmax = NA  -> density normalises that variable to its OBSERVED range instead of
# a theoretical one (used for subscales whose theoretical max is not well defined).

Q_INFO <- tibble::tribble(
  ~var,                ~label,          ~family, ~color,     ~scale,        ~tmin, ~tmax,
  "spq_total",         "SPQ",           "SPQ",   "#2630F5",  "Likert 0–4", 0,   296,
  "spq_b_total",       "SPQ-B",         "SPQ",   "#0B1E7A",  "Binary",          0,   22,
  "spq_cogper",        "SPQ Cog-Per",   "SPQ",   "#3A0CA3",  "Likert 0–4", 0,   132,
  "spq_interp",        "SPQ Interp",    "SPQ",   "#7209B7",  "Likert 0–4", 0,   132,
  "spq_disorg",        "SPQ Disorg",    "SPQ",   "#4CC9F0",  "Likert 0–4", 0,   64,
  "spq_ideas",         "SPQ Ideas",     "SPQ",   "#4895EF",  "Likert 0–4", 0,   36,
  "spq_magic",         "SPQ Magic",     "SPQ",   "#4361EE",  "Likert 0–4", 0,   28,
  "spq_perception",    "SPQ Percept.",  "SPQ",   "#5179E0",  "Likert 0–4", 0,   36,
  "spq_suspicion",     "SPQ Suspic.",   "SPQ",   "#3F37C9",  "Likert 0–4", 0,   32,
  "spq_social",        "SPQ Social",    "SPQ",   "#560BAD",  "Likert 0–4", 0,   32,
  "spq_friends",       "SPQ Friends",   "SPQ",   "#6A0DAD",  "Likert 0–4", 0,   36,
  "spq_affect",        "SPQ Affect",    "SPQ",   "#8338EC",  "Likert 0–4", 0,   32,
  "spq_behavior",      "SPQ Behavior",  "SPQ",   "#48BFE3",  "Likert 0–4", 0,   28,
  "spq_speech",        "SPQ Speech",    "SPQ",   "#56CFE1",  "Likert 0–4", 0,   36,
  "caps_total",        "CAPS",          "CAPS",  "#E76F51",  "Binary",          0,   32,
  "caps_distress",     "CAPS Distress", "CAPS",  "#C44536",  "Likert 1–5", 0,   NA,
  "caps_intrusiveness","CAPS Intrus.",  "CAPS",  "#F4A261",  "Likert 1–5", 0,   NA,
  "caps_frequency",    "CAPS Freq.",    "CAPS",  "#E9967A",  "Likert 1–5", 0,   NA,
  "pdi_total",         "PDI",           "PDI",   "#FF006E",  "Binary",          0,   21,
  "pdi_distress",      "PDI Distress",  "PDI",   "#C9184A",  "Likert 1–5", 0,   NA,
  "pdi_frequency",     "PDI Freq.",     "PDI",   "#FF5C8A",  "Likert 1–5", 0,   NA,
  "pdi_conviction",    "PDI Convict.",  "PDI",   "#A4133C",  "Likert 1–5", 0,   NA,
  "mdq_total",         "MDQ",           "MDQ",   "#06D6A0",  "Binary",          0,   13,
  "phq_total",         "PHQ-9",         "PHQ",   "#FFB703",  "Likert 0–3", 0,   27,
  "ipip_total",        "IPIP",          "IPIP",  "#6C757D",  "Likert 1–5", 11,  55,
  "ppgm_total",        "PPGM",          "PPGM",  "#8C5E2A",  "Mixed",           0,   NA
)

# label -> colour, for scale_*_manual
Q_PAL <- setNames(Q_INFO$color, Q_INFO$label)

# convenient var groupings
Q_GROUPS <- list(
  paper_core   = c("spq_total", "spq_cogper", "spq_interp", "spq_disorg",
                   "caps_total", "pdi_total"),
  all_totals   = c("spq_total", "caps_total", "pdi_total", "mdq_total",
                   "phq_total"),
  totals_plus  = c("spq_total", "spq_b_total", "caps_total", "pdi_total",
                   "mdq_total", "phq_total", "ipip_total", "ppgm_total"),
  spq_all      = c("spq_total", "spq_b_total", "spq_cogper", "spq_interp",
                   "spq_disorg", "spq_ideas", "spq_magic", "spq_perception",
                   "spq_suspicion", "spq_social", "spq_friends", "spq_affect",
                   "spq_behavior", "spq_speech"),
  caps_pdi     = c("caps_total", "caps_distress", "caps_intrusiveness",
                   "caps_frequency", "pdi_total", "pdi_distress",
                   "pdi_frequency", "pdi_conviction"),
  everything   = c("spq_total", "spq_cogper", "spq_interp", "spq_disorg",
                   "spq_ideas", "spq_magic", "spq_perception", "spq_suspicion",
                   "spq_social", "spq_friends", "spq_affect", "spq_behavior",
                   "spq_speech", "caps_total", "caps_distress",
                   "caps_intrusiveness", "caps_frequency", "pdi_total",
                   "pdi_distress", "pdi_frequency", "pdi_conviction",
                   "mdq_total", "phq_total")
)

# SPQ factors (drawn dashed in density plots)
SPQ_FACTOR_LABELS_DENS <- c("SPQ Cog-Per", "SPQ Interp", "SPQ Disorg")
