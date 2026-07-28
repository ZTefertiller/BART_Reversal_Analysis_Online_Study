# Dataset and model fits

**N = 118** participants (online Prolific sample).

## Published data (`data_published/`)
De-identified; `sub_id` (1–118) is the only participant key. All files are
ordered by `sub_id` (trial-level files by `sub_id`, `session`, `trial_number`).
- `bart_rl_online_trials.csv` — all trials (sub_id, session, trial_number,
  balloon_color, behavioural columns, questionnaire totals).
- `bart_rl_online_summary.csv` — one row per participant.
- `bart_rl_online_ewmv_params.csv` — per-participant EWMV parameter means per
  task phase (drives the trait-correlation figures).

## Model fits
- 9 fits = 3 task phases (blue = pre-reversal, orange = post-reversal,
  pink = control) × 3 models (STL, 4PAR, EWMV).
- 4 chains × 1000 warmup / 2000 sampling; subjects indexed by `sub_id`.
- Fit with `modeling/fit_models.R`, assembled into `mcmc/model_fits_list.rds`
  by `modeling/build_model_fits_list.R`. The fit objects are large and stay out
  of the repo; the draws are available on OSF.

## Parameter recovery
- `modeling/ewmv/run_ewmv_recovery.R` simulates from each condition's fitted
  population and refits, writing `mcmc/ewmv_recovery/recovery_<condition>.csv`.

## Model comparison
- `output/model_comparison/` — LOOIC across models (EWMV preferred in all three
  conditions) and the group-level MCMC diagnostics table.
