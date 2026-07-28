# BART_RL_Method



(N = 118). 

Main analysis script to regenerate results from the paper is:
`reporting/BART_RL_Online_Analysis.Rmd`

## Layout

| Path | What it is |
|---|---|
| `data_published/` | De-identified CSVs — everything needed to replicate. `bart_rl_online_trials.csv` (24,780 trials × 118 participants), `bart_rl_online_summary.csv` (one row per participant), `bart_rl_online_ewmv_params.csv` (per-participant EWMV parameter means per task phase). Public key: `sub_id` (1–118). No `participant_id` codes, no raw data. |
| `reporting/BART_RL_Online_Analysis.Rmd` | The single analysis report. Reads only `data_published/` for data; model-fit objects (large `.rds`, not in the repo - find on OSF) are needed only for the LOOIC table, parameter recovery, and Bayesian group contrasts. |
| `analysis/` | The plotting/table functions the report sources — nothing else. |
| `modeling/` | Stan models (STL, 4PAR, EWMV) and the scripts that produced the fits (`fit_models.R`, `build_model_fits_list.R`) plus EWMV parameter recovery (`ewmv/run_ewmv_recovery.R`, `ewmv/simulate_ewmv.R`, `parameter_recovery.Rmd`). Fit outputs stay local (`mcmc/`, gitignored). |

## Render

```r
rmarkdown::render("reporting/BART_RL_Online_Analysis.Rmd")
```

Chunks that need the local fits (`mcmc/model_fits_list.rds`) are the only parts
a data-only clone cannot reproduce; every behavioural and trait-correlation
figure renders from `data_published/` alone.

Raw data and data processing steps are removed to preserve anonymity of participants. 
