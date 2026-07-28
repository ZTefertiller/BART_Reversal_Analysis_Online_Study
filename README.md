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
| `modeling/` | Stan models (STL, 4PAR, EWMV) and the scripts that produced the fits (`fit_models.R`, `build_model_fits_list.R`, shared Stan-data builder `stan_data.R`) plus parameter recovery for all three models (`run_recovery.R`, the per-model simulators `stl/simulate_stl.R`, `fourpar/simulate_fourpar.R`, `ewmv/simulate_ewmv.R`, and the report `parameter_recovery.Rmd`). Fit outputs stay local (`mcmc/`, gitignored). |

## Parameter recovery

```r
Rscript modeling/run_recovery.R            # all 9 model x condition cells
Rscript modeling/run_recovery.R ewmv       # one model, all three conditions
Rscript modeling/run_recovery.R stl blue   # a single cell
```

True parameters are drawn from each cell's population distribution, data are
simulated from that model's generative process on the condition's real fixed
breakpoints, and the model is refit. Results go to
`mcmc/recovery/recovery_<model>_<condition>.csv`; render
`modeling/parameter_recovery.Rmd` to see the grids and correlation tables.

This runs on a clean clone: with no fit objects present it moment-matches the
population to the per-participant posterior means in `data_published/`
(`RECOVERY_TRUTH=published`). Those means are shrunk toward the group, so the
implied population SD is slightly narrow and recovery correlations are, if
anything, conservative. With the fit objects from OSF in `mcmc/`, it uses the
fitted group mean/SD instead (`RECOVERY_TRUTH=fit`); the default `auto` picks
whichever is available.

## Render

```r
rmarkdown::render("reporting/BART_RL_Online_Analysis.Rmd")
```

Chunks that need the local fits (`mcmc/model_fits_list.rds`) are the only parts
a data-only clone cannot reproduce; every behavioural and trait-correlation
figure renders from `data_published/` alone.

Raw data and data processing steps are removed to preserve anonymity of participants. 
