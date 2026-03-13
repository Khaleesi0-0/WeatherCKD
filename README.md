# WeatherCKD

## Layout

- `code/`: active R scripts.
- `data/raw/`: source CSVs and the downloaded weather archive.
- `data/processed/`: cleaned and merged analysis-ready datasets.
- `results/tables/`: model summaries and tabular outputs.
- `results/figures/`: plots and exported figures.

## Removed

The unused two-stage DLNM workflow and its outputs were removed:

- `code/ckd_two_stage_dlnm.R`
- all `two_stage_*` files

## N18 State Method

- `code/n18_state_cluster_meta_gam.R` retains the state-level clustering step for the N18 CKD query and now calculates weather impact with a cluster-period DLNM workflow using conditional Poisson estimation with case-crossover-style strata, restricting temperature to the 1st-99th percentile and reporting pooled cumulative relative-risk-scale contrasts with 95 % confidence intervals.
- The adaptation uses states instead of cities and weekly temperatures/counts instead of daily series, because those are the inputs available in this repo.
