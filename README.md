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

- `code/n18_state_cluster_meta_gam.R` applies a state-level adaptation of the cluster-plus-meta GAM approach for the N18 CKD query already reflected in the weekly state files.
- The adaptation uses states instead of cities and weekly temperatures/counts instead of daily series, because those are the inputs available in this repo.
