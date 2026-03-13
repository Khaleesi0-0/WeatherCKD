suppressPackageStartupMessages({
  library(dplyr)
  library(dlnm)
  library(gnm)
  library(lubridate)
  library(mvmeta)
  library(splines)
  library(tidyr)
})

ckd_file <- "data/processed/ckd_2018_2022_cleaned.csv"
weather_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"

cluster_features_out <- "results/tables/n18_state_cluster_features.csv"
cluster_assignments_out <- "results/tables/n18_state_cluster_assignments.csv"
state_model_out <- "results/tables/n18_state_state_model_coefficients.csv"
cluster_meta_out <- "results/tables/n18_state_cluster_period_effects.csv"
cluster_curve_out <- "results/tables/n18_state_cluster_period_curves.csv"
method_notes_out <- "results/tables/n18_state_cluster_method_notes.txt"

n_clusters_target <- 4
min_cluster_size <- 3
clustering_method <- "ward.D2"
period_breaks <- c(2018L, 2022L, 2026L)
period_labels <- c("2018-2021", "2022-2025")
max_lag_weeks <- 3
var_knots_probs <- c(0.10, 0.75, 0.90)
temp_restriction_probs <- c(0.01, 0.99)

if (!file.exists(ckd_file)) {
  stop("Input file not found: ", ckd_file)
}
if (!file.exists(weather_file)) {
  stop("Input file not found: ", weather_file)
}

dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

season_from_month <- function(month_num) {
  dplyr::case_when(
    month_num %in% c(12L, 1L, 2L) ~ "winter",
    month_num %in% c(3L, 4L, 5L) ~ "spring",
    month_num %in% c(6L, 7L, 8L) ~ "summer",
    TRUE ~ "fall"
  )
}

period_from_year <- function(year_num) {
  cut(
    year_num,
    breaks = period_breaks,
    labels = period_labels,
    right = FALSE,
    include.lowest = TRUE
  ) %>%
    as.character()
}

safe_qtest_value <- function(q_obj, field) {
  if (is.null(q_obj) || !field %in% names(q_obj)) {
    return(NA_real_)
  }
  as.numeric(q_obj[[field]][1])
}

make_pos_def <- function(mat, min_eigen = 1e-6) {
  if (is.null(mat) || any(!is.finite(mat))) {
    return(NULL)
  }

  mat <- (mat + t(mat)) / 2
  eig <- eigen(mat, symmetric = TRUE)
  values <- eig$values

  if (any(!is.finite(values))) {
    return(NULL)
  }

  values[values < min_eigen] <- min_eigen
  repaired <- eig$vectors %*% diag(values, nrow = length(values)) %*% t(eig$vectors)
  repaired <- (repaired + t(repaired)) / 2
  colnames(repaired) <- colnames(mat)
  rownames(repaired) <- rownames(mat)
  repaired
}

build_var_basis <- function(x, argvar) {
  basis <- dlnm::onebasis(
    x,
    fun = argvar$fun,
    knots = argvar$knots,
    Boundary.knots = argvar$Boundary.knots
  )
  as.matrix(basis)
}

make_argvar <- function(x, probs) {
  bounds <- as.numeric(stats::quantile(x, probs = temp_restriction_probs, na.rm = TRUE, names = FALSE))
  knots <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  knots <- unique(knots)

  if (length(knots) == 0) {
    knots <- mean(bounds)
  } else if (length(knots) == 1) {
    knots <- knots[1]
  }

  list(
    fun = "ns",
    knots = knots,
    Boundary.knots = bounds
  )
}

predict_cumulative_rr <- function(temp_value, ref_value, beta, vcov_mat, argvar) {
  if (!is.finite(temp_value) || !is.finite(ref_value)) {
    return(c(rr = NA_real_, rr_low_95 = NA_real_, rr_high_95 = NA_real_))
  }

  bounds <- argvar$Boundary.knots
  if (temp_value < bounds[1] || temp_value > bounds[2] || ref_value < bounds[1] || ref_value > bounds[2]) {
    return(c(rr = NA_real_, rr_low_95 = NA_real_, rr_high_95 = NA_real_))
  }

  basis_temp <- as.numeric(build_var_basis(temp_value, argvar))
  basis_ref <- as.numeric(build_var_basis(ref_value, argvar))
  diff_vec <- basis_temp - basis_ref
  eta <- sum(diff_vec * beta)
  se <- sqrt(drop(diff_vec %*% vcov_mat %*% diff_vec))

  c(
    rr = exp(eta),
    rr_low_95 = exp(eta - 1.96 * se),
    rr_high_95 = exp(eta + 1.96 * se)
  )
}

ckd <- read.csv(ckd_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    year = as.integer(year),
    week = suppressWarnings(as.integer(week)),
    deaths = suppressWarnings(as.numeric(deaths))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(year), !is.na(week), !is.na(deaths), deaths >= 0)

weather <- read.csv(weather_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    tmean_C_weekly = suppressWarnings(as.numeric(tmean_C_weekly))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(tmean_C_weekly))

analysis_df <- ckd %>%
  group_by(state, week_end_date, year, week) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    weather %>%
      select(state, week_end_date, tmean_C_weekly),
    by = c("state", "week_end_date")
  ) %>%
  filter(!is.na(tmean_C_weekly)) %>%
  mutate(
    month = month(week_end_date),
    season = season_from_month(month),
    period = period_from_year(year)
  ) %>%
  filter(!is.na(period)) %>%
  arrange(state, week_end_date)

cluster_features <- analysis_df %>%
  group_by(state, season) %>%
  summarise(
    seasonal_mean_temp = mean(tmean_C_weekly, na.rm = TRUE),
    seasonal_sd_temp = sd(tmean_C_weekly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = season,
    values_from = c(seasonal_mean_temp, seasonal_sd_temp)
  ) %>%
  arrange(state)

feature_matrix <- cluster_features %>%
  select(-state) %>%
  as.matrix()
rownames(feature_matrix) <- cluster_features$state

scaled_features <- scale(feature_matrix)
hc <- hclust(dist(scaled_features), method = clustering_method)
cluster_sizes_ok <- FALSE
n_clusters <- n_clusters_target

while (n_clusters >= 2 && !cluster_sizes_ok) {
  cluster_ids <- cutree(hc, k = n_clusters)
  cluster_sizes_ok <- min(as.integer(table(cluster_ids))) >= min_cluster_size
  if (!cluster_sizes_ok) {
    n_clusters <- n_clusters - 1L
  }
}

if (!cluster_sizes_ok) {
  stop("Unable to construct clusters meeting the minimum size requirement.")
}

cluster_assignments <- tibble(
  state = cluster_features$state,
  cluster = paste0("cluster_", cutree(hc, k = n_clusters))
)

write.csv(cluster_features, cluster_features_out, row.names = FALSE)
write.csv(cluster_assignments, cluster_assignments_out, row.names = FALSE)

analysis_df <- analysis_df %>%
  inner_join(cluster_assignments, by = "state")

fit_state_model <- function(df_state, argvar, arglag, temp_bounds) {
  model_df <- df_state %>%
    arrange(week_end_date) %>%
    filter(
      is.finite(tmean_C_weekly),
      tmean_C_weekly >= temp_bounds[1],
      tmean_C_weekly <= temp_bounds[2]
    ) %>%
    mutate(
      stratum = as.factor(interaction(state, year, month, drop = TRUE))
    )

  valid_strata <- model_df %>%
    count(stratum, name = "n_obs") %>%
    filter(n_obs >= 2) %>%
    pull(stratum)

  model_df <- model_df %>% filter(stratum %in% valid_strata)

  if (nrow(model_df) < 52 || length(unique(model_df$stratum)) < 10) {
    return(NULL)
  }

  cb_temp <- dlnm::crossbasis(
    model_df$tmean_C_weekly,
    lag = max_lag_weeks,
    argvar = argvar,
    arglag = arglag
  )

  fit <- tryCatch(
    gnm::gnm(
      deaths ~ cb_temp,
      eliminate = stratum,
      family = quasipoisson(),
      data = model_df,
      na.action = na.exclude
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  pred0 <- tryCatch(
    dlnm::crosspred(cb_temp, fit, by = 0.1),
    error = function(e) NULL
  )

  if (is.null(pred0)) {
    return(NULL)
  }

  centering_temp <- pred0$predvar[which.min(pred0$allRRfit)]
  red <- tryCatch(
    dlnm::crossreduce(cb_temp, fit, type = "overall", cen = centering_temp),
    error = function(e) NULL
  )

  if (is.null(red)) {
    return(NULL)
  }

  coef_vec <- stats::coef(red)
  vcov_mat <- make_pos_def(stats::vcov(red))

  if (is.null(vcov_mat) || any(!is.finite(coef_vec))) {
    return(NULL)
  }

  list(
    state = unique(model_df$state)[1],
    n_obs = nrow(model_df),
    n_years = length(unique(model_df$year)),
    n_strata = length(unique(model_df$stratum)),
    centering_temp = centering_temp,
    coef = coef_vec,
    vcov = vcov_mat
  )
}

state_model_rows <- list()
cluster_meta_rows <- list()
cluster_curve_rows <- list()
state_counter <- 1L
meta_counter <- 1L
curve_counter <- 1L

cluster_period_keys <- analysis_df %>%
  distinct(cluster, period) %>%
  arrange(cluster, period)

for (i in seq_len(nrow(cluster_period_keys))) {
  cluster_i <- cluster_period_keys$cluster[i]
  period_i <- cluster_period_keys$period[i]

  subset_df <- analysis_df %>%
    filter(cluster == cluster_i, period == period_i) %>%
    arrange(state, week_end_date)

  state_counts <- subset_df %>%
    distinct(state) %>%
    nrow()

  if (state_counts < 2) {
    next
  }

  temp_bounds <- as.numeric(stats::quantile(
    subset_df$tmean_C_weekly,
    probs = temp_restriction_probs,
    na.rm = TRUE,
    names = FALSE
  ))

  restricted_df <- subset_df %>%
    filter(
      is.finite(tmean_C_weekly),
      tmean_C_weekly >= temp_bounds[1],
      tmean_C_weekly <= temp_bounds[2]
    )

  if (nrow(restricted_df) < 100) {
    next
  }

  argvar <- make_argvar(restricted_df$tmean_C_weekly, var_knots_probs)
  arglag <- list(
    fun = "ns",
    knots = dlnm::logknots(max_lag_weeks, nk = 2)
  )

  state_fits <- lapply(
    split(subset_df, subset_df$state),
    fit_state_model,
    argvar = argvar,
    arglag = arglag,
    temp_bounds = temp_bounds
  )
  state_fits <- Filter(Negate(is.null), state_fits)

  if (length(state_fits) < 2) {
    next
  }

  for (state_fit in state_fits) {
    state_model_rows[[state_counter]] <- data.frame(
      cluster = cluster_i,
      period = period_i,
      state = state_fit$state,
      n_obs = state_fit$n_obs,
      n_years = state_fit$n_years,
      n_strata = state_fit$n_strata,
      centering_temp_c = state_fit$centering_temp,
      term = names(state_fit$coef),
      coefficient = as.numeric(state_fit$coef),
      stringsAsFactors = FALSE
    )
    state_counter <- state_counter + 1L
  }

  coef_matrix <- do.call(rbind, lapply(state_fits, function(x) x$coef))
  vcov_list <- lapply(state_fits, function(x) x$vcov)
  meta_fit <- tryCatch(
    mvmeta::mvmeta(coef_matrix, vcov_list, method = "reml", bscov = "diag"),
    error = function(e) {
      tryCatch(
        mvmeta::mvmeta(coef_matrix, vcov_list, method = "fixed"),
        error = function(e2) NULL
      )
    }
  )

  if (is.null(meta_fit)) {
    next
  }

  pooled_beta <- stats::coef(meta_fit)
  pooled_vcov <- make_pos_def(stats::vcov(meta_fit))
  if (is.null(pooled_vcov)) {
    next
  }
  q_info <- tryCatch(mvmeta::qtest(meta_fit), error = function(e) NULL)

  state_centering_values <- vapply(state_fits, function(x) x$centering_temp, numeric(1))
  pooled_centering_temp <- mean(state_centering_values, na.rm = TRUE)

  temp_percentiles <- as.numeric(stats::quantile(
    restricted_df$tmean_C_weekly,
    probs = c(0.01, 0.50, 0.99),
    na.rm = TRUE,
    names = FALSE
  ))

  heat_effect <- predict_cumulative_rr(
    temp_percentiles[3],
    pooled_centering_temp,
    pooled_beta,
    pooled_vcov,
    argvar
  )
  cold_effect <- predict_cumulative_rr(
    temp_percentiles[1],
    pooled_centering_temp,
    pooled_beta,
    pooled_vcov,
    argvar
  )
  relative_heat <- predict_cumulative_rr(
    temp_percentiles[3],
    temp_percentiles[2],
    pooled_beta,
    pooled_vcov,
    argvar
  )
  relative_cold <- predict_cumulative_rr(
    temp_percentiles[1],
    temp_percentiles[2],
    pooled_beta,
    pooled_vcov,
    argvar
  )

  cluster_meta_rows[[meta_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    n_states = length(state_fits),
    temp_p01_c = temp_percentiles[1],
    temp_p50_c = temp_percentiles[2],
    temp_p99_c = temp_percentiles[3],
    centering_temp_c = pooled_centering_temp,
    heat_rr_p99_vs_center = heat_effect["rr"],
    heat_rr_low_95 = heat_effect["rr_low_95"],
    heat_rr_high_95 = heat_effect["rr_high_95"],
    cold_rr_p01_vs_center = cold_effect["rr"],
    cold_rr_low_95 = cold_effect["rr_low_95"],
    cold_rr_high_95 = cold_effect["rr_high_95"],
    relative_heat_rr_p99_vs_p50 = relative_heat["rr"],
    relative_heat_rr_low_95 = relative_heat["rr_low_95"],
    relative_heat_rr_high_95 = relative_heat["rr_high_95"],
    relative_cold_rr_p01_vs_p50 = relative_cold["rr"],
    relative_cold_rr_low_95 = relative_cold["rr_low_95"],
    relative_cold_rr_high_95 = relative_cold["rr_high_95"],
    q_stat = safe_qtest_value(q_info, "qstat"),
    q_df = safe_qtest_value(q_info, "df"),
    q_pvalue = safe_qtest_value(q_info, "pvalue"),
    i2 = safe_qtest_value(q_info, "i2"),
    stringsAsFactors = FALSE
  )
  meta_counter <- meta_counter + 1L

  temp_grid <- seq(temp_bounds[1], temp_bounds[2], length.out = 100)
  curve_vals <- t(vapply(temp_grid, function(x) {
    predict_cumulative_rr(
      x,
      pooled_centering_temp,
      pooled_beta,
      pooled_vcov,
      argvar
    )
  }, numeric(3)))

  cluster_curve_rows[[curve_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    temperature_c = temp_grid,
    rr_cumulative = curve_vals[, "rr"],
    rr_low_95 = curve_vals[, "rr_low_95"],
    rr_high_95 = curve_vals[, "rr_high_95"],
    centering_temp_c = pooled_centering_temp,
    temp_lower_bound_c = temp_bounds[1],
    temp_upper_bound_c = temp_bounds[2],
    stringsAsFactors = FALSE
  )
  curve_counter <- curve_counter + 1L
}

if (length(state_model_rows) == 0 || length(cluster_meta_rows) == 0) {
  stop("No cluster-period combinations had enough state-level data to fit the cluster DLNM/meta-analysis workflow.")
}

write.csv(bind_rows(state_model_rows), state_model_out, row.names = FALSE)
write.csv(bind_rows(cluster_meta_rows), cluster_meta_out, row.names = FALSE)
write.csv(bind_rows(cluster_curve_rows), cluster_curve_out, row.names = FALSE)

writeLines(
  c(
    "Adapted N18 state-level cluster/meta-DLNM workflow",
    "",
    "Inputs:",
    paste0("- ", ckd_file),
    paste0("- ", weather_file),
    "",
    "Method notes:",
    "- Uses states as the study units in place of cities.",
    "- Uses weekly state mortality counts and weekly mean temperature because daily state-level N18 files are not present in this repo.",
    paste0("- Uses hierarchical clustering with method ", clustering_method, " on state seasonal temperature mean and SD."),
    paste0("- Targets ", n_clusters_target, " clusters with a minimum cluster size of ", min_cluster_size, "; if needed, the script reduces the cluster count until that requirement is met."),
    "- Uses two multi-year periods based on available state-week data: 2018-2021 and 2022-2025.",
    "- Within each cluster-period, the weather-impact model uses a DLNM with conditional Poisson estimation via stratum elimination, matching the case-crossover stratification structure used for conditional logistic analyses on aggregated weekly counts.",
    "- Temperature is restricted to the 1st-99th percentile within each cluster-period before model fitting and effect estimation.",
    paste0("- Maximum lag is ", max_lag_weeks, " weeks, with a natural spline temperature basis and log-spaced lag knots."),
    "- State-specific cumulative temperature associations are reduced to the overall temperature basis and then pooled within each cluster-period using multivariate meta-analysis.",
    "- Effect estimates are reported as cumulative relative-risk-scale contrasts with 95 % confidence intervals for the pooled cluster-period curves.",
    "- This is not a literal replication of the published city-daily conditional logistic pipeline; it is the closest state-week adaptation supported by the local files while keeping the clustering method unchanged."
  ),
  method_notes_out
)

cat("Wrote:", cluster_features_out, "\n")
cat("Wrote:", cluster_assignments_out, "\n")
cat("Wrote:", state_model_out, "\n")
cat("Wrote:", cluster_meta_out, "\n")
cat("Wrote:", cluster_curve_out, "\n")
cat("Wrote:", method_notes_out, "\n")
