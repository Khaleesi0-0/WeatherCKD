suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(splines)
  library(mvmeta)
})

ckd_file <- "data/processed/ckd_2018_2022_cleaned.csv"
weather_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"

cluster_features_out <- "results/tables/n18_state_cluster_features.csv"
cluster_assignments_out <- "results/tables/n18_state_cluster_assignments.csv"
state_model_out <- "results/tables/n18_state_state_model_coefficients.csv"
cluster_meta_out <- "results/tables/n18_state_cluster_period_effects.csv"
cluster_curve_out <- "results/tables/n18_state_cluster_period_curves.csv"
method_notes_out <- "results/tables/n18_state_cluster_method_notes.txt"

n_clusters <- 4
absolute_reference_c <- 15.6
absolute_heat_c <- 26.7
absolute_cold_c <- 4.4
period_breaks <- c(2018L, 2022L, 2026L)
period_labels <- c("2018-2021", "2022-2025")

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

make_basis_spec <- function(x, probs = c(0.10, 0.50, 0.90), degree = 2) {
  x <- x[is.finite(x)]
  rng <- range(x, na.rm = TRUE)
  knots <- as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  knots <- unique(knots)

  if (length(knots) < length(probs)) {
    knots <- seq(rng[1], rng[2], length.out = length(probs) + 2)[2:(length(probs) + 1)]
  }

  list(
    knots = knots,
    boundary = rng,
    degree = degree
  )
}

build_basis_matrix <- function(x, spec, prefix) {
  mat <- splines::bs(
    x,
    degree = spec$degree,
    knots = spec$knots,
    Boundary.knots = spec$boundary,
    intercept = FALSE
  )
  colnames(mat) <- paste0(prefix, seq_len(ncol(mat)))
  as.data.frame(mat)
}

predict_component_rr <- function(temp_value, ref_value, beta, vcov_mat, spec) {
  if (!is.finite(temp_value) || !is.finite(ref_value)) {
    return(c(rr = NA_real_, rr_low_95 = NA_real_, rr_high_95 = NA_real_))
  }

  lower <- spec$boundary[1]
  upper <- spec$boundary[2]
  if (temp_value < lower || temp_value > upper || ref_value < lower || ref_value > upper) {
    return(c(rr = NA_real_, rr_low_95 = NA_real_, rr_high_95 = NA_real_))
  }

  basis_temp <- as.numeric(build_basis_matrix(temp_value, spec, "tmp"))
  basis_ref <- as.numeric(build_basis_matrix(ref_value, spec, "tmp"))
  diff_vec <- basis_temp - basis_ref
  eta <- sum(diff_vec * beta)
  se <- sqrt(drop(diff_vec %*% vcov_mat %*% diff_vec))

  c(
    rr = exp(eta),
    rr_low_95 = exp(eta - 1.96 * se),
    rr_high_95 = exp(eta + 1.96 * se)
  )
}

safe_qtest_value <- function(q_obj, field) {
  if (is.null(q_obj)) {
    return(NA_real_)
  }
  if (!field %in% names(q_obj)) {
    return(NA_real_)
  }
  as.numeric(q_obj[[field]][1])
}

ckd <- read.csv(ckd_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    year = as.integer(year),
    deaths = suppressWarnings(as.numeric(deaths))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(year), !is.na(deaths), deaths >= 0)

weather <- read.csv(weather_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    tmean_C_weekly = suppressWarnings(as.numeric(tmean_C_weekly))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(tmean_C_weekly))

analysis_df <- ckd %>%
  group_by(state, week_end_date, year) %>%
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
  arrange(state, week_end_date) %>%
  group_by(state) %>%
  mutate(
    temp_lag0 = tmean_C_weekly,
    temp_lag1 = lag(tmean_C_weekly, 1),
    temp_lag2 = lag(tmean_C_weekly, 2),
    temp_lag3 = lag(tmean_C_weekly, 3),
    temp_lag4 = lag(tmean_C_weekly, 4),
    temp_lag5 = lag(tmean_C_weekly, 5),
    temp_ma15 = ifelse(
      rowSums(is.na(cbind(temp_lag1, temp_lag2, temp_lag3, temp_lag4, temp_lag5))) == 0,
      rowMeans(cbind(temp_lag1, temp_lag2, temp_lag3, temp_lag4, temp_lag5)),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(period))

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
hc <- hclust(dist(scaled_features), method = "single")
cluster_assignments <- tibble(
  state = cluster_features$state,
  cluster = paste0("cluster_", cutree(hc, k = n_clusters))
)

write.csv(cluster_features, cluster_features_out, row.names = FALSE)
write.csv(cluster_assignments, cluster_assignments_out, row.names = FALSE)

analysis_df <- analysis_df %>%
  inner_join(cluster_assignments, by = "state")

fit_state_model <- function(df_state, lag0_spec, ma15_spec) {
  model_df <- df_state %>%
    filter(!is.na(temp_lag0), !is.na(temp_ma15)) %>%
    arrange(week_end_date) %>%
    mutate(time_num = as.numeric(week_end_date))

  if (nrow(model_df) < 52) {
    return(NULL)
  }

  n_years <- length(unique(model_df$year))
  time_df <- max(4, min(6 * n_years, floor(nrow(model_df) / 8)))

  lag0_basis <- build_basis_matrix(model_df$temp_lag0, lag0_spec, "lag0_b")
  ma15_basis <- build_basis_matrix(model_df$temp_ma15, ma15_spec, "ma15_b")
  basis_names <- c(colnames(lag0_basis), colnames(ma15_basis))

  model_frame <- bind_cols(model_df, lag0_basis, ma15_basis)
  formula_txt <- paste(
    "deaths ~ ns(time_num, df = ", time_df, ") + ",
    paste(basis_names, collapse = " + "),
    sep = ""
  )

  fit <- glm(
    as.formula(formula_txt),
    family = quasipoisson(),
    data = model_frame,
    na.action = na.exclude
  )

  list(
    state = unique(model_df$state)[1],
    n_obs = nrow(model_df),
    n_years = n_years,
    coef = coef(fit)[basis_names],
    vcov = vcov(fit)[basis_names, basis_names, drop = FALSE]
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

  lag0_spec <- make_basis_spec(subset_df$temp_lag0)
  ma15_spec <- make_basis_spec(subset_df$temp_ma15[is.finite(subset_df$temp_ma15)])

  state_fits <- lapply(split(subset_df, subset_df$state), fit_state_model, lag0_spec = lag0_spec, ma15_spec = ma15_spec)
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
      term = names(state_fit$coef),
      coefficient = as.numeric(state_fit$coef),
      stringsAsFactors = FALSE
    )
    state_counter <- state_counter + 1L
  }

  coef_matrix <- do.call(rbind, lapply(state_fits, function(x) x$coef))
  vcov_list <- lapply(state_fits, function(x) x$vcov)
  meta_fit <- mvmeta::mvmeta(coef_matrix, vcov_list, method = "reml")

  pooled_beta <- coef(meta_fit)
  pooled_vcov <- vcov(meta_fit)
  q_info <- tryCatch(mvmeta::qtest(meta_fit), error = function(e) NULL)

  lag0_names <- grep("^lag0_b", names(pooled_beta), value = TRUE)
  ma15_names <- grep("^ma15_b", names(pooled_beta), value = TRUE)

  lag0_beta <- pooled_beta[lag0_names]
  lag0_vcov <- pooled_vcov[lag0_names, lag0_names, drop = FALSE]
  ma15_beta <- pooled_beta[ma15_names]
  ma15_vcov <- pooled_vcov[ma15_names, ma15_names, drop = FALSE]

  heat_abs <- predict_component_rr(absolute_heat_c, absolute_reference_c, lag0_beta, lag0_vcov, lag0_spec)
  cold_abs <- predict_component_rr(absolute_cold_c, absolute_reference_c, ma15_beta, ma15_vcov, ma15_spec)

  lag0_pct <- as.numeric(stats::quantile(subset_df$temp_lag0, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE))
  ma15_pct <- as.numeric(stats::quantile(subset_df$temp_ma15, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE))
  heat_rel <- predict_component_rr(lag0_pct[3], lag0_pct[2], lag0_beta, lag0_vcov, lag0_spec)
  cold_rel <- predict_component_rr(ma15_pct[1], ma15_pct[2], ma15_beta, ma15_vcov, ma15_spec)

  cluster_meta_rows[[meta_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    n_states = length(state_fits),
    heat_rr_26_7_vs_15_6 = heat_abs["rr"],
    heat_rr_low_95 = heat_abs["rr_low_95"],
    heat_rr_high_95 = heat_abs["rr_high_95"],
    cold_rr_4_4_vs_15_6 = cold_abs["rr"],
    cold_rr_low_95 = cold_abs["rr_low_95"],
    cold_rr_high_95 = cold_abs["rr_high_95"],
    relative_heat_rr_p99_vs_p50 = heat_rel["rr"],
    relative_heat_rr_low_95 = heat_rel["rr_low_95"],
    relative_heat_rr_high_95 = heat_rel["rr_high_95"],
    relative_cold_rr_p01_vs_p50 = cold_rel["rr"],
    relative_cold_rr_low_95 = cold_rel["rr_low_95"],
    relative_cold_rr_high_95 = cold_rel["rr_high_95"],
    q_stat = safe_qtest_value(q_info, "qstat"),
    q_df = safe_qtest_value(q_info, "df"),
    q_pvalue = safe_qtest_value(q_info, "pvalue"),
    i2 = safe_qtest_value(q_info, "i2"),
    stringsAsFactors = FALSE
  )
  meta_counter <- meta_counter + 1L

  lag0_grid <- seq(lag0_spec$boundary[1], lag0_spec$boundary[2], length.out = 100)
  ma15_grid <- seq(ma15_spec$boundary[1], ma15_spec$boundary[2], length.out = 100)

  lag0_curve <- t(sapply(lag0_grid, function(x) {
    predict_component_rr(x, absolute_reference_c, lag0_beta, lag0_vcov, lag0_spec)
  }))
  ma15_curve <- t(sapply(ma15_grid, function(x) {
    predict_component_rr(x, absolute_reference_c, ma15_beta, ma15_vcov, ma15_spec)
  }))

  cluster_curve_rows[[curve_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    component = "lag0",
    temperature_c = lag0_grid,
    rr = lag0_curve[, "rr"],
    rr_low_95 = lag0_curve[, "rr_low_95"],
    rr_high_95 = lag0_curve[, "rr_high_95"],
    stringsAsFactors = FALSE
  )
  curve_counter <- curve_counter + 1L

  cluster_curve_rows[[curve_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    component = "ma15",
    temperature_c = ma15_grid,
    rr = ma15_curve[, "rr"],
    rr_low_95 = ma15_curve[, "rr_low_95"],
    rr_high_95 = ma15_curve[, "rr_high_95"],
    stringsAsFactors = FALSE
  )
  curve_counter <- curve_counter + 1L
}

if (length(state_model_rows) == 0 || length(cluster_meta_rows) == 0) {
  stop("No cluster-period combinations had enough state-level data to fit the adapted GAM/meta-analysis workflow.")
}

write.csv(bind_rows(state_model_rows), state_model_out, row.names = FALSE)
write.csv(bind_rows(cluster_meta_rows), cluster_meta_out, row.names = FALSE)
write.csv(bind_rows(cluster_curve_rows), cluster_curve_out, row.names = FALSE)

writeLines(
  c(
    "Adapted N18 state-level cluster/meta-smoothing workflow",
    "",
    "Inputs:",
    paste0("- ", ckd_file),
    paste0("- ", weather_file),
    "",
    "Method notes:",
    "- Uses states as the study units in place of cities.",
    "- Uses weekly state mortality counts and weekly mean temperature because daily state-level N18 files are not present in this repo.",
    "- Uses hierarchical single-linkage clustering on state seasonal temperature mean and SD.",
    "- Uses two multi-year periods based on available state-week data: 2018-2021 and 2022-2025.",
    "- Uses quasi-Poisson GAMs with a natural spline for time plus quadratic B-spline terms for current-week temperature and lag 1-5 week mean temperature.",
    "- Cluster-specific temperature basis knots are shared within each cluster-period to preserve coefficient interpretation across states inside that cluster-period.",
    "- This is not a literal replication of the published city-daily pipeline; it is the closest state-week adaptation supported by the local files."
  ),
  method_notes_out
)

cat("Wrote:", cluster_features_out, "\n")
cat("Wrote:", cluster_assignments_out, "\n")
cat("Wrote:", state_model_out, "\n")
cat("Wrote:", cluster_meta_out, "\n")
cat("Wrote:", cluster_curve_out, "\n")
cat("Wrote:", method_notes_out, "\n")
