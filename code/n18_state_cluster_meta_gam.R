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

n_clusters_target <- 4
min_cluster_size <- 3
clustering_method <- "ward.D2"
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

make_pos_def <- function(mat, min_eigen = 1e-6) {
  if (any(!is.finite(mat))) {
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
    temp_lag13 = ifelse(
      rowSums(is.na(cbind(temp_lag1, temp_lag2, temp_lag3))) == 0,
      rowMeans(cbind(temp_lag1, temp_lag2, temp_lag3)),
      NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(period))

absolute_reference_c <- as.numeric(stats::quantile(analysis_df$tmean_C_weekly, probs = 0.50, na.rm = TRUE, names = FALSE))
state_percentile_summary <- analysis_df %>%
  group_by(state) %>%
  summarise(
    temp_p10 = as.numeric(stats::quantile(tmean_C_weekly, probs = 0.10, na.rm = TRUE, names = FALSE)),
    temp_p90 = as.numeric(stats::quantile(tmean_C_weekly, probs = 0.90, na.rm = TRUE, names = FALSE)),
    .groups = "drop"
  )
absolute_cold_c <- mean(state_percentile_summary$temp_p10, na.rm = TRUE)
absolute_heat_c <- mean(state_percentile_summary$temp_p90, na.rm = TRUE)
reference_label <- gsub("\\.", "_", format(round(absolute_reference_c, 1), nsmall = 1))
heat_label <- gsub("\\.", "_", format(round(absolute_heat_c, 1), nsmall = 1))
cold_label <- gsub("\\.", "_", format(round(absolute_cold_c, 1), nsmall = 1))

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
    n_clusters <- n_clusters - 1
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

fit_state_model <- function(df_state, lag0_spec, lag13_spec) {
  model_df <- df_state %>%
    filter(!is.na(temp_lag0), !is.na(temp_lag13)) %>%
    arrange(week_end_date) %>%
    mutate(time_num = as.numeric(week_end_date))

  if (nrow(model_df) < 52) {
    return(NULL)
  }

  n_years <- length(unique(model_df$year))
  time_df <- max(4, min(6 * n_years, floor(nrow(model_df) / 8)))

  lag0_basis <- build_basis_matrix(model_df$temp_lag0, lag0_spec, "lag0_b")
  lag13_basis <- build_basis_matrix(model_df$temp_lag13, lag13_spec, "lag13_b")
  basis_names <- c(colnames(lag0_basis), colnames(lag13_basis))

  model_frame <- bind_cols(model_df, lag0_basis, lag13_basis)
  formula_txt <- paste(
    "deaths ~ ns(time_num, df = ", time_df, ") + ",
    paste(basis_names, collapse = " + "),
    sep = ""
  )

  fit <- tryCatch(
    glm(
      as.formula(formula_txt),
      family = quasipoisson(),
      data = model_frame,
      na.action = na.exclude
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(NULL)
  }

  coef_vec <- coef(fit)[basis_names]
  vcov_mat <- vcov(fit)[basis_names, basis_names, drop = FALSE]
  vcov_mat <- make_pos_def(vcov_mat)

  if (is.null(vcov_mat) || any(!is.finite(coef_vec))) {
    return(NULL)
  }

  list(
    state = unique(model_df$state)[1],
    n_obs = nrow(model_df),
    n_years = n_years,
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

  lag0_spec <- make_basis_spec(subset_df$temp_lag0)
  lag13_spec <- make_basis_spec(subset_df$temp_lag13[is.finite(subset_df$temp_lag13)])

  state_fits <- lapply(split(subset_df, subset_df$state), fit_state_model, lag0_spec = lag0_spec, lag13_spec = lag13_spec)
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

  pooled_beta <- coef(meta_fit)
  pooled_vcov <- vcov(meta_fit)
  q_info <- tryCatch(mvmeta::qtest(meta_fit), error = function(e) NULL)

  lag0_names <- grep("^lag0_b", names(pooled_beta), value = TRUE)
  lag13_names <- grep("^lag13_b", names(pooled_beta), value = TRUE)

  lag0_beta <- pooled_beta[lag0_names]
  lag0_vcov <- pooled_vcov[lag0_names, lag0_names, drop = FALSE]
  lag13_beta <- pooled_beta[lag13_names]
  lag13_vcov <- pooled_vcov[lag13_names, lag13_names, drop = FALSE]

  heat_abs <- predict_component_rr(absolute_heat_c, absolute_reference_c, lag0_beta, lag0_vcov, lag0_spec)
  cold_abs <- predict_component_rr(absolute_cold_c, absolute_reference_c, lag13_beta, lag13_vcov, lag13_spec)

  lag0_pct <- as.numeric(stats::quantile(subset_df$temp_lag0, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE))
  lag13_pct <- as.numeric(stats::quantile(subset_df$temp_lag13, probs = c(0.01, 0.50, 0.99), na.rm = TRUE, names = FALSE))
  heat_rel <- predict_component_rr(lag0_pct[3], lag0_pct[2], lag0_beta, lag0_vcov, lag0_spec)
  cold_rel <- predict_component_rr(lag13_pct[1], lag13_pct[2], lag13_beta, lag13_vcov, lag13_spec)

  cluster_meta_rows[[meta_counter]] <- data.frame(
    cluster = cluster_i,
    period = period_i,
    n_states = length(state_fits),
    heat_rr = heat_abs["rr"],
    heat_rr_low_95 = heat_abs["rr_low_95"],
    heat_rr_high_95 = heat_abs["rr_high_95"],
    cold_rr = cold_abs["rr"],
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
  names(cluster_meta_rows[[meta_counter]])[names(cluster_meta_rows[[meta_counter]]) == "heat_rr"] <- paste0("heat_rr_", heat_label, "_vs_", reference_label)
  names(cluster_meta_rows[[meta_counter]])[names(cluster_meta_rows[[meta_counter]]) == "cold_rr"] <- paste0("cold_rr_", cold_label, "_vs_", reference_label)
  meta_counter <- meta_counter + 1L

  lag0_grid <- seq(lag0_spec$boundary[1], lag0_spec$boundary[2], length.out = 100)
  lag13_grid <- seq(lag13_spec$boundary[1], lag13_spec$boundary[2], length.out = 100)

  lag0_curve <- t(sapply(lag0_grid, function(x) {
    predict_component_rr(x, absolute_reference_c, lag0_beta, lag0_vcov, lag0_spec)
  }))
  lag13_curve <- t(sapply(lag13_grid, function(x) {
    predict_component_rr(x, absolute_reference_c, lag13_beta, lag13_vcov, lag13_spec)
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
    component = "lag13",
    temperature_c = lag13_grid,
    rr = lag13_curve[, "rr"],
    rr_low_95 = lag13_curve[, "rr_low_95"],
    rr_high_95 = lag13_curve[, "rr_high_95"],
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
    paste0("- Uses hierarchical clustering with method ", clustering_method, " on state seasonal temperature mean and SD."),
    paste0("- Targets ", n_clusters_target, " clusters with a minimum cluster size of ", min_cluster_size, "; if needed, the script reduces the cluster count until that requirement is met."),
    "- Uses two multi-year periods based on available state-week data: 2018-2021 and 2022-2025.",
    "- Uses quasi-Poisson GAMs with a natural spline for time plus quadratic B-spline terms for current-week temperature and the lag 1-3 week mean temperature.",
    paste0("- Uses the empirical median weekly temperature as the centering temperature (MMT/reference): ", round(absolute_reference_c, 3), " C."),
    paste0("- Defines the absolute cold contrast as the average state-specific 10th percentile (", round(absolute_cold_c, 3), " C) versus ", round(absolute_reference_c, 3), " C."),
    paste0("- Defines the absolute heat contrast as the average state-specific 90th percentile (", round(absolute_heat_c, 3), " C) versus ", round(absolute_reference_c, 3), " C."),
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
