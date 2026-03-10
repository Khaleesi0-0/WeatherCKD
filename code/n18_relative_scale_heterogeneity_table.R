suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(splines)
})

ckd_file <- "data/processed/ckd_2018_2022_cleaned.csv"
weather_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"
cluster_assignments_file <- "results/tables/n18_state_cluster_assignments.csv"
cluster_features_file <- "results/tables/n18_state_cluster_features.csv"

state_effects_out <- "results/tables/n18_state_relative_scale_state_effects.csv"
heterogeneity_out <- "results/tables/n18_state_relative_scale_heterogeneity_summary.csv"

period_breaks <- c(2018L, 2022L, 2026L)
period_labels <- c("2018-2021", "2022-2025")

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

make_pos_def <- function(mat, min_eigen = 1e-6) {
  if (any(!is.finite(mat))) {
    return(NULL)
  }

  mat <- (mat + t(mat)) / 2
  eig <- eigen(mat, symmetric = TRUE)
  vals <- eig$values
  if (any(!is.finite(vals))) {
    return(NULL)
  }

  vals[vals < min_eigen] <- min_eigen
  repaired <- eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
  repaired <- (repaired + t(repaired)) / 2
  colnames(repaired) <- colnames(mat)
  rownames(repaired) <- rownames(mat)
  repaired
}

contrast_stats <- function(temp_value, ref_value, beta, vcov_mat, spec) {
  lower <- spec$boundary[1]
  upper <- spec$boundary[2]

  if (!is.finite(temp_value) || !is.finite(ref_value) ||
      temp_value < lower || temp_value > upper ||
      ref_value < lower || ref_value > upper) {
    return(c(log_rr = NA_real_, se = NA_real_, rr = NA_real_, low = NA_real_, high = NA_real_))
  }

  basis_temp <- as.numeric(build_basis_matrix(temp_value, spec, "tmp"))
  basis_ref <- as.numeric(build_basis_matrix(ref_value, spec, "tmp"))
  diff_vec <- basis_temp - basis_ref
  log_rr <- sum(diff_vec * beta)
  se <- sqrt(drop(diff_vec %*% vcov_mat %*% diff_vec))

  c(
    log_rr = log_rr,
    se = se,
    rr = exp(log_rr),
    low = exp(log_rr - 1.96 * se),
    high = exp(log_rr + 1.96 * se)
  )
}

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
    coef = coef_vec,
    vcov = vcov_mat,
    lag0_values = model_df$temp_lag0,
    lag13_values = model_df$temp_lag13
  )
}

fit_heterogeneity_model <- function(df, scale_label, model_label, formula_rhs) {
  model_df <- df %>%
    filter(is.finite(log_rr), is.finite(var_log_rr), var_log_rr > 0)

  if (nrow(model_df) < 3) {
    return(NULL)
  }

  X <- model.matrix(as.formula(paste("~", formula_rhs)), data = model_df)
  y <- model_df$log_rr
  w <- 1 / model_df$var_log_rr

  XtW <- t(X) %*% (w * X)
  beta <- solve(XtW, t(X) %*% (w * y))
  resid <- y - drop(X %*% beta)
  q_stat <- sum(w * resid^2)
  df_q <- nrow(model_df) - ncol(X)
  p_value <- if (df_q > 0) pchisq(q_stat, df = df_q, lower.tail = FALSE) else NA_real_
  i2 <- if (is.finite(q_stat) && q_stat > 0 && df_q > 0) max(0, 100 * (q_stat - df_q) / q_stat) else NA_real_

  data.frame(
    scale = scale_label,
    model = model_label,
    q = q_stat,
    df = df_q,
    p = p_value,
    i2_percent = i2,
    n_effects = nrow(model_df),
    stringsAsFactors = FALSE
  )
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

cluster_assignments <- read.csv(cluster_assignments_file, stringsAsFactors = FALSE)
cluster_features <- read.csv(cluster_features_file, stringsAsFactors = FALSE)

analysis_df <- ckd %>%
  group_by(state, week_end_date, year) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    weather %>% select(state, week_end_date, tmean_C_weekly),
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
  filter(!is.na(period)) %>%
  inner_join(cluster_assignments, by = "state")

state_effect_rows <- list()
row_idx <- 1L

cluster_period_keys <- analysis_df %>%
  distinct(cluster, period) %>%
  arrange(cluster, period)

for (i in seq_len(nrow(cluster_period_keys))) {
  cluster_i <- cluster_period_keys$cluster[i]
  period_i <- cluster_period_keys$period[i]

  subset_df <- analysis_df %>%
    filter(cluster == cluster_i, period == period_i)

  lag0_spec <- make_basis_spec(subset_df$temp_lag0)
  lag13_spec <- make_basis_spec(subset_df$temp_lag13[is.finite(subset_df$temp_lag13)])

  state_fits <- lapply(split(subset_df, subset_df$state), fit_state_model, lag0_spec = lag0_spec, lag13_spec = lag13_spec)
  state_fits <- Filter(Negate(is.null), state_fits)

  if (!length(state_fits)) {
    next
  }

  for (state_fit in state_fits) {
    lag0_names <- grep("^lag0_b", names(state_fit$coef), value = TRUE)
    lag13_names <- grep("^lag13_b", names(state_fit$coef), value = TRUE)

    lag0_beta <- state_fit$coef[lag0_names]
    lag0_vcov <- state_fit$vcov[lag0_names, lag0_names, drop = FALSE]
    lag13_beta <- state_fit$coef[lag13_names]
    lag13_vcov <- state_fit$vcov[lag13_names, lag13_names, drop = FALSE]

    lag0_pct <- as.numeric(stats::quantile(state_fit$lag0_values, probs = c(0.50, 0.99), na.rm = TRUE, names = FALSE))
    lag13_pct <- as.numeric(stats::quantile(state_fit$lag13_values, probs = c(0.01, 0.50), na.rm = TRUE, names = FALSE))

    heat_stats <- contrast_stats(lag0_pct[2], lag0_pct[1], lag0_beta, lag0_vcov, lag0_spec)
    cold_stats <- contrast_stats(lag13_pct[1], lag13_pct[2], lag13_beta, lag13_vcov, lag13_spec)

    state_effect_rows[[row_idx]] <- data.frame(
      state = state_fit$state,
      cluster = cluster_i,
      period = period_i,
      scale = "Heat effect",
      contrast = "Relative scale",
      log_rr = heat_stats["log_rr"],
      se = heat_stats["se"],
      var_log_rr = heat_stats["se"]^2,
      rr = heat_stats["rr"],
      rr_low_95 = heat_stats["low"],
      rr_high_95 = heat_stats["high"],
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L

    state_effect_rows[[row_idx]] <- data.frame(
      state = state_fit$state,
      cluster = cluster_i,
      period = period_i,
      scale = "Cold effect",
      contrast = "Relative scale",
      log_rr = cold_stats["log_rr"],
      se = cold_stats["se"],
      var_log_rr = cold_stats["se"]^2,
      rr = cold_stats["rr"],
      rr_low_95 = cold_stats["low"],
      rr_high_95 = cold_stats["high"],
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L
  }
}

state_effects <- bind_rows(state_effect_rows) %>%
  left_join(cluster_features, by = "state")

write.csv(state_effects, state_effects_out, row.names = FALSE)

heat_df <- state_effects %>% filter(scale == "Heat effect")
cold_df <- state_effects %>% filter(scale == "Cold effect")

heterogeneity_rows <- bind_rows(
  fit_heterogeneity_model(heat_df, "Relative scale", "Intercept only", "1"),
  fit_heterogeneity_model(heat_df, "Relative scale", "Clusters, time", "cluster + period"),
  fit_heterogeneity_model(heat_df, "Relative scale", "Clusters, time, summer mean", "cluster + period + seasonal_mean_temp_summer"),
  fit_heterogeneity_model(cold_df, "Relative scale", "Intercept only", "1"),
  fit_heterogeneity_model(cold_df, "Relative scale", "Clusters, time", "cluster + period"),
  fit_heterogeneity_model(cold_df, "Relative scale", "Clusters, time, winter mean", "cluster + period + seasonal_mean_temp_winter")
) %>%
  mutate(
    effect = c("Heat effect", "Heat effect", "Heat effect", "Cold effect", "Cold effect", "Cold effect")
  ) %>%
  select(effect, scale, model, q, df, p, i2_percent, n_effects)

write.csv(heterogeneity_rows, heterogeneity_out, row.names = FALSE)

cat("Wrote:", state_effects_out, "\n")
cat("Wrote:", heterogeneity_out, "\n")
