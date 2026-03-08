# Two-stage DLNM analysis for CKD deaths and temperature
# Stage 1: state-specific time-series models
# Stage 2: random-effects meta-analysis + BLUP

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(splines)
  library(dlnm)
  library(mvmeta)
})

input_file <- "data/ckd_2018_2022_with_weather_weekly.csv"
out_dir <- "data"

# Settings that mirror the referenced framework
max_lag <- 7
lag_df <- 3
calendar_df_per_year <- 3
doy_df <- 3
trim_temp_quantiles <- c(0.01, 0.99)

# Optional confounder names (used only if present)
humidity_col <- "humidity"
pm10_col <- "PM10"
o3_col <- "O3"
holiday_col <- "PH"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

required_pkgs <- c("dlnm", "mvmeta")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

dat <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths)),
    temp = suppressWarnings(as.numeric(tmean_C_weekly))
  ) %>%
  filter(!is.na(state), state != "", !is.na(date), !is.na(deaths), !is.na(temp), deaths >= 0)

# Restrict to 1st-99th percentile to reduce instability at extremes.
q <- quantile(dat$temp, probs = trim_temp_quantiles, na.rm = TRUE)
dat <- dat %>% filter(temp >= q[[1]], temp <= q[[2]])

# Time controls following the referenced design.
dat <- dat %>%
  arrange(state, date) %>%
  mutate(
    cal_num = as.numeric(date - min(date)) + 1,
    doy = yday(date),
    dow = factor(
      wday(date, week_start = 1),
      levels = 1:7,
      labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")
    ),
    year = year(date)
  )

if (!(holiday_col %in% names(dat))) {
  dat[[holiday_col]] <- 0L
}

if (humidity_col %in% names(dat)) {
  dat[[humidity_col]] <- suppressWarnings(as.numeric(dat[[humidity_col]]))
}
if (pm10_col %in% names(dat)) {
  dat[[pm10_col]] <- suppressWarnings(as.numeric(dat[[pm10_col]]))
}
if (o3_col %in% names(dat)) {
  dat[[o3_col]] <- suppressWarnings(as.numeric(dat[[o3_col]]))
}

# Build formula pieces based on available data.
base_terms <- c(
  "cb_temp",
  sprintf("ns(cal_num, df = %d)", max(3, calendar_df_per_year * length(unique(dat$year)))),
  sprintf("ns(doy, df = %d)", doy_df)
)

if (dplyr::n_distinct(dat$dow) > 1) {
  base_terms <- c(base_terms, "dow")
}
if (holiday_col %in% names(dat) && dplyr::n_distinct(dat[[holiday_col]]) > 1) {
  base_terms <- c(base_terms, holiday_col)
}

# Humidity cross-basis (only if available)
use_humidity <- humidity_col %in% names(dat) && any(!is.na(dat[[humidity_col]]))
if (use_humidity) {
  base_terms <- c(base_terms, "cb_humidity")
}

# Air pollutants as lag-matched moving averages (only if available)
if (pm10_col %in% names(dat) && any(!is.na(dat[[pm10_col]]))) {
  base_terms <- c(base_terms, "pm10_ma")
}
if (o3_col %in% names(dat) && any(!is.na(dat[[o3_col]]))) {
  base_terms <- c(base_terms, "o3_ma")
}

model_formula <- as.formula(paste("deaths ~", paste(base_terms, collapse = " + ")))

# Stage 1: fit state-level models and extract cumulative log-RR per +1 C.
states <- sort(unique(dat$state))
coef_list <- list()
vcov_list <- list()
first_stage_rows <- list()
af_rows <- list()

roll_mean_forward <- function(x, k) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    j <- min(n, i + k)
    out[i] <- mean(x[i:j], na.rm = TRUE)
  }
  out
}

for (s in states) {
  ds <- dat %>% filter(state == s) %>% arrange(date)
  if (nrow(ds) < (max_lag + 15)) next

  # State-specific pollutant moving averages with same lag span.
  if (pm10_col %in% names(ds)) {
    ds$pm10_ma <- stats::filter(ds[[pm10_col]], rep(1 / (max_lag + 1), max_lag + 1), sides = 1)
  }
  if (o3_col %in% names(ds)) {
    ds$o3_ma <- stats::filter(ds[[o3_col]], rep(1 / (max_lag + 1), max_lag + 1), sides = 1)
  }

  # Temperature cross-basis: linear exposure-response + ns lag-response.
  cb_temp <- crossbasis(
    ds$temp,
    lag = max_lag,
    argvar = list(fun = "lin"),
    arglag = list(fun = "ns", df = lag_df)
  )

  if (use_humidity) {
    cb_humidity <- crossbasis(
      ds[[humidity_col]],
      lag = max_lag,
      argvar = list(fun = "ns", df = 3),
      arglag = list(fun = "ns", df = 3)
    )
  }

  fit <- try(glm(model_formula, family = quasipoisson(), data = ds, na.action = na.exclude), silent = TRUE)
  if (inherits(fit, "try-error")) next

  red <- try(crossreduce(cb_temp, fit, type = "overall"), silent = TRUE)
  if (inherits(red, "try-error") || any(!is.finite(red$coef))) next

  coef_list[[s]] <- red$coef
  vcov_list[[s]] <- red$vcov

  first_stage_rows[[s]] <- data.frame(
    state = s,
    log_rr_per_1C = as.numeric(red$coef),
    rr_per_1C = exp(as.numeric(red$coef)),
    se = sqrt(diag(red$vcov))[1],
    stringsAsFactors = FALSE
  )
}

if (length(coef_list) < 3) {
  stop("Too few successful first-stage fits for meta-analysis.")
}

coef_mat <- do.call(rbind, lapply(coef_list, function(x) matrix(x, nrow = 1)))
rownames(coef_mat) <- names(coef_list)
S_list <- unname(vcov_list[names(coef_list)])

# Stage 2: random-effects meta-analysis (no predictors)
meta_fit <- mvmeta(coef_mat ~ 1, S = S_list, method = "reml")
bl <- blup(meta_fit, vcov = TRUE)

pooled_log_rr <- as.numeric(coef(meta_fit)[1])
pooled_se <- sqrt(vcov(meta_fit)[1, 1])
pooled_rr <- exp(pooled_log_rr)
pooled_ci <- exp(c(pooled_log_rr - 1.96 * pooled_se, pooled_log_rr + 1.96 * pooled_se))

first_stage_df <- bind_rows(first_stage_rows)
first_stage_df <- first_stage_df %>%
  mutate(
    rr_low_95 = exp(log_rr_per_1C - 1.96 * se),
    rr_high_95 = exp(log_rr_per_1C + 1.96 * se)
  )

if (is.list(bl) && length(bl) > 0 && is.list(bl[[1]]) && all(c("blup", "vcov") %in% names(bl[[1]]))) {
  blup_df <- data.frame(
    state = names(bl),
    blup_log_rr_per_1C = vapply(bl, function(x) as.numeric(x$blup), numeric(1)),
    blup_se = vapply(bl, function(x) sqrt(diag(x$vcov))[1], numeric(1)),
    stringsAsFactors = FALSE
  )
} else {
  if (is.matrix(bl) && ncol(bl) >= 2 && all(c("blup", "vcov") %in% colnames(bl))) {
    state_ids <- rownames(bl)
    if (is.null(state_ids)) state_ids <- rownames(coef_mat)
    blup_vals <- as.numeric(bl[, "blup"])
    blup_se <- sqrt(as.numeric(bl[, "vcov"]))
  } else {
    blup_vals <- as.numeric(bl)
    state_ids <- names(blup_vals)
    if (is.null(state_ids) || any(state_ids == "")) {
      state_ids <- rownames(coef_mat)
    }
    bl_vcov <- attr(bl, "vcov")
    if (!is.null(bl_vcov) && length(dim(bl_vcov)) == 3) {
      blup_se <- sqrt(vapply(seq_along(state_ids), function(i) bl_vcov[1, 1, i], numeric(1)))
    } else if (!is.null(bl_vcov) && is.matrix(bl_vcov)) {
      blup_se <- sqrt(diag(bl_vcov))
    } else {
      blup_se <- rep(NA_real_, length(state_ids))
    }
  }
  blup_df <- data.frame(
    state = state_ids,
    blup_log_rr_per_1C = blup_vals,
    blup_se = blup_se,
    stringsAsFactors = FALSE
  )
}

blup_df <- blup_df %>%
  mutate(
    blup_rr_per_1C = exp(blup_log_rr_per_1C),
    rr_low_95 = exp(blup_log_rr_per_1C - 1.96 * blup_se),
    rr_high_95 = exp(blup_log_rr_per_1C + 1.96 * blup_se)
  )

# Attributable cases/fraction using BLUP city(state)-specific cumulative association.
# RR_i is relative to state-specific minimum-risk temperature over observed range.
for (s in blup_df$state) {
  ds <- dat %>% filter(state == s) %>% arrange(date)
  beta_s <- blup_df %>%
    filter(state == s) %>%
    slice(1) %>%
    pull(blup_log_rr_per_1C)
  if (length(beta_s) == 0 || !is.finite(beta_s)) next

  tmin <- min(ds$temp, na.rm = TRUE)
  tmax <- max(ds$temp, na.rm = TRUE)
  mmt_s <- if (beta_s >= 0) tmin else tmax

  rr_i <- exp(beta_s * (ds$temp - mmt_s))
  Ni <- roll_mean_forward(ds$deaths, max_lag)
  ac_i <- Ni * (rr_i - 1) / rr_i

  af_rows[[s]] <- data.frame(
    state = s,
    date = ds$date,
    temp = ds$temp,
    deaths = ds$deaths,
    Ni = Ni,
    rr_i = rr_i,
    ac_i = ac_i,
    stringsAsFactors = FALSE
  )
}

af_daily <- bind_rows(af_rows)
af_state <- af_daily %>%
  group_by(state) %>%
  summarise(
    total_ac = sum(ac_i, na.rm = TRUE),
    total_deaths = sum(deaths, na.rm = TRUE),
    af = total_ac / total_deaths,
    .groups = "drop"
  )

af_national <- af_state %>%
  summarise(
    total_ac = sum(total_ac, na.rm = TRUE),
    total_deaths = sum(total_deaths, na.rm = TRUE),
    af = total_ac / total_deaths
  )

write.csv(first_stage_df, file.path(out_dir, "two_stage_first_stage_state_estimates.csv"), row.names = FALSE)
write.csv(blup_df, file.path(out_dir, "two_stage_blup_state_estimates.csv"), row.names = FALSE)
write.csv(af_daily, file.path(out_dir, "two_stage_attributable_daily.csv"), row.names = FALSE)
write.csv(af_state, file.path(out_dir, "two_stage_attributable_state.csv"), row.names = FALSE)
write.csv(af_national, file.path(out_dir, "two_stage_attributable_national.csv"), row.names = FALSE)

summary_txt <- file.path(out_dir, "two_stage_model_summary.txt")
con <- file(summary_txt, open = "wt")
writeLines(c(
  "Two-stage DLNM summary",
  paste0("Input file: ", input_file),
  paste0("Rows used after temperature trimming: ", nrow(dat)),
  paste0("Temperature trimming (1st-99th percentile): ", round(q[[1]], 3), " to ", round(q[[2]], 3), " C"),
  paste0("Successful first-stage units (states): ", nrow(first_stage_df)),
  paste0("Max lag: ", max_lag),
  paste0("Lag df: ", lag_df),
  "",
  "Pooled national cumulative association (per +1 C)",
  paste0("log(RR): ", round(pooled_log_rr, 5)),
  paste0("RR: ", round(pooled_rr, 5), " (95% CI ", round(pooled_ci[1], 5), ", ", round(pooled_ci[2], 5), ")"),
  "",
  "National attributable fraction",
  paste0("AF: ", round(af_national$af, 5)),
  paste0("Total attributable cases: ", round(af_national$total_ac, 2)),
  paste0("Total observed deaths: ", round(af_national$total_deaths, 2)),
  "",
  "Notes:",
  "1) This script follows the two-stage framework but uses state-level weekly deaths (not daily city hospitalizations).",
  "2) humidity/PM10/O3 terms are included only if those columns exist in the input data.",
  "3) DOW and PH are retained for structural consistency with the reference method."
), con)
close(con)

cat("Done. Wrote outputs in data/:\n")
cat("- two_stage_model_summary.txt\n")
cat("- two_stage_first_stage_state_estimates.csv\n")
cat("- two_stage_blup_state_estimates.csv\n")
cat("- two_stage_attributable_daily.csv\n")
cat("- two_stage_attributable_state.csv\n")
cat("- two_stage_attributable_national.csv\n")
