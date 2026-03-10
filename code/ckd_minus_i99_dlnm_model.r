suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(dlnm)
  library(gnm)
})

ckd_clean_file <- "data/processed/ckd_2018_2022_cleaned.csv"
i99_clean_file <- "data/processed/i99_n00n99_cleaned.csv"
weather_weekly_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"

out_clean_diff <- "data/processed/ckd_minus_i99_cleaned.csv"
out_model_input <- "data/processed/ckd_minus_i99_with_weather_weekly.csv"
out_summary <- "results/tables/ckd_minus_i99_dlnm_model_summary.txt"
out_curve <- "results/tables/ckd_minus_i99_dlnm_cumulative_rr_curve.csv"
out_pct <- "results/tables/ckd_minus_i99_dlnm_cumulative_rr_percentiles.csv"
out_plot <- "results/figures/ckd_minus_i99_dlnm_cumulative_rr_plot.png"

max_lag_weeks <- 3

ckd <- read.csv(ckd_clean_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(deaths), deaths >= 0) %>%
  group_by(state, week_end_date, mmwr_week_code) %>%
  summarise(ckd_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

i99 <- read.csv(i99_clean_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(deaths), deaths >= 0) %>%
  group_by(state, week_end_date, mmwr_week_code) %>%
  summarise(i99_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

# Keep overlapping state-week rows only, then compute difference.
diff_clean <- inner_join(
  ckd,
  i99,
  by = c("state", "week_end_date", "mmwr_week_code")
) %>%
  mutate(
    deaths_diff = ckd_deaths - i99_deaths,
    year = suppressWarnings(as.integer(sub("/.*", "", mmwr_week_code))),
    week = suppressWarnings(as.integer(sub(".*/", "", mmwr_week_code)))
  ) %>%
  filter(!is.na(deaths_diff), deaths_diff >= 0) %>%
  arrange(state, week_end_date)

write.csv(diff_clean, out_clean_diff, row.names = FALSE, na = "")

weather_weekly <- read.csv(weather_weekly_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    tmean_C_weekly = suppressWarnings(as.numeric(tmean_C_weekly)),
    n_days_weather = suppressWarnings(as.numeric(n_days_weather))
  ) %>%
  group_by(state, week_end_date) %>%
  summarise(
    tmean_C_weekly = first(tmean_C_weekly),
    n_days_weather = first(n_days_weather),
    .groups = "drop"
  )

model_input <- diff_clean %>%
  left_join(weather_weekly, by = c("state", "week_end_date")) %>%
  filter(!is.na(tmean_C_weekly))

write.csv(model_input, out_model_input, row.names = FALSE, na = "")

# Same model approach as prior script.
model_df <- model_input %>%
  mutate(
    deaths = as.numeric(deaths_diff),
    week_end_date = as.Date(week_end_date)
  ) %>%
  filter(!is.na(deaths), !is.na(tmean_C_weekly), deaths >= 0)

temp_bounds <- quantile(model_df$tmean_C_weekly, probs = c(0.01, 0.99), na.rm = TRUE)
model_df <- model_df %>%
  filter(tmean_C_weekly >= temp_bounds[[1]], tmean_C_weekly <= temp_bounds[[2]]) %>%
  mutate(
    year = year(week_end_date),
    month = month(week_end_date),
    stratum = as.factor(interaction(state, year, month, drop = TRUE))
  )

valid_strata <- model_df %>%
  count(stratum, name = "n_obs") %>%
  filter(n_obs >= 2) %>%
  pull(stratum)
model_df <- model_df %>% filter(stratum %in% valid_strata)

var_knots <- quantile(model_df$tmean_C_weekly, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
lag_knots <- logknots(max_lag_weeks, nk = 2)

cb_temp <- crossbasis(
  model_df$tmean_C_weekly,
  lag = max_lag_weeks,
  argvar = list(fun = "ns", knots = var_knots),
  arglag = list(fun = "ns", knots = lag_knots)
)

fit <- gnm(
  deaths ~ cb_temp,
  eliminate = stratum,
  family = quasipoisson(),
  data = model_df,
  na.action = na.exclude
)

pred0 <- crosspred(cb_temp, fit, by = 0.1)
mmt <- pred0$predvar[which.min(pred0$allRRfit)]
pred <- crosspred(cb_temp, fit, cen = mmt, by = 0.1)

curve_df <- data.frame(
  temperature_C = pred$predvar,
  rr_cumulative = pred$allRRfit,
  rr_low_95 = pred$allRRlow,
  rr_high_95 = pred$allRRhigh
)

pct_vals <- quantile(model_df$tmean_C_weekly, probs = c(0.01, 0.10, 0.50, 0.90, 0.99), na.rm = TRUE)
pct_df <- data.frame(
  percentile = c("p01", "p10", "p50", "p90", "p99"),
  temperature_C = as.numeric(pct_vals)
)
nearest_idx <- sapply(pct_df$temperature_C, function(x) which.min(abs(curve_df$temperature_C - x)))
pct_out <- data.frame(
  percentile = pct_df$percentile,
  temperature_C = pct_df$temperature_C,
  rr_cumulative = curve_df$rr_cumulative[nearest_idx],
  rr_low_95 = curve_df$rr_low_95[nearest_idx],
  rr_high_95 = curve_df$rr_high_95[nearest_idx]
)

write.csv(curve_df, out_curve, row.names = FALSE)
write.csv(pct_out, out_pct, row.names = FALSE)

png(out_plot, width = 1200, height = 760, res = 120)
plot(
  curve_df$temperature_C, curve_df$rr_cumulative,
  type = "l", lwd = 2, col = "#6a1b9a",
  xlab = "Weekly Mean Temperature (C)",
  ylab = "Cumulative Relative Risk (RR)",
  main = "DLNM: (CKD deaths - I99 N00-N99 deaths) vs Temperature"
)
lines(curve_df$temperature_C, curve_df$rr_low_95, col = "#8e24aa", lty = 2)
lines(curve_df$temperature_C, curve_df$rr_high_95, col = "#8e24aa", lty = 2)
abline(h = 1, lty = 3, col = "gray40")
abline(v = mmt, lty = 3, col = "#37474f")
grid(col = "gray85")
dev.off()

sink(out_summary)
cat("DLNM model summary: CKD minus I99 deaths\n")
cat("Matched overlap rows (before weather join):", nrow(diff_clean), "\n")
cat("Rows used in model:", nrow(model_df), "\n")
cat("States used:", length(unique(model_df$state)), "\n")
cat("Temperature restriction (1st-99th percentile): ",
    round(temp_bounds[[1]], 3), " to ", round(temp_bounds[[2]], 3), " C\n", sep = "")
cat("Max lag (weeks):", max_lag_weeks, "\n")
cat("Centering temperature (MMT):", round(mmt, 3), "C\n\n")
print(summary(fit))
cat("\nCumulative RR at selected percentiles\n")
print(pct_out)
sink()

cat("Done. Wrote:\n")
cat("-", out_clean_diff, "\n")
cat("-", out_model_input, "\n")
cat("-", out_summary, "\n")
cat("-", out_curve, "\n")
cat("-", out_pct, "\n")
cat("-", out_plot, "\n")
