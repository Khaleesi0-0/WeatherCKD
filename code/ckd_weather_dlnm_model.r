library(dplyr)
library(lubridate)
library(dlnm)
library(splines)
library(gnm)

input_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"
out_summary <- "results/tables/dlnm_model_summary.txt"
out_curve <- "results/tables/dlnm_cumulative_rr_curve.csv"
out_pct <- "results/tables/dlnm_cumulative_rr_percentiles.csv"
out_plot <- "results/figures/dlnm_cumulative_rr_plot.png"

max_lag_weeks <- 3

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths)),
    tmean_C_weekly = suppressWarnings(as.numeric(tmean_C_weekly))
  ) %>%
  filter(!is.na(week_end_date), !is.na(deaths), !is.na(tmean_C_weekly), deaths >= 0)

# Restrict exposure to 1st-99th percentile, per requested method.
temp_bounds <- quantile(df$tmean_C_weekly, probs = c(0.01, 0.99), na.rm = TRUE)
df <- df %>%
  filter(tmean_C_weekly >= temp_bounds[[1]], tmean_C_weekly <= temp_bounds[[2]])

# Time-stratified matching set: same state, year, and month.
df <- df %>%
  mutate(
    year = year(week_end_date),
    month = month(week_end_date),
    stratum = as.factor(interaction(state, year, month, drop = TRUE))
  )

# Keep strata with at least 2 observations so within-stratum contrasts exist.
valid_strata <- df %>%
  count(stratum, name = "n_obs") %>%
  filter(n_obs >= 2) %>%
  pull(stratum)

df <- df %>% filter(stratum %in% valid_strata)

if (nrow(df) == 0) {
  stop("No rows available after filtering; cannot fit model.")
}

var_knots <- quantile(df$tmean_C_weekly, probs = c(0.10, 0.75, 0.90), na.rm = TRUE)
lag_knots <- logknots(max_lag_weeks, nk = 2)

cb_temp <- crossbasis(
  df$tmean_C_weekly,
  lag = max_lag_weeks,
  argvar = list(fun = "ns", knots = var_knots),
  arglag = list(fun = "ns", knots = lag_knots)
)

# Conditional count model with stratum fixed effects (weekly aggregated analogue).
# `eliminate = stratum` is equivalent to conditioning on stratum and is faster.
model <- gnm(
  deaths ~ cb_temp,
  eliminate = stratum,
  family = quasipoisson(),
  data = df,
  na.action = na.exclude
)

# First pass: find minimum-risk temperature for centering.
pred0 <- crosspred(cb_temp, model, by = 0.1)
mmt <- pred0$predvar[which.min(pred0$allRRfit)]

# Final cumulative association centered at MMT.
pred <- crosspred(cb_temp, model, cen = mmt, by = 0.1)

curve_df <- data.frame(
  temperature_C = pred$predvar,
  rr_cumulative = pred$allRRfit,
  rr_low_95 = pred$allRRlow,
  rr_high_95 = pred$allRRhigh,
  or_cumulative = pred$allRRfit,
  or_low_95 = pred$allRRlow,
  or_high_95 = pred$allRRhigh
)

pct_vals <- quantile(df$tmean_C_weekly, probs = c(0.01, 0.10, 0.50, 0.90, 0.99), na.rm = TRUE)
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
  rr_high_95 = curve_df$rr_high_95[nearest_idx],
  or_cumulative = curve_df$or_cumulative[nearest_idx],
  or_low_95 = curve_df$or_low_95[nearest_idx],
  or_high_95 = curve_df$or_high_95[nearest_idx]
)

write.csv(curve_df, out_curve, row.names = FALSE)
write.csv(pct_out, out_pct, row.names = FALSE)

png(out_plot, width = 1200, height = 760, res = 120)
plot(
  curve_df$temperature_C, curve_df$rr_cumulative,
  type = "l", lwd = 2, col = "#1b5e20",
  xlab = "Weekly Mean Temperature (C)",
  ylab = "Cumulative Relative Risk (RR)",
  main = "DLNM Cumulative Temperature-CKD Mortality Association"
)
lines(curve_df$temperature_C, curve_df$rr_low_95, col = "#43a047", lty = 2)
lines(curve_df$temperature_C, curve_df$rr_high_95, col = "#43a047", lty = 2)
abline(h = 1, lty = 3, col = "gray40")
abline(v = mmt, lty = 3, col = "#37474f")
grid(col = "gray85")
dev.off()

sink(out_summary)
cat("DLNM model summary\n")
cat("Input rows used:", nrow(df), "\n")
cat("States used:", length(unique(df$state)), "\n")
cat("Temperature restriction (1st-99th percentile): ",
    round(temp_bounds[[1]], 3), " to ", round(temp_bounds[[2]], 3), " C\n", sep = "")
cat("Max lag (weeks):", max_lag_weeks, "\n")
cat("Centering temperature (MMT):", round(mmt, 3), "C\n\n")
print(summary(model))
cat("\nCumulative RR at selected percentiles\n")
print(pct_out %>% select(percentile, temperature_C, rr_cumulative, rr_low_95, rr_high_95))
cat("\nCumulative OR at selected percentiles\n")
print(pct_out %>% select(percentile, temperature_C, or_cumulative, or_low_95, or_high_95))
sink()

cat("Modeling complete.\n")
cat("Wrote:", out_summary, "\n")
cat("Wrote:", out_curve, "\n")
cat("Wrote:", out_pct, "\n")
cat("Wrote:", out_plot, "\n")
