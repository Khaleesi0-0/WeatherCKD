suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(httr)
  library(jsonlite)
  library(dlnm)
  library(gnm)
})

# -----------------------------
# Inputs / outputs
# -----------------------------
input_files <- c(
  "data/raw/I99ckd2018.csv",
  "data/raw/I99ckd2020.csv",
  "data/raw/I99ckd2022.csv",
  "data/raw/I99ckd2024.csv"
)
weather_file <- "data/raw/us_state_daily_tmean_2004_2024.csv"

clean_out <- "data/processed/i99_n00n99_cleaned.csv"
suppression_out <- "results/tables/i99_n00n99_state_suppression_summary.csv"
merged_out <- "data/processed/i99_n00n99_with_weather_weekly.csv"

model_summary_out <- "results/tables/i99_n00n99_dlnm_model_summary.txt"
model_curve_out <- "results/tables/i99_n00n99_dlnm_cumulative_rr_curve.csv"
model_pct_out <- "results/tables/i99_n00n99_dlnm_cumulative_rr_percentiles.csv"
model_plot_out <- "results/figures/i99_n00n99_dlnm_cumulative_rr_plot.png"

suppression_threshold <- 0.20
max_lag_weeks <- 3

state_coords <- data.frame(
  state = c(
    "Alabama", "Alaska", "Arizona", "Arkansas", "California",
    "Colorado", "Connecticut", "Delaware", "Florida", "Georgia",
    "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa",
    "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland",
    "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri",
    "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey",
    "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio",
    "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina",
    "South Dakota", "Tennessee", "Texas", "Utah", "Vermont",
    "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming"
  ),
  lat = c(
    32.37, 58.30, 33.45, 34.75, 38.58,
    39.74, 41.76, 39.74, 30.44, 33.75,
    21.31, 43.62, 39.80, 39.77, 41.59,
    39.04, 38.20, 30.46, 44.31, 39.29,
    42.36, 42.73, 44.95, 32.30, 38.57,
    46.59, 40.81, 39.16, 43.21, 40.22,
    35.69, 42.65, 35.78, 46.81, 39.96,
    35.47, 44.94, 40.27, 41.82, 34.00,
    44.37, 36.17, 30.27, 40.76, 44.26,
    37.54, 47.04, 38.35, 43.07, 41.14
  ),
  lon = c(
    -86.30, -134.42, -112.07, -92.29, -121.49,
    -104.99, -72.68, -75.55, -84.28, -84.39,
    -157.86, -116.21, -89.65, -86.16, -93.62,
    -95.69, -84.87, -91.14, -69.78, -76.61,
    -71.06, -84.56, -93.09, -90.18, -92.17,
    -112.02, -96.70, -119.77, -71.54, -74.76,
    -105.94, -73.76, -78.64, -100.78, -82.99,
    -97.52, -123.03, -76.88, -71.41, -81.03,
    -100.35, -86.78, -97.74, -111.89, -72.58,
    -77.43, -122.90, -81.63, -89.40, -104.82
  ),
  stringsAsFactors = FALSE
)

parse_week_end <- function(x) {
  as.Date(sub(".* ending ", "", x), format = "%B %d, %Y")
}

fetch_state_daily <- function(state_name, lat, lon, start_date, end_date, retries = 3) {
  for (attempt in seq_len(retries)) {
    resp <- GET(
      url = "https://archive-api.open-meteo.com/v1/archive",
      query = list(
        latitude = lat,
        longitude = lon,
        start_date = start_date,
        end_date = end_date,
        daily = "temperature_2m_mean",
        timezone = "auto"
      ),
      timeout(120)
    )

    if (!http_error(resp)) {
      dat <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      if (!is.null(dat$daily)) {
        return(data.frame(
          state = state_name,
          date = as.Date(dat$daily$time),
          tmean_C = dat$daily$temperature_2m_mean,
          stringsAsFactors = FALSE
        ))
      }
    }

    Sys.sleep(0.5 * attempt)
  }

  warning("Failed weather fetch for ", state_name)
  NULL
}

# -----------------------------
# 1) Clean I99 files, keep N00-N99
# -----------------------------
raw <- bind_rows(lapply(input_files, function(f) {
  if (!file.exists(f)) stop("Missing input file: ", f)
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}))

ckd <- raw %>%
  filter(`MCD - ICD Chapter Code` == "N00-N99") %>%
  mutate(
    state = trimws(`Residence State`),
    residence_state_code = trimws(as.character(`Residence State Code`)),
    residence_state_code = ifelse(
      grepl("^[0-9]+$", residence_state_code),
      sprintf("%02d", as.integer(residence_state_code)),
      residence_state_code
    ),
    mmwr_week = `MMWR Week`,
    mmwr_week_code = `MMWR Week Code`,
    icd_chapter = `MCD - ICD Chapter`,
    icd_chapter_code = `MCD - ICD Chapter Code`,
    deaths_raw = trimws(as.character(Deaths)),
    is_suppressed = tolower(deaths_raw) == "suppressed",
    deaths = suppressWarnings(as.numeric(deaths_raw)),
    week_end_date = parse_week_end(mmwr_week),
    year = suppressWarnings(as.integer(sub("/.*", "", mmwr_week_code))),
    week = suppressWarnings(as.integer(sub(".*/", "", mmwr_week_code))),
    suppression_or_missing = is_suppressed | is.na(deaths)
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date))

suppression_summary <- ckd %>%
  group_by(state) %>%
  summarise(
    rows = n(),
    suppressed_rows = sum(is_suppressed, na.rm = TRUE),
    missing_or_suppressed_rows = sum(suppression_or_missing, na.rm = TRUE),
    suppression_share = mean(suppression_or_missing, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(suppression_share), desc(missing_or_suppressed_rows), state)

states_to_drop <- suppression_summary %>%
  filter(suppression_share > suppression_threshold) %>%
  pull(state)

ckd_clean <- ckd %>%
  filter(!state %in% states_to_drop, !suppression_or_missing) %>%
  group_by(
    state, residence_state_code, mmwr_week, mmwr_week_code,
    year, week, week_end_date, icd_chapter, icd_chapter_code
  ) %>%
  summarise(deaths = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
  arrange(state, week_end_date)

# -----------------------------
# 2) Ensure weather coverage, then merge
# -----------------------------
weather_start_date <- as.character(min(ckd_clean$week_end_date, na.rm = TRUE) - 6)
weather_end_date <- as.character(max(ckd_clean$week_end_date, na.rm = TRUE))
required_states <- sort(unique(ckd_clean$state))

if (file.exists(weather_file)) {
  weather_daily <- read.csv(weather_file, stringsAsFactors = FALSE) %>%
    mutate(state = trimws(state), date = as.Date(date))
} else {
  weather_daily <- data.frame(state = character(), date = as.Date(character()), tmean_C = numeric())
}

have_coverage <- weather_daily %>%
  filter(state %in% required_states) %>%
  group_by(state) %>%
  summarise(min_date = min(date, na.rm = TRUE), max_date = max(date, na.rm = TRUE), .groups = "drop")

coverage_needed <- data.frame(state = required_states, stringsAsFactors = FALSE) %>%
  left_join(have_coverage, by = "state") %>%
  mutate(
    needs_fetch = is.na(min_date) | is.na(max_date) |
      min_date > as.Date(weather_start_date) | max_date < as.Date(weather_end_date)
  )

states_to_fetch <- coverage_needed %>% filter(needs_fetch) %>% pull(state)
if (length(states_to_fetch) > 0) {
  coords_to_fetch <- state_coords %>% filter(state %in% states_to_fetch)
  fetched <- lapply(seq_len(nrow(coords_to_fetch)), function(i) {
    fetch_state_daily(
      coords_to_fetch$state[i],
      coords_to_fetch$lat[i],
      coords_to_fetch$lon[i],
      weather_start_date,
      weather_end_date
    )
  })
  fetched_df <- bind_rows(fetched)
  if (nrow(fetched_df) > 0) {
    weather_daily <- bind_rows(weather_daily, fetched_df) %>%
      arrange(state, date) %>%
      distinct(state, date, .keep_all = TRUE)
    write.csv(weather_daily, weather_file, row.names = FALSE, na = "")
  }
}

weather_weekly <- weather_daily %>%
  mutate(week_end_date = date + ((7 - wday(date, week_start = 7)) %% 7)) %>%
  group_by(state, week_end_date) %>%
  summarise(
    tmean_C_weekly = mean(tmean_C, na.rm = TRUE),
    n_days_weather = sum(!is.na(tmean_C)),
    .groups = "drop"
  )

merged <- ckd_clean %>%
  left_join(weather_weekly, by = c("state", "week_end_date"))

if (any(is.na(merged$tmean_C_weekly))) {
  warning("Missing weather rows after merge: ", sum(is.na(merged$tmean_C_weekly)))
}

write.csv(suppression_summary, suppression_out, row.names = FALSE, na = "")
write.csv(ckd_clean, clean_out, row.names = FALSE, na = "")
write.csv(merged, merged_out, row.names = FALSE, na = "")

# -----------------------------
# 3) DLNM model (same approach as prior run)
# -----------------------------
model_df <- merged %>%
  mutate(
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths)),
    tmean_C_weekly = suppressWarnings(as.numeric(tmean_C_weekly))
  ) %>%
  filter(!is.na(week_end_date), !is.na(deaths), !is.na(tmean_C_weekly), deaths >= 0)

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

write.csv(curve_df, model_curve_out, row.names = FALSE)
write.csv(pct_out, model_pct_out, row.names = FALSE)

png(model_plot_out, width = 1200, height = 760, res = 120)
plot(
  curve_df$temperature_C, curve_df$rr_cumulative,
  type = "l", lwd = 2, col = "#0d47a1",
  xlab = "Weekly Mean Temperature (C)",
  ylab = "Cumulative Relative Risk (RR)",
  main = "I99 N00-N99: DLNM Cumulative Temperature-Mortality Association"
)
lines(curve_df$temperature_C, curve_df$rr_low_95, col = "#1976d2", lty = 2)
lines(curve_df$temperature_C, curve_df$rr_high_95, col = "#1976d2", lty = 2)
abline(h = 1, lty = 3, col = "gray40")
abline(v = mmt, lty = 3, col = "#37474f")
grid(col = "gray85")
dev.off()

sink(model_summary_out)
cat("I99 N00-N99 DLNM model summary\n")
cat("Rows used in model:", nrow(model_df), "\n")
cat("States used:", length(unique(model_df$state)), "\n")
cat("Suppression threshold:", suppression_threshold, "\n")
cat("States removed:", if (length(states_to_drop)) paste(states_to_drop, collapse = ", ") else "None", "\n")
cat("Temperature restriction (1st-99th percentile):", round(temp_bounds[[1]], 3), "to", round(temp_bounds[[2]], 3), "C\n")
cat("Max lag (weeks):", max_lag_weeks, "\n")
cat("Centering temperature (MMT):", round(mmt, 3), "C\n\n")
print(summary(fit))
cat("\nCumulative RR at selected percentiles\n")
print(pct_out)
sink()

cat("Done. Outputs written:\n")
cat("-", clean_out, "\n")
cat("-", suppression_out, "\n")
cat("-", merged_out, "\n")
cat("-", model_summary_out, "\n")
cat("-", model_curve_out, "\n")
cat("-", model_pct_out, "\n")
cat("-", model_plot_out, "\n")
