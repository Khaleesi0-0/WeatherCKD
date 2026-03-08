library(dplyr)
library(lubridate)
library(httr)
library(jsonlite)

ckd_files <- c("data/ckd2018.csv", "data/ckd2022.csv")
weather_file <- "data/us_state_daily_tmean_2004_2024.csv"

clean_ckd_out <- "data/ckd_2018_2022_cleaned.csv"
state_suppression_out <- "data/ckd_state_suppression_summary.csv"
merged_out <- "data/ckd_2018_2022_with_weather_weekly.csv"

# States above this share of suppressed/non-numeric Deaths are removed.
suppression_threshold <- 0.20

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

ckd_raw <- bind_rows(lapply(ckd_files, function(f) {
  read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}))

ckd <- ckd_raw %>%
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
    icd_subchapter = `MCD - ICD Sub-Chapter`,
    icd_subchapter_code = `MCD - ICD Sub-Chapter Code`,
    deaths_raw = trimws(as.character(Deaths)),
    is_suppressed = tolower(deaths_raw) == "suppressed",
    deaths = suppressWarnings(as.numeric(deaths_raw)),
    week_end_date = parse_week_end(mmwr_week),
    year = suppressWarnings(as.integer(sub("/.*", "", mmwr_week_code))),
    week = suppressWarnings(as.integer(sub(".*/", "", mmwr_week_code))),
    suppression_or_missing = is_suppressed | is.na(deaths)
  ) %>%
  filter(!is.na(state), state != "")

weather_start_date <- as.character(min(ckd$week_end_date, na.rm = TRUE) - 6)
weather_end_date <- as.character(max(ckd$week_end_date, na.rm = TRUE))

if (file.exists(weather_file)) {
  weather_daily <- read.csv(weather_file, stringsAsFactors = FALSE) %>%
    mutate(
      state = trimws(state),
      date = as.Date(date)
    )
} else {
  weather_daily <- data.frame(
    state = character(),
    date = as.Date(character()),
    tmean_C = numeric(),
    stringsAsFactors = FALSE
  )
}

required_states <- sort(unique(ckd$state))

have_coverage <- weather_daily %>%
  filter(state %in% required_states) %>%
  group_by(state) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    .groups = "drop"
  )

coverage_needed <- data.frame(state = required_states, stringsAsFactors = FALSE) %>%
  left_join(have_coverage, by = "state") %>%
  mutate(
    needs_fetch = is.na(min_date) | is.na(max_date) |
      min_date > as.Date(weather_start_date) |
      max_date < as.Date(weather_end_date)
  )

states_to_fetch <- coverage_needed %>%
  filter(needs_fetch) %>%
  pull(state)

if (length(states_to_fetch) > 0) {
  message("Fetching weather for ", length(states_to_fetch), " states to complete coverage...")
  coords_to_fetch <- state_coords %>% filter(state %in% states_to_fetch)

  fetched <- lapply(seq_len(nrow(coords_to_fetch)), function(i) {
    s <- coords_to_fetch$state[i]
    message("Weather: ", s)
    out <- fetch_state_daily(
      s,
      coords_to_fetch$lat[i],
      coords_to_fetch$lon[i],
      weather_start_date,
      weather_end_date
    )
    Sys.sleep(0.2)
    out
  })

  fetched_df <- bind_rows(fetched)
  if (nrow(fetched_df) > 0) {
    weather_daily <- bind_rows(weather_daily, fetched_df) %>%
      arrange(state, date) %>%
      distinct(state, date, .keep_all = TRUE)
    write.csv(weather_daily, weather_file, row.names = FALSE, na = "")
  }
}

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
  filter(!state %in% states_to_drop) %>%
  filter(!suppression_or_missing) %>%
  select(
    state,
    residence_state_code,
    mmwr_week,
    mmwr_week_code,
    year,
    week,
    week_end_date,
    icd_subchapter,
    icd_subchapter_code,
    deaths
  )

weather_daily <- weather_daily %>%
  mutate(
    state = trimws(state),
    date = as.Date(date),
    week_end_date = date + ((7 - wday(date, week_start = 7)) %% 7)
  )

weather_weekly <- weather_daily %>%
  group_by(state, week_end_date) %>%
  summarise(
    tmean_C_weekly = mean(tmean_C, na.rm = TRUE),
    n_days_weather = sum(!is.na(tmean_C)),
    .groups = "drop"
  )

ckd_weekly <- ckd_clean %>%
  group_by(
    state,
    residence_state_code,
    mmwr_week,
    mmwr_week_code,
    year,
    week,
    week_end_date,
    icd_subchapter,
    icd_subchapter_code
  ) %>%
  summarise(
    deaths = sum(deaths, na.rm = TRUE),
    .groups = "drop"
  )

merged <- ckd_weekly %>%
  left_join(weather_weekly, by = c("state", "week_end_date"))

write.csv(suppression_summary, state_suppression_out, row.names = FALSE, na = "")
write.csv(ckd_clean, clean_ckd_out, row.names = FALSE, na = "")
write.csv(merged, merged_out, row.names = FALSE, na = "")

cat("Suppression threshold:", suppression_threshold, "\n")
cat("States removed:", if (length(states_to_drop)) paste(states_to_drop, collapse = ", ") else "None", "\n")
cat("Clean CKD rows:", nrow(ckd_clean), "\n")
cat("Merged rows:", nrow(merged), "\n")
cat("Rows missing weekly weather after merge:", sum(is.na(merged$tmean_C_weekly)), "\n")
cat("Wrote:", clean_ckd_out, "\n")
cat("Wrote:", state_suppression_out, "\n")
cat("Wrote:", merged_out, "\n")
