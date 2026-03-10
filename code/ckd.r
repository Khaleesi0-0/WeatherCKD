
# -------------------------
# Weather data collection
# -------------------------
library(noaa)
library(httr)
library(jsonlite)
library(dplyr)
library(lubridate)
library(tidyr)

weather_output_file <- "data/raw/us_state_daily_tmean_2004_2024.csv"
start_date <- "2018-01-01"
end_date <- "2025-12-31"

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

get_state_daily <- function(state_name, lat, lon) {
  message("Fetching daily data for ", state_name, " ...")

  resp <- GET(
    url = "https://archive-api.open-meteo.com/v1/archive",
    query = list(
      latitude = lat,
      longitude = lon,
      start_date = start_date,
      end_date = end_date,
      daily = "temperature_2m_mean",
      timezone = "auto"
    )
  )

  if (http_error(resp)) {
    stop(
      "Weather error for ", state_name, ": ",
      status_code(resp), " - ", content(resp, "text")
    )
  }

  dat <- fromJSON(content(resp, "text", encoding = "UTF-8"))

  if (is.null(dat$daily)) {
    warning("No daily data for ", state_name)
    return(NULL)
  }

  tibble(
    state = state_name,
    date = as.Date(dat$daily$time),
    tmean_C = dat$daily$temperature_2m_mean
  )
}

weather_list <- lapply(seq_len(nrow(state_coords)), function(i) {
  tryCatch(
    get_state_daily(state_coords$state[i], state_coords$lat[i], state_coords$lon[i]),
    error = function(e) {
      warning(conditionMessage(e))
      NULL
    }
  )
})

weather_df <- bind_rows(weather_list)
weather_df <- weather_df %>% arrange(state, date)

write.csv(weather_df, weather_output_file, row.names = FALSE, na = "")
print(paste("Weather file written to:", weather_output_file))
