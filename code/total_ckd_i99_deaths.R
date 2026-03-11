suppressPackageStartupMessages({
  library(dplyr)
})

ckd_file <- "data/processed/ckd_2018_2022_cleaned.csv"
i99_file <- "data/processed/i99_n00n99_cleaned.csv"

summary_out <- "results/tables/ckd_i99_total_deaths_summary.csv"
yearly_out <- "results/tables/ckd_i99_total_deaths_by_year.csv"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

ckd <- read.csv(ckd_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths)),
    year = suppressWarnings(as.integer(year))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(deaths), deaths >= 0) %>%
  group_by(state, week_end_date, mmwr_week_code, year) %>%
  summarise(ckd_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

i99 <- read.csv(i99_file, stringsAsFactors = FALSE) %>%
  mutate(
    state = trimws(state),
    week_end_date = as.Date(week_end_date),
    deaths = suppressWarnings(as.numeric(deaths)),
    year = suppressWarnings(as.integer(year))
  ) %>%
  filter(!is.na(state), state != "", !is.na(week_end_date), !is.na(deaths), deaths >= 0) %>%
  group_by(state, week_end_date, mmwr_week_code, year) %>%
  summarise(i99_deaths = sum(deaths, na.rm = TRUE), .groups = "drop")

matched <- inner_join(
  ckd,
  i99,
  by = c("state", "week_end_date", "mmwr_week_code", "year")
) %>%
  mutate(
    ckd_minus_i99_deaths = ckd_deaths - i99_deaths
  )

summary_df <- data.frame(
  metric = c(
    "ckd_total_deaths_raw",
    "i99_total_deaths_raw",
    "ckd_total_deaths_matched_overlap",
    "i99_total_deaths_matched_overlap",
    "ckd_minus_i99_total_deaths",
    "matched_state_weeks"
  ),
  value = c(
    sum(ckd$ckd_deaths, na.rm = TRUE),
    sum(i99$i99_deaths, na.rm = TRUE),
    sum(matched$ckd_deaths, na.rm = TRUE),
    sum(matched$i99_deaths, na.rm = TRUE),
    sum(matched$ckd_minus_i99_deaths, na.rm = TRUE),
    nrow(matched)
  )
)

yearly_df <- matched %>%
  group_by(year) %>%
  summarise(
    ckd_deaths = sum(ckd_deaths, na.rm = TRUE),
    i99_deaths = sum(i99_deaths, na.rm = TRUE),
    ckd_minus_i99_deaths = sum(ckd_minus_i99_deaths, na.rm = TRUE),
    matched_state_weeks = n(),
    .groups = "drop"
  ) %>%
  arrange(year)

write.csv(summary_df, summary_out, row.names = FALSE)
write.csv(yearly_df, yearly_out, row.names = FALSE)

cat("Wrote:", summary_out, "\n")
cat("Wrote:", yearly_out, "\n")
