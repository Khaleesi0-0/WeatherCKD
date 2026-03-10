suppressPackageStartupMessages({
  library(dplyr)
})

curves_file <- "results/tables/n18_state_cluster_period_curves.csv"
assignments_file <- "results/tables/n18_state_cluster_assignments.csv"
ckd_file <- "data/processed/ckd_2018_2022_cleaned.csv"
weather_file <- "data/processed/ckd_2018_2022_with_weather_weekly.csv"
png_out <- "results/figures/n18_pooled_temperature_rr_by_cluster.png"
pdf_out <- "results/figures/n18_pooled_temperature_rr_by_cluster.pdf"

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

curves <- read.csv(curves_file, stringsAsFactors = FALSE)
assignments <- read.csv(assignments_file, stringsAsFactors = FALSE)
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
reference_temp_c <- ckd %>%
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
  summarise(ref = as.numeric(stats::quantile(tmean_C_weekly, probs = 0.50, na.rm = TRUE, names = FALSE))) %>%
  pull(ref)

cluster_labels <- assignments %>%
  group_by(cluster) %>%
  summarise(states = paste(sort(state), collapse = ", "), .groups = "drop") %>%
  mutate(
    label = dplyr::case_when(
      cluster == "cluster_1" ~ "Cluster 1: hot southern/subtropical",
      cluster == "cluster_2" ~ "Cluster 2: warm transitional/mixed",
      cluster == "cluster_3" ~ "Cluster 3: cool western/coastal-mountain",
      cluster == "cluster_4" ~ "Cluster 4: cold continental/Midwest-Appalachian",
      TRUE ~ cluster
    )
  )

curves <- curves %>%
  left_join(cluster_labels %>% select(cluster, label), by = "cluster") %>%
  mutate(
    component = factor(component, levels = c("lag0", "lag13")),
    period = factor(period, levels = c("2018-2021", "2022-2025")),
    cluster = factor(cluster, levels = paste0("cluster_", 1:4))
  )

cluster_colors <- c(
  "cluster_1" = "#C44E52",
  "cluster_2" = "#DD8452",
  "cluster_3" = "#4C72B0",
  "cluster_4" = "#55A868"
)

period_lty <- c(
  "2018-2021" = 1,
  "2022-2025" = 2
)

plot_component <- function(component_name, panel_title) {
  df <- curves %>% filter(component == component_name)

  if (nrow(df) == 0) {
    plot.new()
    title(panel_title)
    text(0.5, 0.5, paste("No data for", component_name))
    return(invisible(NULL))
  }

  y_limits <- range(c(df$rr_low_95, df$rr_high_95), finite = TRUE)
  x_limits <- range(df$temperature_c, finite = TRUE)

  plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Temperature (C)",
    ylab = "Relative Risk (RR)",
    main = panel_title,
    las = 1
  )
  grid(col = "gray88")
  abline(h = 1, lty = 3, col = "gray45")
  abline(v = reference_temp_c, lty = 3, col = "gray45")

  for (cl in levels(df$cluster)) {
    for (pd in levels(df$period)) {
      sub_df <- df %>%
        filter(cluster == cl, period == pd) %>%
        arrange(temperature_c)

      if (nrow(sub_df) == 0) {
        next
      }

      lines(
        sub_df$temperature_c,
        sub_df$rr,
        col = cluster_colors[[cl]],
        lty = period_lty[[pd]],
        lwd = 2
      )
    }
  }
}

draw_legends <- function() {
  legend(
    "topleft",
    legend = cluster_labels$label[match(names(cluster_colors), cluster_labels$cluster)],
    col = unname(cluster_colors),
    lty = 1,
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
  legend(
    "topright",
    legend = names(period_lty),
    col = "black",
    lty = unname(period_lty),
    lwd = 2,
    bty = "n",
    cex = 0.8,
    title = "Period"
  )
}

make_figure <- function(device_fun, file, width, height) {
  device_fun(file, width = width, height = height)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 4, 1) + 0.1, oma = c(0, 0, 2, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  plot_component("lag0", "A. Pooled lag 0 relationship")
  draw_legends()
  plot_component("lag13", "B. Pooled lag 1-3 week average relationship")
  mtext(paste0("Temperature and N18 mortality pooled curves by cluster (reference = ", round(reference_temp_c, 1), " C)"), outer = TRUE, cex = 1.1)
}

make_figure(
  function(file, width, height) png(file, width = width, height = height, res = 140),
  png_out,
  width = 2200,
  height = 900
)

make_figure(
  pdf,
  pdf_out,
  width = 14,
  height = 6
)

cat("Wrote:", png_out, "\n")
cat("Wrote:", pdf_out, "\n")
