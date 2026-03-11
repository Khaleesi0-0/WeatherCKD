suppressPackageStartupMessages({
  library(dplyr)
})

curves_file <- "results/tables/n18_state_cluster_period_curves.csv"
effects_file <- "results/tables/n18_state_cluster_period_effects.csv"
assignments_file <- "results/tables/n18_state_cluster_assignments.csv"
png_out <- "results/figures/n18_pooled_temperature_or_by_cluster.png"
pdf_out <- "results/figures/n18_pooled_temperature_or_by_cluster.pdf"

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

curves <- read.csv(curves_file, stringsAsFactors = FALSE) %>%
  mutate(
    cluster = trimws(cluster),
    period = trimws(period)
  )

effects <- read.csv(effects_file, stringsAsFactors = FALSE) %>%
  mutate(
    cluster = trimws(cluster),
    period = trimws(period)
  )

assignments <- read.csv(assignments_file, stringsAsFactors = FALSE) %>%
  mutate(state = trimws(state), cluster = trimws(cluster))

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
    period = factor(period, levels = c("2018-2021", "2022-2025")),
    cluster = factor(cluster, levels = sort(unique(cluster))),
    label = factor(label, levels = cluster_labels$label[match(levels(cluster), cluster_labels$cluster)])
  )

effects <- effects %>%
  mutate(
    period = factor(period, levels = c("2018-2021", "2022-2025")),
    cluster = factor(cluster, levels = levels(curves$cluster))
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

period_fill <- c(
  "2018-2021" = "#D9E7F5",
  "2022-2025" = "#FBE5D6"
)

plot_cluster_panel <- function(cluster_name) {
  df <- curves %>%
    filter(cluster == cluster_name) %>%
    filter(is.finite(or_cumulative), is.finite(or_low_95), is.finite(or_high_95)) %>%
    arrange(period, temperature_c)

  meta_df <- effects %>% filter(cluster == cluster_name)
  cluster_label <- cluster_labels$label[match(as.character(cluster_name), cluster_labels$cluster)]

  if (nrow(df) == 0) {
    plot.new()
    title(main = cluster_label)
    text(0.5, 0.5, "No finite pooled curve data")
    return(invisible(NULL))
  }

  x_limits <- range(df$temperature_c, finite = TRUE)
  y_limits <- range(c(df$or_low_95, df$or_high_95), finite = TRUE)
  y_limits[1] <- min(y_limits[1], 1)
  y_limits[2] <- max(y_limits[2], 1)

  plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Temperature (C)",
    ylab = "Cumulative OR",
    main = cluster_label,
    las = 1
  )
  grid(col = "gray88")
  abline(h = 1, lty = 3, col = "gray45")

  for (pd in levels(df$period)) {
    band_df <- df %>%
      filter(period == pd) %>%
      arrange(temperature_c)

    if (nrow(band_df) == 0) {
      next
    }

    fill_col <- grDevices::adjustcolor(period_fill[[pd]], alpha.f = 0.45)
    polygon(
      x = c(band_df$temperature_c, rev(band_df$temperature_c)),
      y = c(band_df$or_low_95, rev(band_df$or_high_95)),
      col = fill_col,
      border = NA
    )
  }

  for (pd in levels(df$period)) {
    line_df <- df %>%
      filter(period == pd) %>%
      arrange(temperature_c)

    if (nrow(line_df) == 0) {
      next
    }

    lines(
      line_df$temperature_c,
      line_df$or_cumulative,
      col = cluster_colors[[as.character(cluster_name)]],
      lty = period_lty[[pd]],
      lwd = 2.5
    )

    center_x <- meta_df$centering_temp_c[meta_df$period == pd][1]
    if (is.finite(center_x)) {
      abline(v = center_x, lty = period_lty[[pd]], col = "gray55")
    }
  }

  legend(
    "topright",
    legend = c(
      paste0("States: ", nrow(assignments %>% filter(cluster == cluster_name))),
      paste0("2018-2021 center: ", round(meta_df$centering_temp_c[meta_df$period == "2018-2021"][1], 1), " C"),
      paste0("2022-2025 center: ", round(meta_df$centering_temp_c[meta_df$period == "2022-2025"][1], 1), " C")
    ),
    bty = "n",
    cex = 0.75
  )
}

draw_legends <- function() {
  legend(
    "bottomleft",
    legend = names(period_lty),
    col = "black",
    lty = unname(period_lty),
    lwd = 2.5,
    bty = "n",
    cex = 0.85,
    title = "Period"
  )
  legend(
    "bottomright",
    legend = c("95% CI: 2018-2021", "95% CI: 2022-2025"),
    fill = grDevices::adjustcolor(unname(period_fill), alpha.f = 0.45),
    border = NA,
    bty = "n",
    cex = 0.85
  )
}

make_figure <- function(device_fun, file, width, height) {
  device_fun(file, width = width, height = height)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 3, 0))
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  for (cl in levels(curves$cluster)) {
    plot_cluster_panel(cl)
  }

  draw_legends()
  mtext(
    "N18 pooled cluster-period temperature curves\nDLNM cumulative odds-ratio-scale associations within cluster-specific 1st-99th percentile ranges",
    outer = TRUE,
    cex = 1.05
  )
}

make_figure(
  function(file, width, height) png(file, width = width, height = height, res = 140),
  png_out,
  width = 2200,
  height = 1600
)

make_figure(
  pdf,
  pdf_out,
  width = 14,
  height = 10
)

cat("Wrote:", png_out, "\n")
cat("Wrote:", pdf_out, "\n")
