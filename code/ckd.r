file_1999 <- "data/N181999.csv"
file_2018 <- "data/N182018.csv"
output_file <- "data/ckd_joined_1999_2018.csv"
plot_file <- "data/n18_aamr_trend_by_sex.png"

df_1999 <- read.csv(file_1999, check.names = FALSE, stringsAsFactors = FALSE)
df_2018 <- read.csv(file_2018, check.names = FALSE, stringsAsFactors = FALSE)

df_1999$source_file <- "N181999"
df_2018$source_file <- "N182018"

all_cols <- union(names(df_1999), names(df_2018))
missing_1999 <- setdiff(all_cols, names(df_1999))
missing_2018 <- setdiff(all_cols, names(df_2018))

for (col in missing_1999) df_1999[[col]] <- NA
for (col in missing_2018) df_2018[[col]] <- NA

df_1999 <- df_1999[, all_cols]
df_2018 <- df_2018[, all_cols]

joined_df <- rbind(df_1999, df_2018)

write.csv(joined_df, output_file, row.names = FALSE, na = "")

print(paste("Joined file written to:", output_file))

# Build AAMR trend data by sex (Female/Male only).
trend_df <- joined_df[
  joined_df$Sex %in% c("Female", "Male"),
  c("Year", "Sex", "Age Adjusted Rate")
]

trend_df$Year <- suppressWarnings(as.numeric(trend_df$Year))
trend_df$`Age Adjusted Rate` <- suppressWarnings(as.numeric(trend_df$`Age Adjusted Rate`))
trend_df <- trend_df[!is.na(trend_df$Year) & !is.na(trend_df$`Age Adjusted Rate`), ]

# If overlapping years exist across the two source files, average AAMR by Year + Sex.
trend_summary <- aggregate(
  `Age Adjusted Rate` ~ Year + Sex,
  data = trend_df,
  FUN = mean
)
trend_summary <- trend_summary[order(trend_summary$Year, trend_summary$Sex), ]

female <- trend_summary[trend_summary$Sex == "Female", ]
male <- trend_summary[trend_summary$Sex == "Male", ]

png(plot_file, width = 1200, height = 700, res = 120)
plot(
  female$Year, female$`Age Adjusted Rate`,
  type = "o", pch = 16, lwd = 2, col = "#D55E00",
  ylim = range(trend_summary$`Age Adjusted Rate`, na.rm = TRUE),
  xlab = "Year",
  ylab = "Age-Adjusted Mortality Rate (AAMR)",
  main = "N18 AAMR Trend by Sex"
)
lines(
  male$Year, male$`Age Adjusted Rate`,
  type = "o", pch = 17, lwd = 2, col = "#0072B2"
)
legend(
  "topleft",
  legend = c("Female", "Male"),
  col = c("#D55E00", "#0072B2"),
  pch = c(16, 17),
  lwd = 2,
  bty = "n"
)
grid(col = "gray85")
dev.off()

print(paste("AAMR trend plot written to:", plot_file))
