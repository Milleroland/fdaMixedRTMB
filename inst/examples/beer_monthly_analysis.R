# ============================================================
# beer_monthly_analysis.R
#
# Monthly beer gain analysis.
#
# The curves are cumulative beer gain within each month.  Since the current
# data end on the 28th, every month is cut at day 28 so the months are compared
# on the same time grid.
# ============================================================

library(tidyverse)
library(lubridate)
library(googlesheets4)
library(emmeans)
library(fdaMixedRTMB)

theme_set(
  theme_minimal(base_size = 20) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 28),
      plot.subtitle = element_text(color = "grey35", size = 20),
      strip.text = element_text(face = "bold", size = 20),
      axis.title = element_text(color = "grey20", size = 22),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)
    )
)

sex_cols <- c(
  F = "#D55E00",
  M = "#0072B2"
)

sheet_url <- "https://docs.google.com/spreadsheets/d/1Nz0i7Q9p61VT15Ia4E5muW4npkYCniaKVtAWBzTnzP8/edit?gid=0#gid=0"

sex_lookup <- tibble(
  navn = c(
    "Mille", "Milla", "Augusta", "Lise", "Frigg",
    "Baldur", "Julle", "Malthe", "Nikolai", "William"
  ),
  sex = c(rep("F", 5), rep("M", 5))
)

parse_tracker_date <- function(x) {
  as.Date(parse_date_time(
    as.character(x),
    orders = c("dmy HMS", "dmy HM", "dmy", "ymd HMS", "ymd HM", "ymd"),
    quiet = TRUE
  ))
}

parse_beer_count <- function(x) {
  as.numeric(gsub(",", ".", as.character(x), fixed = TRUE))
}

# -- 1. Data preparation -----------------------------------------------------

gs4_deauth()

beer_raw <- read_sheet(sheet_url) %>%
  rename_with(tolower) %>%
  mutate(
    date = parse_tracker_date(tid),
    beer = replace_na(parse_beer_count(antal), 1)
  ) %>%
  select(-any_of("sex")) %>%
  left_join(sex_lookup, by = "navn") %>%
  transmute(
    navn = as.character(navn),
    date,
    beer,
    sex = factor(sex, levels = c("F", "M"))
  )

analysis_year <- as.integer(format(max(beer_raw$date, na.rm = TRUE), "%Y"))
analysis_start_date <- as.Date(paste0(analysis_year, "-01-01"))
analysis_end_date <- as.Date(paste0(analysis_year, "-05-28"))
n_days_month <- 28L
months_ordered <- month.name[1:5]

date_grid <- tibble(
  date = seq(analysis_start_date, analysis_end_date, by = "day")
) %>%
  mutate(
    day = row_number(),
    month_num = month(date),
    month_name = factor(month.name[month_num], levels = month.name),
    day_within_month = mday(date)
  )

beer_daily <- beer_raw %>%
  filter(date >= analysis_start_date, date <= analysis_end_date) %>%
  group_by(navn, sex, date) %>%
  summarise(beer_daily = sum(beer), .groups = "drop") %>%
  complete(
    nesting(navn, sex),
    date = date_grid$date,
    fill = list(beer_daily = 0)
  ) %>%
  left_join(date_grid, by = "date") %>%
  arrange(navn, day) %>%
  group_by(navn) %>%
  mutate(beer_gain = cumsum(beer_daily)) %>%
  ungroup()

subject_lookup <- beer_daily %>%
  distinct(name_key = navn, sex) %>%
  arrange(name_key) %>%
  mutate(
    subject = row_number(),
    name = paste("Friend", subject),
    sex = factor(sex, levels = c("F", "M"))
  )

subjects <- subject_lookup$name

beer_monthly_df <- beer_daily %>%
  transmute(
    name_key = navn,
    date,
    month_num,
    month_name,
    day_within_month,
    beer_daily
  ) %>%
  left_join(subject_lookup, by = "name_key") %>%
  filter(
    as.character(month_name) %in% months_ordered,
    day_within_month <= n_days_month
  ) %>%
  arrange(subject, month_num, day_within_month) %>%
  group_by(subject, name, sex, month_num, month_name) %>%
  mutate(monthly_cumsum = cumsum(beer_daily)) %>%
  ungroup() %>%
  mutate(
    name = factor(name, levels = subjects),
    sex = factor(sex, levels = c("F", "M")),
    month_name = factor(as.character(month_name), levels = months_ordered),
    curve = paste(name, month_name, sep = ": ")
  )

curve_meta <- beer_monthly_df %>%
  distinct(curve, subject, name, sex, month_name, month_num) %>%
  arrange(subject, month_num) %>%
  mutate(
    curve = factor(curve, levels = curve),
    name = factor(name, levels = subjects),
    sex = factor(sex, levels = c("F", "M")),
    month_name = factor(month_name, levels = months_ordered)
  )

curve_order <- as.character(curve_meta$curve)

beer_monthly_df <- beer_monthly_df %>%
  mutate(
    curve = factor(curve, levels = curve_order),
    month_name = factor(month_name, levels = months_ordered)
  )

data_summary <- list(
  date_range = range(beer_daily$date),
  months = months_ordered,
  days_per_month = n_days_month,
  subjects = n_distinct(beer_monthly_df$name),
  curves = n_distinct(beer_monthly_df$curve)
)
print(data_summary)

# -- 2. Exploratory plot -----------------------------------------------------

p_monthly_gain <- ggplot(
  beer_monthly_df,
  aes(day_within_month, monthly_cumsum, group = curve, colour = sex)
) +
  geom_line(alpha = 0.25, linewidth = 0.35) +
  stat_summary(
    aes(group = sex),
    fun = mean,
    geom = "line",
    linewidth = 1,
    alpha = 0.95
  ) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  facet_wrap(~ month_name, ncol = 3) +
  scale_colour_manual(values = sex_cols) +
  labs(
    title = "Monthly beer gain",
    subtitle = "Cumulative beer gain within each month, using the first 28 days",
    x = "Day within month",
    y = "Cumulative beers",
    colour = NULL
  ) +
  theme(
    plot.title = element_text(size = 28),
    plot.subtitle = element_text(size = 20),
    strip.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 20)
  )

print(p_monthly_gain)

# -- 3. Monthly boundary-value models ---------------------------------------

fit_month_end_rate <- fdaLm_rtmb(
  monthly_cumsum | curve ~ 0,
  data = beer_monthly_df,
  design = curve_meta,
  boundary_formula = ~ month_name,
  boundary_value = TRUE,
  operator = operator_laplace(),
  right_boundary = right_boundary_operator_k1(alpha = c(0, 1)),
  time_variable = "day_within_month"
)

fit_month_total <- fdaLm_rtmb(
  monthly_cumsum | curve ~ 0,
  data = beer_monthly_df,
  design = curve_meta,
  boundary_formula = ~ month_name,
  boundary_value = TRUE,
  operator = operator_laplace(),
  right_boundary = right_boundary_operator_k1(alpha = c(1, 0)),
  time_variable = "day_within_month"
)

# -- 4. Package methods ------------------------------------------------------

# End-rate boundary model
print(fit_month_end_rate)
print(summary(fit_month_end_rate))
plot(fit_month_end_rate, which = c("fit", "residual", "rho", "boundary"))

# Total-gain boundary model
print(fit_month_total)
print(summary(fit_month_total))
plot(fit_month_total, which = c("fit", "residual", "rho", "boundary"))

# A few standard methods
print(coef(fit_month_end_rate, component = "boundary"))
print(coef(fit_month_end_rate, component = "variance"))
print(head(predict(fit_month_end_rate)))
print(head(residuals(fit_month_end_rate)))
print(logLik(fit_month_end_rate))
print(nobs(fit_month_end_rate))
print(df.residual(fit_month_end_rate))

# emmeans uses the boundary-value methods in fdaMixedRTMB.
emm_month_end_rate <- emmeans(
  fit_month_end_rate,
  ~ month_name,
  mode = "boundary"
)

emm_month_total <- emmeans(
  fit_month_total,
  ~ month_name,
  mode = "boundary"
)

print(emm_month_end_rate)
print(pairs(emm_month_end_rate))
print(emm_month_total)
