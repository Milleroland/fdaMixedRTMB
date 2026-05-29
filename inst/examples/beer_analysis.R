

# Real-data application for the beer tracker data.
#
# The script uses the same data in two curve formats:
#   1. one full-period cumulative curve per person, for a sex comparison;
#   2. one 28-day cumulative curve per person-month, for month comparisons
#      while keeping sex in the boundary-value model.
library(devtools)

devtools::load_all()


library(tidyverse)
library(lubridate)
library(googlesheets4)
library(patchwork)
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

month_cols <- c(
  January  = "#3366CC",
  February = "#FF9900",
  March    = "#0099C6",
  April    = "#DD4477",
  May      = "#B0B0B0"
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

last_beer_date <- max(beer_raw$date, na.rm = TRUE)
latest_month_day_28 <- as.Date(paste0(format(last_beer_date, "%Y-%m"), "-28"))
analysis_end_date <- max(last_beer_date, latest_month_day_28)
analysis_start_date <- as.Date(paste0(format(analysis_end_date, "%Y"), "-01-01"))
analysis_month <- as.integer(format(analysis_end_date, "%m"))
months_ordered <- month.name[seq_len(analysis_month)]
month_cols <- month_cols[months_ordered]
n_days_month <- min(28L, as.integer(format(analysis_end_date, "%d")))

date_grid <- tibble(
  date = seq(analysis_start_date, analysis_end_date, by = "day")
) %>%
  mutate(
    day = row_number(),
    week = (day - 1) / 7,
    month_num = month(date),
    month_name = factor(month.name[month_num], levels = months_ordered),
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

subject_info <- subject_lookup %>%
  select(subject, name, sex) %>%
  mutate(name = factor(name, levels = name))

subjects <- as.character(subject_info$name)

beer_df <- beer_daily %>%
  transmute(
    name_key = navn,
    date,
    day,
    week,
    month_num,
    month_name,
    day_within_month,
    beer_daily,
    beer_gain
  ) %>%
  left_join(subject_lookup, by = "name_key") %>%
  arrange(subject, day) %>%
  select(
    subject,
    name,
    sex,
    date,
    day,
    week,
    month_num,
    month_name,
    day_within_month,
    beer_daily,
    beer_gain
  ) %>%
  mutate(
    name = factor(name, levels = subjects),
    sex = factor(sex, levels = c("F", "M"))
  )

weeks <- beer_df %>%
  distinct(day, week) %>%
  arrange(day) %>%
  pull(week)

beer_gain_wide <- beer_df %>%
  select(week, name, beer_gain) %>%
  pivot_wider(names_from = name, values_from = beer_gain) %>%
  arrange(week)

beer_gain <- beer_gain_wide %>%
  select(all_of(subjects)) %>%
  as.matrix()
rownames(beer_gain) <- beer_gain_wide$week

beer_weekly <- beer_df %>%
  mutate(week_id = floor((day - min(day)) / 7) + 1L) %>%
  group_by(subject, name, sex, week_id) %>%
  summarise(
    beer_weekly = sum(beer_daily),
    n_days = n_distinct(day),
    week_midpoint = mean(week),
    .groups = "drop"
  ) %>%
  filter(n_days == 7L) %>%
  mutate(
    name = factor(name, levels = subjects),
    sex = factor(sex, levels = c("F", "M"))
  )

beer_monthly_df <- beer_df %>%
  filter(
    month_name %in% months_ordered,
    day_within_month <= n_days_month
  ) %>%
  mutate(month_name = factor(as.character(month_name), levels = months_ordered)) %>%
  arrange(subject, month_num, day_within_month) %>%
  group_by(subject, name, sex, month_num, month_name) %>%
  mutate(monthly_cumsum = cumsum(beer_daily)) %>%
  ungroup() %>%
  mutate(curve = paste(name, month_name, sep = ": "))

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
    name = factor(name, levels = subjects),
    sex = factor(sex, levels = c("F", "M")),
    month_name = factor(month_name, levels = months_ordered)
  )

beer_gain_monthly <- matrix(
  beer_monthly_df %>%
    arrange(curve, day_within_month) %>%
    pull(monthly_cumsum),
  nrow = n_days_month,
  ncol = nrow(curve_meta),
  dimnames = list(seq_len(n_days_month), curve_order)
)

data_summary <- list(
  subjects = subject_info,
  date_range = range(beer_df$date),
  months = months_ordered,
  days_per_month = n_days_month,
  n_full_curves = n_distinct(beer_df$name),
  n_monthly_curves = n_distinct(beer_monthly_df$curve)
)
print(data_summary)

# -- 2. Exploratory plots ----------------------------------------------------

eda_gain <- beer_df %>%
  transmute(
    subject,
    sex,
    week,
    response = "Cumulative beer gain",
    value = beer_gain
  )

eda_weekly <- beer_weekly %>%
  transmute(
    subject,
    sex,
    week = week_midpoint,
    response = "Weekly beers",
    value = beer_weekly
  )

eda_df <- bind_rows(eda_weekly, eda_gain) %>%
  mutate(
    response = factor(response, levels = c("Weekly beers", "Cumulative beer gain"))
  )

p_eda <- ggplot(eda_df, aes(week, value, group = subject, colour = sex)) +
  geom_line(alpha = 0.4, linewidth = 0.35) +
  stat_summary(
    aes(group = sex),
    fun = mean,
    geom = "line",
    linewidth = 1.05,
    alpha = 0.95
  ) +
  geom_hline(
    data = tibble(
      response = factor("Cumulative beer gain", levels = levels(eda_df$response))
    ),
    aes(yintercept = 0),
    inherit.aes = FALSE,
    colour = "grey55",
    linewidth = 0.35
  ) +
  facet_wrap(~ response, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = sex_cols) +
  labs(
    title = "Beer tracker trajectories",
    subtitle = "Weekly totals and cumulative beer gain with sex-specific mean curves",
    x = "Week since Jan 1",
    y = NULL,
    colour = NULL
  )

p_monthly_gain <- ggplot(
  beer_monthly_df,
  aes(day_within_month, monthly_cumsum, group = curve, colour = sex)
) +
  geom_line(alpha = 0.24, linewidth = 0.35) +
  stat_summary(
    aes(group = sex),
    fun = mean,
    geom = "line",
    linewidth = 1.05,
    alpha = 0.95
  ) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  facet_wrap(~ month_name, ncol = 3) +
  scale_colour_manual(values = sex_cols) +
  labs(
    title = "Monthly beer gain",
    subtitle = "Person-month curves use the first 28 days of each month",
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

print(p_eda)
print(p_monthly_gain)

# -- 3. Shared helpers -------------------------------------------------------

boundary_mean_grid <- function(fit, grid, keep_cols) {
  boundary_mf <- model.frame(
    fit$boundary_terms,
    data = grid,
    na.action = na.pass,
    xlev = fit$boundary_xlev
  )
  boundary_X <- model.matrix(
    fit$boundary_terms,
    boundary_mf,
    contrasts.arg = fit$boundary_contrasts
  )
  bc <- coef(fit, component = "boundary")
  grid[, keep_cols, drop = FALSE] %>%
    mutate(rho_model = as.numeric(boundary_X %*% bc))
}

variance_table <- function(fit) {
  vc <- coef(fit, component = "variance")
  tibble(
    term = c("tau", "sigma2", "sigma"),
    estimate = c(vc[["tau"]], vc[["sigma2"]], sqrt(vc[["sigma2"]]))
  )
}

observed_right_operator <- function(response_matrix, time, alpha) {
  n <- nrow(response_matrix)
  y_end <- response_matrix[n, ]
  dy_end <- (response_matrix[n, ] - response_matrix[n - 1L, ]) /
    (time[n] - time[n - 1L])
  alpha[1] * y_end + alpha[2] * dy_end
}

# -- 4. Full-period sex endpoint trend --------------------------------------

fit_sex_end_rate <- fdaLm_rtmb(
  beer_gain | name ~ 0,
  data = beer_df,
  boundary_formula = ~ sex,
  boundary_value = TRUE,
  operator = operator_laplace(),
  right_boundary = right_boundary_operator_k1(alpha = c(0, 1)),
  time_variable = "week"
)

sex_boundary_table <- function(fit) {
  rho_model <- boundary_mean_grid(
    fit,
    tibble(sex = factor(c("F", "M"), levels = c("F", "M"))),
    "sex"
  )

  tibble(
    term = c("Females", "Males", "Males minus females"),
    estimate = c(
      rho_model$rho_model[rho_model$sex == "F"],
      rho_model$rho_model[rho_model$sex == "M"],
      rho_model$rho_model[rho_model$sex == "M"] -
        rho_model$rho_model[rho_model$sex == "F"]
    )
  )
}

make_sex_diagnostics <- function(fit) {
  alpha <- fit$right_boundary$alpha
  rho_hat <- coef(fit, component = "rho")
  observed_operator <- observed_right_operator(beer_gain, weeks, alpha)

  subject_fit_df <- subject_info %>%
    mutate(
      rho_hat = as.numeric(rho_hat[as.character(name)]),
      observed_operator = as.numeric(observed_operator[as.character(name)])
    )

  beer_fit_df <- beer_df %>%
    mutate(
      fitted = as.numeric(fit$fitted),
      residual = beer_gain - fitted,
      serial_effect = as.numeric(fit$serial_effect),
      boundary_effect = as.numeric(fit$boundary_effect)
    )

  serial_mean_df <- beer_fit_df %>%
    group_by(sex, week) %>%
    summarise(serial_effect = mean(serial_effect), .groups = "drop")

  boundary_mean_df <- beer_fit_df %>%
    group_by(sex, week) %>%
    summarise(boundary_effect = mean(boundary_effect), .groups = "drop")

  rho_model_df <- boundary_mean_grid(
    fit,
    tibble(sex = factor(c("F", "M"), levels = c("F", "M"))),
    "sex"
  )

  y_lab_rho <- sprintf(
    "Right boundary operator: %.0f * value + %.0f * weekly derivative",
    alpha[1],
    alpha[2]
  )

  p_rho <- ggplot(subject_fit_df, aes(sex, observed_operator, colour = sex)) +
    geom_jitter(width = 0.09, height = 0, alpha = 0.62, size = 2.2) +
    geom_point(
      data = rho_model_df,
      aes(sex, rho_model, colour = sex),
      inherit.aes = FALSE,
      shape = 95,
      size = 13,
      stroke = 1.5
    ) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Full-period endpoint trend by sex",
      x = NULL,
      y = y_lab_rho,
      colour = NULL
    ) +
    theme(legend.position = "none")

  p_fit <- ggplot(beer_fit_df, aes(week, fitted, group = name, colour = sex)) +
    geom_line(alpha = 0.25, linewidth = 0.35) +
    stat_summary(aes(group = sex), fun = mean, geom = "line", linewidth = 1.15) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Fitted cumulative beer gain",
      x = "Week since Jan 1",
      y = "Fitted cumulative beers",
      colour = NULL
    )

  p_serial_mean <- ggplot(serial_mean_df, aes(week, serial_effect, colour = sex)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Mean serial effect by sex",
      x = "Week since Jan 1",
      y = "Mean serial effect",
      colour = NULL
    )

  p_boundary <- ggplot(boundary_mean_df, aes(week, boundary_effect, colour = sex)) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Mean boundary effect by sex",
      x = "Week since Jan 1",
      y = "Boundary effect",
      colour = NULL
    )

  p_resid_fitted <- ggplot(beer_fit_df, aes(fitted, residual, colour = sex)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_point(alpha = 0.28, size = 0.9) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.75,
      colour = "grey20"
    ) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Residuals versus fitted",
      x = "Fitted cumulative beers",
      y = "Residual",
      colour = NULL
    )

  list(
    boundary_results = sex_boundary_table(fit),
    variance_results = variance_table(fit),
    plots = list(
      rho = p_rho,
      fitted = p_fit,
      serial_mean = p_serial_mean,
      boundary = p_boundary,
      residual_vs_fitted = p_resid_fitted
    )
  )
}

sex_diagnostics <- make_sex_diagnostics(fit_sex_end_rate)

print(fit_sex_end_rate)
print(summary(fit_sex_end_rate))
plot(fit_sex_end_rate, which = c("fit", "residual", "rho", "boundary"))
print(sex_diagnostics$boundary_results)
print(sex_diagnostics$variance_results)
print(sex_diagnostics$plots$rho)
print(sex_diagnostics$plots$fitted)
print(sex_diagnostics$plots$serial_mean)
print(sex_diagnostics$plots$boundary)
print(sex_diagnostics$plots$residual_vs_fitted)

emm_sex_end_rate <- emmeans(
  fit_sex_end_rate,
  ~ sex,
  mode = "boundary"
)

emm_sex_end_rate_diff <- contrast(
  emm_sex_end_rate,
  method = list("M - F" = c(-1, 1))
)

print(emm_sex_end_rate)
print(emm_sex_end_rate_diff)
print(confint(emm_sex_end_rate_diff))
print(summary(emm_sex_end_rate_diff, infer = TRUE))

# -- 5. Monthly endpoint trends adjusted for sex -----------------------------

fit_month_end_rate <- fdaLm_rtmb(
  monthly_cumsum | curve ~ 0,
  data = beer_monthly_df,
  design = curve_meta,
  boundary_formula = ~ month_name + sex,
  boundary_value = TRUE,
  operator = operator_laplace(),
  right_boundary = right_boundary_operator_k1(alpha = c(0, 1)),
  time_variable = "day_within_month"
)

month_grid <- expand_grid(
  month_name = factor(months_ordered, levels = months_ordered),
  sex = factor(c("F", "M"), levels = c("F", "M"))
)

month_boundary_table <- function(fit) {
  boundary_mean_grid(fit, month_grid, c("month_name", "sex")) %>%
    group_by(month_name) %>%
    summarise(estimate = mean(rho_model), .groups = "drop") %>%
    transmute(month = as.character(month_name), estimate)
}

make_month_diagnostics <- function(fit) {
  alpha <- fit$right_boundary$alpha
  rho_hat <- coef(fit, component = "rho")
  observed_operator <- observed_right_operator(
    beer_gain_monthly,
    seq_len(n_days_month),
    alpha
  )

  curve_fit_df <- curve_meta %>%
    mutate(
      rho_hat = as.numeric(rho_hat[as.character(curve)]),
      observed_operator = as.numeric(observed_operator[as.character(curve)])
    )

  monthly_fit_df <- beer_monthly_df %>%
    mutate(
      fitted = as.numeric(fit$fitted),
      residual = monthly_cumsum - fitted,
      serial_effect = as.numeric(fit$serial_effect),
      boundary_effect = as.numeric(fit$boundary_effect)
    )

  serial_mean_df <- monthly_fit_df %>%
    group_by(month_name, day_within_month) %>%
    summarise(serial_effect = mean(serial_effect), .groups = "drop")

  boundary_mean_df <- monthly_fit_df %>%
    group_by(month_name, day_within_month) %>%
    summarise(boundary_effect = mean(boundary_effect), .groups = "drop")

  rho_model_df <- boundary_mean_grid(
    fit,
    month_grid,
    c("month_name", "sex")
  )

  y_lab_rho <- sprintf(
    "Right boundary operator: %.0f * value + %.0f * daily derivative",
    alpha[1],
    alpha[2]
  )

  p_rho <- ggplot(curve_fit_df, aes(month_name, observed_operator, colour = sex)) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.45),
      alpha = 0.62,
      size = 2
    ) +
    geom_point(
      data = rho_model_df,
      aes(month_name, rho_model, colour = sex),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.45),
      shape = 95,
      size = 11,
      stroke = 1.4
    ) +
    scale_colour_manual(values = sex_cols) +
    labs(
      title = "Monthly endpoint trends",
      subtitle = "Points are person-month endpoint operators; bars are fitted boundary means",
      x = NULL,
      y = y_lab_rho,
      colour = NULL
    )

  p_fit <- ggplot(
    monthly_fit_df,
    aes(day_within_month, fitted, group = curve, colour = month_name)
  ) +
    geom_line(alpha = 0.28, linewidth = 0.35) +
    stat_summary(aes(group = month_name), fun = mean, geom = "line", linewidth = 1.55) +
    scale_colour_manual(values = month_cols) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.8))) +
    labs(
      title = "Fitted monthly beer-gain curves",
      x = "Day within month",
      y = "Fitted cumulative beers",
      colour = NULL
    ) +
    theme(
      plot.title = element_text(size = 28),
      axis.title = element_text(size = 22),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 20)
    )

  p_serial <- ggplot(
    serial_mean_df,
    aes(day_within_month, serial_effect, colour = month_name)
  ) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = month_cols) +
    labs(
      title = "Mean serial effect by month",
      x = "Day within month",
      y = "Mean serial effect",
      colour = NULL
    )

  p_boundary <- ggplot(
    boundary_mean_df,
    aes(day_within_month, boundary_effect, colour = month_name)
  ) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = month_cols) +
    labs(
      title = "Mean boundary effect by month",
      x = "Day within month",
      y = "Boundary effect",
      colour = NULL
    )

  p_resid_fitted <- ggplot(monthly_fit_df, aes(fitted, residual, colour = month_name)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
    geom_point(alpha = 0.28, size = 0.9) +
    geom_smooth(
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.75,
      colour = "grey20"
    ) +
    scale_colour_manual(values = month_cols) +
    labs(
      title = "Monthly residuals versus fitted",
      x = "Fitted cumulative beers",
      y = "Residual",
      colour = NULL
    )

  list(
    boundary_results = month_boundary_table(fit),
    variance_results = variance_table(fit),
    plots = list(
      rho = p_rho,
      fitted = p_fit,
      serial_mean = p_serial,
      boundary = p_boundary,
      residual_vs_fitted = p_resid_fitted
    )
  )
}

month_diagnostics <- make_month_diagnostics(fit_month_end_rate)

print(fit_month_end_rate)
print(summary(fit_month_end_rate))
plot(fit_month_end_rate, which = c("fit", "residual", "rho", "boundary"))
print(month_diagnostics$boundary_results)
print(month_diagnostics$variance_results)
print(month_diagnostics$plots$rho)
print(month_diagnostics$plots$fitted)
print(month_diagnostics$plots$serial_mean)
print(month_diagnostics$plots$boundary)
print(month_diagnostics$plots$residual_vs_fitted)

emm_month_end_rate <- emmeans(
  fit_month_end_rate,
  ~ month_name,
  mode = "boundary"
)

emm_month_end_rate_pairs <- contrast(
  emm_month_end_rate,
  method = "pairwise"
)

emm_month_sex_end_rate <- emmeans(
  fit_month_end_rate,
  ~ sex,
  mode = "boundary"
)

print(emm_month_end_rate)
print(emm_month_end_rate_pairs)
print(confint(emm_month_end_rate))
print(summary(emm_month_end_rate_pairs, infer = TRUE))
print(emm_month_sex_end_rate)

predict(fit_month_end_rate)



# -- 6. Methods for the monthly endpoint-rate fit ---------------------------

# Basic print and summary methods
print(fit_month_end_rate)
print(summary(fit_month_end_rate))

# Coefficients
print(coef(fit_month_end_rate, component = "boundary"))
print(coef(fit_month_end_rate, component = "rho"))
print(coef(fit_month_end_rate, component = "variance"))

# Ordinary fixed effects
# This is empty here because the model formula is monthly_cumsum | curve ~ 0.
print(fixef(fit_month_end_rate))

# Variance-covariance matrices
print(vcov(fit_month_end_rate, component = "boundary"))
print(vcov(fit_month_end_rate, component = "all"))

# Fitted values and model components
print(head(predict(fit_month_end_rate), 10))
print(head(predict(fit_month_end_rate, type = "boundary"), 10))
print(head(predict(fit_month_end_rate, type = "serial"), 10))
print(head(predict(fit_month_end_rate, type = "rho"), 10))
print(head(predict(fit_month_end_rate, type = "right_derivative"), 10))

# Residuals
print(head(residuals(fit_month_end_rate), 10))
print(head(residuals(fit_month_end_rate, type = "pearson"), 10))

# Model dimensions and likelihood
print(nobs(fit_month_end_rate))
print(df.residual(fit_month_end_rate))
print(logLik(fit_month_end_rate))

# Model matrices
print(dim(model.matrix(fit_month_end_rate)))
print(dim(model.matrix(fit_month_end_rate, component = "all")))

# Terms object for ordinary fixed effects
# This is NULL here because the ordinary fixed-effect formula is ~ 0.
print(terms(fit_month_end_rate))

# Diagnostic plot method
plot(fit_month_end_rate, which = c("fit", "residual", "qq", "serial", "rho", "boundary"))

# emmeans compatibility for the boundary-value regression
emm_month_end_rate <- emmeans(
  fit_month_end_rate,
  ~ month_name,
  mode = "boundary"
)

emm_month_end_rate_pairs <- contrast(
  emm_month_end_rate,
  method = "pairwise"
)

emm_month_sex_end_rate <- emmeans(
  fit_month_end_rate,
  ~ sex,
  mode = "boundary"
)

emm_month_sex_end_rate_diff <- contrast(
  emm_month_sex_end_rate,
  method = "pairwise"
)

print(emm_month_end_rate)
print(emm_month_end_rate_pairs)
print(confint(emm_month_end_rate))
print(summary(emm_month_end_rate_pairs, infer = TRUE))
print(emm_month_sex_end_rate)
print(emm_month_sex_end_rate_diff)
