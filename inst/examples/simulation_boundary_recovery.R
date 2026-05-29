library(fdaMixedRTMB)
library(ggplot2)
library(patchwork)

# True parameters
true_gamma  <- c(gamma0 = 0.5, gamma1 = 0.3)
true_tau    <- 1.0
true_sigma2 <- 0.10

left_bc  <- c(1, 0)
right_bc <- c(1, 0)

# Data-generating function
simulate_one <- function(n_curves, n_time, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  time_grid <- seq(0, 1, length.out = n_time)
  x_cov     <- seq(-1, 1, length.out = n_curves)
  true_rho  <- true_gamma["gamma0"] + true_gamma["gamma1"] * x_cov
  h_vals    <- time_grid

  y     <- numeric(n_curves * n_time)
  curve <- factor(rep(seq_len(n_curves), each = n_time))
  time  <- rep(time_grid, n_curves)

  for (i in seq_len(n_curves)) {
    idx <- (i - 1L) * n_time + seq_len(n_time)
    amp <- rnorm(1, sd = sqrt(true_tau) * 0.3)
    y[idx] <- h_vals * true_rho[i] +
              amp * sin(pi * time_grid) +
              rnorm(n_time, sd = sigma)
  }

  data.frame(
    y          = y,
    time       = time,
    curve      = curve,
    boundary_x = rep(x_cov, each = n_time)
  )
}

# Sanity check
example_data <- simulate_one(n_curves = 8L, n_time = 20L,
                              sigma = sqrt(true_sigma2), seed = 42L)

example_fit <- fdaLm_rtmb(
  y | curve ~ 0,
  data             = example_data,
  boundary_formula = ~ boundary_x,
  left_boundary    = left_bc,
  right_boundary   = right_bc,
  operator         = operator_k1(lambda = 1),
  time_variable    = "time"
)
summary(example_fit)

coef(example_fit, component = "boundary")
true_gamma

# Simulation engine
run_scenario <- function(n_curves, n_time, sigma, n_sim = 200L, base_seed = 1000L) {
  results <- vector("list", n_sim)

  for (i in seq_len(n_sim)) {
    dat <- simulate_one(n_curves, n_time, sigma, seed = base_seed + i)

    fit <- tryCatch(
      fdaLm_rtmb(
        y | curve ~ 0,
        data             = dat,
        boundary_formula = ~ boundary_x,
        left_boundary    = left_bc,
        right_boundary   = right_bc,
        operator         = operator_k1(lambda = 1),
        time_variable    = "time",
        control          = list(eval.max = 500L, iter.max = 300L)
      ),
      error = function(e) NULL
    )

    if (is.null(fit) || !isTRUE(fit$converged)) {
      results[[i]] <- NULL
      next
    }

    bc_hat  <- unname(coef(fit, component = "boundary"))
    bc_vcov <- vcov(fit, component = "boundary")
    bc_se   <- sqrt(pmax(diag(bc_vcov), 0))
    ci_lo   <- bc_hat - 1.96 * bc_se
    ci_hi   <- bc_hat + 1.96 * bc_se
    covered <- (ci_lo <= unname(true_gamma)) & (unname(true_gamma) <= ci_hi)

    results[[i]] <- data.frame(
      rep        = i,
      converged  = TRUE,
      gamma0_hat = bc_hat[1L],    gamma1_hat = bc_hat[2L],
      gamma0_se  = bc_se[1L],     gamma1_se  = bc_se[2L],
      gamma0_lo  = ci_lo[1L],     gamma0_hi  = ci_hi[1L],
      gamma1_lo  = ci_lo[2L],     gamma1_hi  = ci_hi[2L],
      gamma0_cov = covered[1L],   gamma1_cov = covered[2L],
      sigma2_hat = fit$sigma2
    )
  }

  do.call(rbind, Filter(Negate(is.null), results))
}

# Scenarios
baseline_sigma <- sqrt(true_sigma2)

scenarios <- data.frame(
  label    = c("nc=6", "nc=10", "nc=20",
               "nt=10", "nt=20", "nt=40",
               "sd=0.1", "sd=0.3", "sd=0.6"),
  n_curves = c(6L, 10L, 20L, 10L, 10L, 10L, 10L, 10L, 10L),
  n_time   = c(20L, 20L, 20L, 10L, 20L, 40L, 20L, 20L, 20L),
  sigma    = c(rep(baseline_sigma, 3L), rep(baseline_sigma, 3L), 0.1, 0.3, 0.6),
  stringsAsFactors = FALSE
)

n_sim_per_scenario <- 500L

sim_results <- vector("list", nrow(scenarios))
for (s in seq_len(nrow(scenarios))) {
  sim_results[[s]] <- run_scenario(
    n_curves  = scenarios$n_curves[s],
    n_time    = scenarios$n_time[s],
    sigma     = scenarios$sigma[s],
    n_sim     = n_sim_per_scenario,
    base_seed = s * 1000L
  )
}
names(sim_results) <- scenarios$label

# Summary statistics
compute_summary <- function(res, scenario_row) {
  data.frame(
    label     = scenario_row$label,
    n_curves  = scenario_row$n_curves,
    n_time    = scenario_row$n_time,
    sigma     = round(scenario_row$sigma, 3L),
    n_conv    = nrow(res),
    bias_g0   = mean(res$gamma0_hat - true_gamma["gamma0"]),
    rmse_g0   = sqrt(mean((res$gamma0_hat - true_gamma["gamma0"])^2)),
    cov_g0    = mean(res$gamma0_cov),
    bias_g1   = mean(res$gamma1_hat - true_gamma["gamma1"]),
    rmse_g1   = sqrt(mean((res$gamma1_hat - true_gamma["gamma1"])^2)),
    cov_g1    = mean(res$gamma1_cov),
    bias_sig2 = mean(res$sigma2_hat - true_sigma2),
    rmse_sig2 = sqrt(mean((res$sigma2_hat - true_sigma2)^2)),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(
  rbind,
  mapply(compute_summary, sim_results, split(scenarios, seq_len(nrow(scenarios))),
         SIMPLIFY = FALSE)
)

summary_table

# Plots
theme_sim <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    legend.position  = "bottom"
  )

col_g0   <- "#3B82F6"
col_g1   <- "#F97316"
col_true <- "#EF4444"
col_hit  <- "#22C55E"
col_miss <- "#94A3B8"

group_labels <- c(
  n_curves = "Number of curves",
  n_time   = "Time points per curve",
  sigma    = "Noise SD (σ)"
)

summary_long <- rbind(
  data.frame(factor = "n_curves", x = summary_table$n_curves[1:3],
             rmse_g0 = summary_table$rmse_g0[1:3], rmse_g1 = summary_table$rmse_g1[1:3],
             cov_g0  = summary_table$cov_g0[1:3],  cov_g1  = summary_table$cov_g1[1:3],
             bias_g0 = summary_table$bias_g0[1:3], bias_g1 = summary_table$bias_g1[1:3]),
  data.frame(factor = "n_time", x = summary_table$n_time[4:6],
             rmse_g0 = summary_table$rmse_g0[4:6], rmse_g1 = summary_table$rmse_g1[4:6],
             cov_g0  = summary_table$cov_g0[4:6],  cov_g1  = summary_table$cov_g1[4:6],
             bias_g0 = summary_table$bias_g0[4:6], bias_g1 = summary_table$bias_g1[4:6]),
  data.frame(factor = "sigma", x = summary_table$sigma[7:9],
             rmse_g0 = summary_table$rmse_g0[7:9], rmse_g1 = summary_table$rmse_g1[7:9],
             cov_g0  = summary_table$cov_g0[7:9],  cov_g1  = summary_table$cov_g1[7:9],
             bias_g0 = summary_table$bias_g0[7:9], bias_g1 = summary_table$bias_g1[7:9])
)
summary_long$factor_label <- group_labels[summary_long$factor]

rmse_long <- rbind(
  data.frame(summary_long[, c("factor", "factor_label", "x")],
             param = "γ0 (intercept)", rmse = summary_long$rmse_g0),
  data.frame(summary_long[, c("factor", "factor_label", "x")],
             param = "γ1 (slope)",     rmse = summary_long$rmse_g1)
)
rmse_long$factor_label <- factor(rmse_long$factor_label, levels = group_labels)

p_rmse <- ggplot(rmse_long, aes(x = x, y = rmse, colour = param, shape = param)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  facet_wrap(~ factor_label, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(col_g0, col_g1), name = NULL) +
  scale_shape_manual(values = c(19, 17), name = NULL) +
  labs(x = NULL, y = "RMSE", title = "RMSE of boundary coefficient estimates") +
  theme_sim

cov_long <- rbind(
  data.frame(summary_long[, c("factor", "factor_label", "x")],
             param = "γ0 (intercept)", coverage = 100 * summary_long$cov_g0),
  data.frame(summary_long[, c("factor", "factor_label", "x")],
             param = "γ1 (slope)",     coverage = 100 * summary_long$cov_g1)
)
cov_long$factor_label <- factor(cov_long$factor_label, levels = group_labels)

p_cov <- ggplot(cov_long, aes(x = x, y = coverage, colour = param, shape = param)) +
  geom_hline(yintercept = 95, linetype = "dashed", colour = col_true, linewidth = 0.7) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  facet_wrap(~ factor_label, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(col_g0, col_g1), name = NULL) +
  scale_shape_manual(values = c(19, 17), name = NULL) +
  scale_y_continuous(limits = c(50, 105)) +
  labs(x = NULL, y = "Coverage (%)", title = "Empirical coverage of 95% Wald CIs",
       caption = "Dashed line: nominal 95%") +
  theme_sim

p_rmse / p_cov 

# Sampling distributions — baseline scenario
baseline_res <- sim_results[["nc=10"]]

dist_data <- rbind(
  data.frame(param = "gamma[0]~(intercept)", value = baseline_res$gamma0_hat,
             true_val = true_gamma["gamma0"]),
  data.frame(param = "gamma[1]~(slope)",     value = baseline_res$gamma1_hat,
             true_val = true_gamma["gamma1"]),
  data.frame(param = "sigma^2~(noise)",      value = baseline_res$sigma2_hat,
             true_val = true_sigma2)
)
dist_data$param <- factor(dist_data$param,
                           levels = c("gamma[0]~(intercept)", "gamma[1]~(slope)",
                                      "sigma^2~(noise)"))

panel_summary <- data.frame(
  param    = levels(dist_data$param),
  mean_val = c(mean(baseline_res$gamma0_hat), mean(baseline_res$gamma1_hat),
               mean(baseline_res$sigma2_hat)),
  true_val = c(true_gamma["gamma0"], true_gamma["gamma1"], true_sigma2),
  col      = c(col_g0, col_g1, "#10B981")
)
panel_summary$param <- factor(panel_summary$param, levels = levels(dist_data$param))

ggplot(dist_data, aes(x = value)) +
  geom_density(aes(fill = param, colour = param), adjust = 1.2, alpha = 0.2, linewidth = 1) +
  geom_vline(data = panel_summary, aes(xintercept = true_val),
             colour = col_true, linetype = "dashed", linewidth = 0.9) +
  geom_vline(data = panel_summary, aes(xintercept = mean_val, colour = param),
             linetype = "dotted", linewidth = 0.9) +
  facet_wrap(~ param, scales = "free", labeller = label_parsed) +
  scale_fill_manual(values  = c(col_g0, col_g1, "#10B981"), guide = "none") +
  scale_colour_manual(values = c(col_g0, col_g1, "#10B981"), guide = "none") +
  labs(x = "Estimate", y = "Density",
       title = "Sampling distributions — baseline scenario (n_curves = 10, n_time = 20)",
       caption = "Dashed red: true value  ·  Dotted: empirical mean") +
  theme_sim +
  theme(legend.position = "none")

# CI strip chart — baseline scenario
n_strip   <- min(nrow(baseline_res), 60L)
strip_res <- baseline_res[seq_len(n_strip), ]
strip_res$rep <- seq_len(n_strip)

ci_long <- rbind(
  data.frame(rep = strip_res$rep, param = "gamma[0]~(intercept)",
             est = strip_res$gamma0_hat, lo = strip_res$gamma0_lo, hi = strip_res$gamma0_hi,
             hit = strip_res$gamma0_cov, true_v = true_gamma["gamma0"]),
  data.frame(rep = strip_res$rep, param = "gamma[1]~(slope)",
             est = strip_res$gamma1_hat, lo = strip_res$gamma1_lo, hi = strip_res$gamma1_hi,
             hit = strip_res$gamma1_cov, true_v = true_gamma["gamma1"])
)
ci_long$param <- factor(ci_long$param,
                         levels = c("gamma[0]~(intercept)", "gamma[1]~(slope)"))

cov_labels <- data.frame(
  param = factor(c("gamma[0]~(intercept)", "gamma[1]~(slope)"),
                 levels = c("gamma[0]~(intercept)", "gamma[1]~(slope)")),
  label = paste0("Coverage: ",
                 round(100 * c(mean(baseline_res$gamma0_cov),
                               mean(baseline_res$gamma1_cov)), 1L), "%"),
  x     = n_strip * 0.85,
  y     = c(max(strip_res$gamma0_hi) * 0.97, max(strip_res$gamma1_hi) * 0.97)
)

ggplot(ci_long, aes(x = rep)) +
  geom_linerange(aes(ymin = lo, ymax = hi, colour = hit), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(y = est, colour = hit), size = 0.9, shape = 20) +
  geom_hline(aes(yintercept = true_v), colour = col_true, linetype = "dashed", linewidth = 0.8) +
  geom_text(data = cov_labels, aes(x = x, y = y, label = label), size = 3.2, colour = "grey30") +
  facet_wrap(~ param, scales = "free_y", labeller = label_parsed) +
  scale_colour_manual(values = c("TRUE" = col_hit, "FALSE" = col_miss),
                      labels = c("TRUE" = "Covers truth", "FALSE" = "Misses"), name = NULL) +
  labs(x = "Replication", y = "Estimate",
       title = "95% Wald confidence intervals — first 60 replications",
       caption = "Dashed red: true value") +
  theme_sim

# Bias across all scenarios
bias_long <- rbind(
  data.frame(label = summary_table$label,
             group = rep(c("n_curves", "n_time", "sigma"), each = 3L),
             param = "γ0 (intercept)", bias = summary_table$bias_g0),
  data.frame(label = summary_table$label,
             group = rep(c("n_curves", "n_time", "sigma"), each = 3L),
             param = "γ1 (slope)",     bias = summary_table$bias_g1)
)
bias_long$label <- factor(bias_long$label, levels = summary_table$label)
bias_long$group <- factor(bias_long$group,
                           levels = c("n_curves", "n_time", "sigma"),
                           labels = group_labels)

ggplot(bias_long, aes(x = label, y = bias, colour = param, shape = param, group = param)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", colour = "grey70") +
  geom_line(linewidth = 0.8) + geom_point(size = 2.8) +
  scale_colour_manual(values = c(col_g0, col_g1), name = NULL) +
  scale_shape_manual(values = c(19, 17), name = NULL) +
  labs(x = "Scenario", y = "Bias",
       title = "Bias of boundary coefficient estimates across all scenarios") +
  theme_sim +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# Numerical summary
round(summary_table[, c("label", "n_curves", "n_time", "sigma", "n_conv",
                         "bias_g0", "rmse_g0", "cov_g0",
                         "bias_g1", "rmse_g1", "cov_g1",
                         "bias_sig2", "rmse_sig2")],
      digits = 3L)

