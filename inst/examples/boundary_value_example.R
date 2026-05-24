make_boundary_value_example <- function(
  n_curves = 4,
  n_time = 16,
  boundary_x = seq(-1, 1, length.out = n_curves),
  boundary_coef = c(0.5, 0.25),
  left_boundary = c(1, 0),
  right_boundary = c(1, 0),
  sigma = 0.15,
  seed = 1
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  left_boundary  <- fdaMixedRTMB::fixed_left_boundary_k1(as.numeric(left_boundary))
  right_boundary <- fdaMixedRTMB::right_boundary_operator_k1(as.numeric(right_boundary))

  time_grid <- seq(0, 1, length.out = n_time)
  curve <- factor(rep(seq_len(n_curves), each = n_time))
  time <- rep(time_grid, n_curves)

  boundary_X <- cbind("(Intercept)" = 1, boundary_x = boundary_x)
  boundary_coef <- rep_len(boundary_coef, ncol(boundary_X))
  curve_rho <- drop(boundary_X %*% boundary_coef)

  basis <- local({
    a <- min(time_grid)
    b <- max(time_grid)
    alpha_a <- left_boundary$boundary_row
    alpha_b <- right_boundary$boundary_row
    system_matrix <- rbind(
      c(alpha_a[1], alpha_a[1] * a + alpha_a[2]),
      c(alpha_b[1], alpha_b[1] * b + alpha_b[2])
    )
    coef <- solve(system_matrix, c(0, 1))
    drop(coef[1] + coef[2] * time_grid)
  })

  boundary_effect <- rep(basis, n_curves) * rep(curve_rho, each = n_time)
  serial_effect <- numeric(length(time))
  for (ii in seq_len(n_curves)) {
    rows <- curve == levels(curve)[ii]
    amp <- stats::rnorm(1, sd = 0.2)
    serial_effect[rows] <- amp * sin(pi * time_grid)
  }
  noise <- stats::rnorm(length(time), sd = sigma)
  y <- boundary_effect + serial_effect + noise

  data <- data.frame(
    y = y,
    time = time,
    curve = curve,
    boundary_x = rep(boundary_x, each = n_time)
  )

  list(
    data = data,
    y = y,
    time = time,
    curve = curve,
    boundary_coef = stats::setNames(boundary_coef, colnames(boundary_X)),
    curve_rho = stats::setNames(curve_rho, levels(curve)),
    boundary_effect = boundary_effect,
    serial_effect = serial_effect,
    left_boundary = left_boundary,
    right_boundary = right_boundary
  )
}

## Example:
## example_data <- make_boundary_value_example()
## fit <- fdaMixedRTMB::fdaLm_rtmb(
##   y | curve ~ 0,
##   data = example_data$data,
##   boundary_formula = ~ 1 + boundary_x,
##   left_boundary  = c(1, 0),
##   right_boundary = c(1, 0),
##   operator = "laplace"
## )
## summary(fit)
