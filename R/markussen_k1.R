# Active RTMB implementation for the K.order = 1 boundary-value model
# using the Markussen (2013) operator-projection / trace-integral REML approach.
#
# Engine: "markussen_reml"
#
# Active fitting path:
#   fdaLm_rtmb()
#     -> .fit_boundary_value_rtmb()
#       -> .fit_boundary_value_markussen_reml_k1()   [this file]
#             -> .markussen_evaluate_k1()            [markussen_projection_legacy.R]
#
# This implementation applies RTMB automatic differentiation to the Markussen
# REML objective.  Specifically:
#
#   - sigma^2 and ALL fixed effects (ordinary beta AND boundary coefficients
#     gamma) are profiled out analytically inside .markussen_evaluate_k1() via
#     forward-backward Green's-function projection sweeps.
#   - The REML log-determinant correction is approximated via fixed
#     Gauss-Legendre quadrature over the Markussen trace integral.
#   - Grouped random effects u are optionally marginalized by an RTMB Laplace
#     approximation (random = "u").  No latent serial-coefficient block s is
#     introduced; the serial effect X_m(t) is recovered via the projection
#     equations.
#   - RTMB differentiates the outer REML objective with respect to the variance
#     parameters (log_lambda_d2, optionally log_lambda_level),
#     and nlminb() optimizes using those gradients.
#
# Supported operators (K.order = 1):
#   Pure Laplace:    L f = -lambda_d2 * f''
#   Shifted Laplace: L f = -lambda_d2 * f'' + lambda_level * f
#
# The boundary-value extension enters as ordinary fixed-effect columns in X:
#   h(t) * boundary_X[m, ]  (the affine carrier satisfying the BCs).
# These columns are profiled out together with ordinary beta.
#
# Projection and trace-integral utility functions are in
# R/markussen_projection_legacy.R  (now the active Markussen engine utilities).


# main fitting function

.fit_boundary_value_markussen_reml_k1 <- function(
  prepared,
  operator,
  left_boundary,
  right_boundary,
  lambda_d2_init,
  sigma2_init,      # retained for API compatibility; sigma2 is profiled internally
  lambda_level_init,
  control,
  trace_quad_n,
  call = NULL
) {
  # basic setup
  estimate_level <- isTRUE(operator$estimate_level)
  n_obs          <- length(prepared$y)
  n_fixed        <- ncol(prepared$X)   # base_X + boundary_value_design columns

  lambda_d2_init <- max(as.numeric(lambda_d2_init), 1e-8)

  # Basic identifiability checks
  if (n_fixed >= n_obs) {
    stop(
      "The Markussen REML engine requires at least one residual degree of freedom ",
      "(n_obs > ncol(X))."
    )
  }
  if (n_fixed > 0L && qr(prepared$X)$rank < n_fixed) {
    stop("The Markussen REML engine requires `X` to have full column rank.")
  }

  # Boundary matrices (plain R constants, captured by the RTMB closure)
  left_row   <- left_boundary$boundary_row   # numeric length-2
  right_rows <- .right_boundary_rows_k1(     # (n_curves x 2) numeric matrix
    right_boundary = right_boundary,
    n_curves       = prepared$n_curves
  )

  # Plain R matrices passed to the RTMB tape as data
  X_full <- unname(prepared$X)   # (n_obs  x n_fixed)
  y_vec  <- prepared$y           # (n_obs)

  # parameter list for RTMB
  # Only variance parameters are free; sigma2, beta, and gamma are profiled.
  par_list <- list(
    log_lambda_d2 = log(lambda_d2_init)
  )
  if (estimate_level) {
    par_list$log_lambda_level <- log(max(lambda_level_init, 1e-8))
  }

  # Data list for RTMB
  data_list <- list(
    y_data = y_vec,
    X_data = X_full
  )

  # RTMB objective
  # Captures (as R constants at tape time):
  #   prepared, left_row, right_rows, estimate_level, trace_quad_n
  objective <- function(parms) {
    RTMB::getAll(data_list, parms, warn = FALSE)

    lambda_d2    <- exp(log_lambda_d2)
    lambda_level <- if (estimate_level) exp(log_lambda_level) else 0

    result <- .markussen_evaluate_k1(
      lambda_d2    = lambda_d2,
      lambda_level = lambda_level,
      y            = y_data,
      X            = X_data,
      prepared     = prepared,
      left_row     = left_row,
      right_rows   = right_rows,
      trace_quad_n = as.integer(trace_quad_n)
    )

    result$nll
  }

  obj <- RTMB::MakeADFun(
    objective,
    par_list,
    silent = TRUE
  )

  opt <- do.call(
    stats::nlminb,
    c(list(start = obj$par, objective = obj$fn, gradient = obj$gr), control)
  )
  obj$par    <- opt$par
  fitted_par <- obj$env$parList(opt$par)

  # extract var estimates
  lambda_d2_hat    <- exp(as.numeric(fitted_par$log_lambda_d2))
  lambda_level_hat <- if (estimate_level) exp(as.numeric(fitted_par$log_lambda_level)) else 0

  # Re-run at the optimum to recover sigma2, Cbeta, serial effect, etc.
  final <- .markussen_evaluate_k1(
    lambda_d2    = lambda_d2_hat,
    lambda_level = lambda_level_hat,
    y            = y_vec,
    X            = X_full,
    prepared     = prepared,
    left_row     = left_row,
    right_rows   = right_rows,
    trace_quad_n = as.integer(trace_quad_n)
  )

  # split profiled fixed effects into beta and gamma
  base_cols     <- prepared$fixed_column_groups$beta
  boundary_cols <- prepared$fixed_column_groups$boundary

  beta_hat_full <- as.numeric(final$beta)   # length n_fixed (or 0)
  beta_hat      <- if (length(base_cols)) beta_hat_full[base_cols] else numeric(0L)
  gamma_hat     <- if (length(boundary_cols)) beta_hat_full[boundary_cols] else numeric(0L)
  names(beta_hat)  <- colnames(prepared$base_X)
  names(gamma_hat) <- colnames(prepared$boundary_value_design)

  beta_full <- c(beta_hat, gamma_hat)
  names(beta_full) <- colnames(prepared$X)

  # fixed-effect covariance
  # Var(beta_hat) = sigma2 * Cbeta  (analytically from the Markussen projection)
  Cbeta_mat <- as.matrix(final$Cbeta)
  if (ncol(Cbeta_mat) > 0L) {
    colnames(Cbeta_mat) <- rownames(Cbeta_mat) <- colnames(prepared$X)
  }
  fixed_vcov <- final$sigma2 * Cbeta_mat

  boundary_coef_vcov <- if (length(boundary_cols)) fixed_vcov[boundary_cols, boundary_cols, drop = FALSE] else matrix(numeric(0L), 0L, 0L)

  # fitted-value decomposition
  # fitted = fixed_effect + boundary_effect + serial_effect
  # All are recoverable from the Markussen evaluate output.
  fixed_effect    <- if (length(base_cols) > 0L) drop(prepared$base_X %*% beta_hat) else rep(0, n_obs)
  boundary_effect <- if (length(boundary_cols) > 0L) drop(prepared$boundary_value_design %*% gamma_hat) else rep(0, n_obs)

  # curve boundary values
  curve_rho        <- .resolve_curve_boundary_values_k1(prepared$boundary_X, gamma_hat)
  names(curve_rho) <- prepared$curve_levels

  # right-boundary quantities for the serial effect 
  right_boundary_deriv <- as.numeric(final$right_boundary_derivative)
  right_serial_value   <- vapply(
    prepared$curve_index,
    function(idx) final$serial_effect[idx[length(idx)]],
    numeric(1L)
  )
  right_serial_operator <- right_boundary$alpha[[1]] * right_serial_value +
    right_boundary$alpha[[2]] * right_boundary_deriv
  names(right_serial_value)    <- prepared$curve_levels
  names(right_serial_operator) <- prepared$curve_levels
  names(right_boundary_deriv)  <- prepared$curve_levels

  # assemble result object 
  fitted_object <- list(
    call                           = if (is.null(call)) match.call() else call,
    method                         = "REML",
    engine                         = "markussen_reml",
    boundary_value                 = length(boundary_cols) > 0L,
    beta                           = beta_hat,
    beta_full                      = beta_full,
    gamma                          = gamma_hat,
    boundary_coef                  = gamma_hat,
    boundary_coef_vcov             = boundary_coef_vcov,
    fixed_vcov                     = fixed_vcov,
    Cbeta                          = Cbeta_mat,
    tau                            = final$sigma2 / lambda_d2_hat,
    sigma2                         = final$sigma2,
    lambda_d2                      = lambda_d2_hat,
    lambda_level                   = lambda_level_hat,
    trace_integral                 = as.numeric(final$trace_integral),
    curve_rho                      = curve_rho,
    right_boundary_rows            = right_rows,
    right_boundary_derivative      = right_boundary_deriv,
    right_boundary_serial_value    = right_serial_value,
    right_boundary_serial_operator = right_serial_operator,
    fitted                         = as.numeric(final$fitted),
    fixed_effect                   = fixed_effect,
    boundary_effect                = boundary_effect,
    serial_effect                  = as.numeric(final$serial_effect),
    serial_derivative              = as.numeric(final$serial_derivative),
    data                           = prepared,
    operator                       = operator,
    left_boundary                  = left_boundary,
    right_boundary                 = right_boundary,
    objective                      = obj,
    opt                            = opt,
    logLik                         = -opt$objective,
    logLik_type                    = "REML",
    converged                      = identical(opt$convergence, 0L)
  )

  class(fitted_object) <- "boundary_value_fit"
  fitted_object
}
