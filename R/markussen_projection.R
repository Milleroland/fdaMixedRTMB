# Markussen operator-projection engine — ACTIVE utility functions
#
# The functions in this file implement the Markussen (2013) operator-projection
# / trace-integral REML approach for K.order = 1 functional mixed models.
# They are called by the active fitting path:
#
#   fdaLm_rtmb()
#     -> .fit_boundary_value_rtmb()
#       -> .fit_boundary_value_markussen_reml_k1()  [markussen_k1.R]
#             -> .markussen_evaluate_k1()           [this file]
#
# The active engine (engine = "markussen_reml"):
#   - profiles out sigma^2 and ALL fixed effects (ordinary beta and boundary
#     coefficients gamma) analytically via forward/backward Green's-function
#     projection sweeps
#   - approximates the REML log-determinant via fixed Gauss-Legendre quadrature
#     of the Markussen trace integral, matching the fdaMixed trace machinery
#   - optionally marginalizes grouped random effects u via an RTMB Laplace
#     approximation; no latent serial-coefficient block s is introduced
#   - RTMB automatic differentiation drives nlminb() over the remaining
#     variance parameters (lambda_d2, optionally lambda_level, sigma_u2)
#
# The helper functions are:
#   .markussen_state_vector_k1   — builds the [1; eta] state vector
#   .markussen_roots_k1          — characteristic roots of the discrete Green's fn
#   .markussen_project_matrix_k1 — forward/backward sweep projecting a matrix
#   .markussen_trace_single_k1   — trace of Green's fn at one quadrature node
#   .markussen_trace_integral_k1 — full Gauss-Legendre trace-integral
#   .markussen_evaluate_k1       — full REML evaluation (profiled sigma2 + beta)

.markussen_state_vector_k1 <- function(eta) {
  RTMB::matrix(c(1 + eta * 0, eta), ncol = 1L)
}

.markussen_roots_k1 <- function(constant_term, quadratic_term) {
  root_scale <- sqrt(-constant_term / quadratic_term)
  list(left = -root_scale, right = root_scale)
}

.markussen_project_matrix_k1 <- function(
  Ymat,
  eta_left,
  eta_right,
  Fleft,
  Fright,
  left,
  right,
  alpha_2k
) {
  NN <- nrow(Ymat)
  Ycols <- ncol(Ymat)
  Delta <- (right - left) / NN

  vleft <- -1 / (eta_right - eta_left)
  vright <- 1 / (eta_right - eta_left)

  xi_left <- (exp(Delta / 2 * eta_left) - 1) / eta_left
  xi0_left <- (1 - (1 - Delta * eta_left) * exp(Delta * eta_left)) /
    (Delta * eta_left^2)
  xi1_left <- (exp(Delta * eta_left) - 1 - Delta * eta_left) /
    (Delta * eta_left^2)

  xi_right <- (1 - exp(-Delta / 2 * eta_right)) / eta_right
  xi0_right <- (exp(-Delta * eta_right) - 1 + Delta * eta_right) /
    (Delta * eta_right^2)
  xi1_right <- (1 - (1 + Delta * eta_right) * exp(-Delta * eta_right)) /
    (Delta * eta_right^2)

  Fleft_W <- drop(Fleft %*% .markussen_state_vector_k1(eta_right)) /
    drop(Fleft %*% .markussen_state_vector_k1(eta_left))
  Fright_W <- drop(Fright %*% .markussen_state_vector_k1(eta_left)) /
    drop(Fright %*% .markussen_state_vector_k1(eta_right))

  discounter_left <- exp(Delta * eta_left)
  half_discounter_left <- exp(0.5 * Delta * eta_left)
  discounter_right <- exp(-Delta * eta_right)
  half_discounter_right <- exp(-0.5 * Delta * eta_right)

  zero_row <- RTMB::matrix(eta_left * 0, nrow = 1L, ncol = Ycols)
  sum1 <- zero_row
  sum2 <- (vleft * xi_left) * RTMB::matrix(Ymat[1, ], nrow = 1L)
  sum3 <- zero_row
  sum4 <- (vright * xi_right) * RTMB::matrix(Ymat[1, ], nrow = 1L)

  prod_discounter_left <- half_discounter_left
  prod_discounter_right <- half_discounter_right

  exp_FW_exp <- exp((0.5 - NN) * Delta * eta_right) *
    Fright_W *
    exp((NN - 0.5) * Delta * eta_left)
  exp_FW <- prod_discounter_left * Fleft_W
  phi_fac <- 1 / (1 - exp_FW * prod_discounter_right * exp_FW_exp)
  phi0 <- (1 - exp_FW_exp) * phi_fac
  phi1 <- (eta_left - eta_right * exp_FW_exp) * phi_fac

  proj <- RTMB::matrix(eta_left * 0, nrow = NN, ncol = 2L * Ycols)
  row_forward <- drop(sum2 + exp_FW * sum4)
  proj[1, ] <- c(phi0 * row_forward, phi1 * row_forward)

  if (NN > 1L) {
    for (nn in 2:NN) {
      sum1 <- discounter_left * sum1 +
        (vleft * xi0_left) * RTMB::matrix(Ymat[nn - 1L, ], nrow = 1L)
      sum2 <- discounter_left * sum2 +
        (vleft * xi1_left) * RTMB::matrix(Ymat[nn, ], nrow = 1L)
      sum3 <- sum3 +
        (prod_discounter_right * vright * xi0_right) *
          RTMB::matrix(Ymat[nn - 1L, ], nrow = 1L)
      sum4 <- sum4 +
        (prod_discounter_right * vright * xi1_right) *
          RTMB::matrix(Ymat[nn, ], nrow = 1L)

      prod_discounter_left <- prod_discounter_left * discounter_left
      prod_discounter_right <- prod_discounter_right * discounter_right

      exp_FW_exp <- exp((0.5 + nn - 1L - NN) * Delta * eta_right) *
        Fright_W *
        exp((NN - (nn - 1L) - 0.5) * Delta * eta_left)
      exp_FW <- prod_discounter_left * Fleft_W
      phi_fac <- 1 / (1 - exp_FW * prod_discounter_right * exp_FW_exp)
      phi0 <- (1 - exp_FW_exp) * phi_fac
      phi1 <- (eta_left - eta_right * exp_FW_exp) * phi_fac

      row_forward <- drop(sum1 + sum2 + exp_FW * (sum3 + sum4))
      proj[nn, ] <- c(phi0 * row_forward, phi1 * row_forward)
    }
  }

  sum1 <- (vright * xi_right) * RTMB::matrix(Ymat[NN, ], nrow = 1L)
  sum2 <- zero_row
  sum3 <- (vleft * xi_left) * RTMB::matrix(Ymat[NN, ], nrow = 1L)
  sum4 <- zero_row

  prod_discounter_left <- half_discounter_left
  prod_discounter_right <- half_discounter_right

  exp_FW_exp <- exp((NN - 0.5) * Delta * eta_left) *
    Fleft_W *
    exp((0.5 - NN) * Delta * eta_right)
  exp_FW <- prod_discounter_right * Fright_W
  phi_fac <- 1 / (1 - exp_FW * prod_discounter_left * exp_FW_exp)
  phi0 <- (1 - exp_FW_exp) * phi_fac
  phi1 <- (eta_right - eta_left * exp_FW_exp) * phi_fac

  row_backward <- drop(sum1 + exp_FW * sum3)
  proj[NN, ] <- proj[NN, ] - c(phi0 * row_backward, phi1 * row_backward)

  if (NN > 1L) {
    for (nn in NN:2) {
      sum1 <- discounter_right * sum1 +
        (vright * xi0_right) * RTMB::matrix(Ymat[nn - 1L, ], nrow = 1L)
      sum2 <- discounter_right * sum2 +
        (vright * xi1_right) * RTMB::matrix(Ymat[nn, ], nrow = 1L)
      sum3 <- sum3 +
        (prod_discounter_left * vleft * xi0_left) *
          RTMB::matrix(Ymat[nn - 1L, ], nrow = 1L)
      sum4 <- sum4 +
        (prod_discounter_left * vleft * xi1_left) *
          RTMB::matrix(Ymat[nn, ], nrow = 1L)

      prod_discounter_left <- prod_discounter_left * discounter_left
      prod_discounter_right <- prod_discounter_right * discounter_right

      exp_FW_exp <- exp((nn - 0.5) * Delta * eta_left) *
        Fleft_W *
        exp((0.5 - nn) * Delta * eta_right)
      exp_FW <- prod_discounter_right * Fright_W
      phi_fac <- 1 / (1 - exp_FW * prod_discounter_left * exp_FW_exp)
      phi0 <- (1 - exp_FW_exp) * phi_fac
      phi1 <- (eta_right - eta_left * exp_FW_exp) * phi_fac

      row_backward <- drop(sum1 + sum2 + exp_FW * (sum3 + sum4))
      proj[nn - 1L, ] <- proj[nn - 1L, ] - c(phi0 * row_backward, phi1 * row_backward)
    }
  }

  proj / alpha_2k
}

.markussen_gauss_legendre_k1 <- function(n) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(nodes = 0.5, weights = 1))
  }

  ii <- seq_len(n - 1L)
  off_diag <- ii / sqrt(4 * ii^2 - 1)
  jacobi <- matrix(0, nrow = n, ncol = n)
  jacobi[cbind(ii, ii + 1L)] <- off_diag
  jacobi[cbind(ii + 1L, ii)] <- off_diag

  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  list(
    nodes = (eig$values[ord] + 1) / 2,
    weights = eig$vectors[1L, ord]^2
  )
}

.markussen_trace_single_k1 <- function(
  left,
  right,
  tau,
  eta_left,
  eta_right,
  Fleft,
  Fright,
  N
) {
  Delta <- (right - left) / N

  vleft <- -1 / (eta_right - eta_left)
  vright <- 1 / (eta_right - eta_left)

  exp_left <- exp((right - left) * eta_left)
  exp_right <- exp((left - right) * eta_right)

  Fleft_W <- drop(Fleft %*% .markussen_state_vector_k1(eta_right)) /
    drop(Fleft %*% .markussen_state_vector_k1(eta_left))
  Fright_W <- drop(Fright %*% .markussen_state_vector_k1(eta_left)) /
    drop(Fright %*% .markussen_state_vector_k1(eta_right))

  Aleft_left <- N * exp_left
  Aright_right <- N * exp_right
  Aleft_right <- (exp((right - left) * (eta_left - eta_right)) - 1) /
    (Delta * (eta_left - eta_right))
  Aright_left <- Aleft_right

  Bmat <- Fleft_W * exp_right * Fright_W /
    (1 - exp_left * Fleft_W * exp_right * Fright_W)

  trace_green <- (
    N * vleft +
      Fleft_W * Aleft_right * vright -
      Fright_W * Aright_left * vleft -
      (Fright_W * exp_left * Fleft_W) * Aright_right * vright +
      Bmat * Aleft_left * vleft +
      (Bmat * exp_left * Fleft_W) * Aleft_right * vright -
      (Fright_W * exp_left * Bmat) * Aright_left * vleft -
      (Fright_W * exp_left * Bmat * exp_left * Fleft_W) * Aright_right * vright
  )

  trace_green / tau
}

.markussen_trace_integral_k1 <- function(
  lambda_d2,
  lambda_level,
  left,
  right,
  Fleft,
  Fright,
  N,
  n_quad
) {
  Delta <- (right - left) / N
  quadratic_term <- -Delta * lambda_d2
  quad <- .markussen_gauss_legendre_k1(n_quad)
  trace_integral <- lambda_d2 * 0

  for (ii in seq_along(quad$nodes)) {
    roots <- .markussen_roots_k1(
      constant_term = quad$nodes[[ii]] + Delta * lambda_level,
      quadratic_term = quadratic_term
    )
    trace_integral <- trace_integral + quad$weights[[ii]] *
      .markussen_trace_single_k1(
        left = left,
        right = right,
        tau = -lambda_d2,
        eta_left = roots$left,
        eta_right = roots$right,
        Fleft = Fleft,
        Fright = Fright,
        N = N
      )
  }

  trace_integral
}

.markussen_evaluate_k1 <- function(
  lambda_d2,
  lambda_level = 0,
  y,
  X,
  Z,
  prepared,
  left_row,
  right_rows,
  sigma_u2 = 0,
  u = numeric(0),
  trace_quad_n = 64L
) {
  NN <- prepared$n_per_curve
  MM <- prepared$n_curves
  n_total <- length(y)
  p0 <- ncol(X)
  q <- ncol(Z)
  has_random <- q > 0L
  left <- min(prepared$time_grid)
  right <- max(prepared$time_grid)
  Delta <- (right - left) / NN
  alpha_2k <- -Delta * lambda_d2

  roots <- .markussen_roots_k1(
    constant_term = 1 + Delta * lambda_level,
    quadratic_term = alpha_2k
  )

  random_part <- if (has_random) drop(Z %*% u) else y * 0
  y_adj <- y - random_part

  proj_list <- vector("list", MM)
  deriv_list <- vector("list", MM)
  gamma_proj_list <- if (p0 > 0L) vector("list", MM) else NULL
  beta_hat <- numeric(0)
  Cbeta <- RTMB::matrix(0, 0, 0)

  if (p0 > 0L) {
    Cbeta_num <- RTMB::matrix(lambda_d2 * 0, nrow = p0, ncol = p0)
    Cbeta_rhs <- RTMB::matrix(lambda_d2 * 0, nrow = p0, ncol = 1L)
  }

  trace_integral <- lambda_d2 * 0
  for (ii in seq_len(MM)) {
    idx <- prepared$curve_index[[ii]]
    Fright <- RTMB::matrix(right_rows[ii, ], nrow = 1L)
    Fleft <- RTMB::matrix(left_row, nrow = 1L)
    trace_integral <- trace_integral + .markussen_trace_integral_k1(
      lambda_d2 = lambda_d2,
      lambda_level = lambda_level,
      left = left,
      right = right,
      Fleft = Fleft,
      Fright = Fright,
      N = NN,
      n_quad = trace_quad_n
    )

    Ymat_ii <- cbind(y_adj[idx], X[idx, , drop = FALSE])
    proj_ii <- .markussen_project_matrix_k1(
      Ymat = Ymat_ii,
      eta_left = roots$left,
      eta_right = roots$right,
      Fleft = Fleft,
      Fright = Fright,
      left = left,
      right = right,
      alpha_2k = alpha_2k
    )

    proj_list[[ii]] <- proj_ii
    deriv_list[[ii]] <- proj_ii[, ncol(Ymat_ii) + seq_len(ncol(Ymat_ii)), drop = FALSE]

    if (p0 > 0L) {
      gamma_ii <- X[idx, , drop = FALSE]
      gamma_proj_ii <- proj_ii[, 1L + seq_len(p0), drop = FALSE]
      gamma_proj_list[[ii]] <- gamma_proj_ii
      Cbeta_num <- Cbeta_num + crossprod(gamma_ii, gamma_ii - gamma_proj_ii)
      Cbeta_rhs <- Cbeta_rhs + crossprod(gamma_ii, y_adj[idx] - proj_ii[, 1])
    }
  }

  if (p0 > 0L) {
    Cbeta <- RTMB::solve(Cbeta_num)
    beta_hat <- drop(Cbeta %*% Cbeta_rhs)
  }

  x_hat <- RTMB::matrix(lambda_d2 * 0, nrow = n_total, ncol = 2L)
  mean_hat <- x_hat[, 1]
  right_boundary_derivative <- RTMB::matrix(lambda_d2 * 0, nrow = MM, ncol = 1L)[, 1]
  cond_resid <- x_hat[, 1]

  quad_resid <- lambda_d2 * 0
  serial_norm <- lambda_d2 * 0
  deriv_norm <- lambda_d2 * 0

  for (ii in seq_len(MM)) {
    idx <- prepared$curve_index[[ii]]
    x_curve_ii <- proj_list[[ii]][, 1]
    if (p0 > 0L) {
      x_curve_ii <- x_curve_ii - drop(gamma_proj_list[[ii]] %*% beta_hat)
      mean_hat[idx] <- drop(X[idx, , drop = FALSE] %*% beta_hat)
    }
    x_hat[idx, 1] <- x_curve_ii

    deriv_curve_ii <- deriv_list[[ii]][, 1]
    if (p0 > 0L) {
      gamma_deriv_proj_ii <- deriv_list[[ii]][, 1L + seq_len(p0), drop = FALSE]
      deriv_curve_ii <- deriv_curve_ii - drop(gamma_deriv_proj_ii %*% beta_hat)
    }
    x_hat[idx, 2] <- deriv_curve_ii

    cond_resid[idx] <- y_adj[idx] - mean_hat[idx] - x_curve_ii
    quad_resid <- quad_resid + sum(cond_resid[idx]^2)
    serial_norm <- serial_norm + Delta * sum(x_curve_ii^2)
    deriv_norm <- deriv_norm + Delta * sum(deriv_curve_ii^2)
    right_boundary_derivative[ii] <- deriv_curve_ii[length(idx)]
  }

  sigma2_hat <- (quad_resid + lambda_level * serial_norm + lambda_d2 * deriv_norm) /
    (n_total - p0)
  log_det_reml <- if (p0 > 0L) determinant(Cbeta_num, logarithm = TRUE)$modulus[[1L]] else lambda_d2 * 0

  nll <- (n_total - p0) * log(sigma2_hat) +
    trace_integral + log_det_reml + n_total - p0
  if (has_random) {
    nll <- nll + q * log(sigma_u2) + sum(u^2) / sigma_u2
  }

  list(
    nll = nll,
    sigma2 = sigma2_hat,
    lambda_d2 = lambda_d2,
    lambda_level = lambda_level,
    tau = sigma2_hat / lambda_d2,
    beta = beta_hat,
    u_hat = if (has_random) u else numeric(0),
    serial_effect = x_hat[, 1],
    serial_derivative = x_hat[, 2],
    random_effect = random_part,
    fitted = mean_hat + random_part + x_hat[, 1],
    right_boundary_derivative = right_boundary_derivative,
    trace_integral = trace_integral,
    Cbeta = Cbeta,
    cond_resid = cond_resid
  )
}
