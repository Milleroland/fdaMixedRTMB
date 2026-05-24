#' Boundary operators for the k = 1 model
#'
#' The model imposes a homogeneous linear boundary condition on the serial
#' effect at the left endpoint:
#' \deqn{\alpha_{a,0}\, X(a) + \alpha_{a,1}\, X'(a) = 0.}{alpha_a0 * X(a) + alpha_a1 * X'(a) = 0.}
#'
#' The right endpoint uses the same linear operator row, but its value is
#' estimated:
#' \deqn{\alpha_{b,0}\, X(b) + \alpha_{b,1}\, X'(b) = \rho.}{alpha_b0 * X(b) + alpha_b1 * X'(b) = rho.}
#'
#' Both boundary operators are specified as a numeric length-two vector
#' `c(alpha0, alpha1)` or a one-row two-column matrix. A Dirichlet (value)
#' condition corresponds to `c(1, 0)` and a Neumann (derivative) condition
#' to `c(0, 1)`.
#'
#' @param alpha Numeric length-two vector `c(alpha0, alpha1)` or a one-row
#'   two-column matrix defining the boundary operator row.
#' @return A boundary-operator object.
#' @export
fixed_left_boundary_k1 <- function(alpha = c(1, 0)) {
  alpha <- .coerce_boundary_alpha_k1(alpha, what = "alpha")
  names(alpha) <- c("alpha0", "alpha1")
  type <- .classify_left_boundary_type_k1(alpha)

  structure(
    list(
      side = "left",
      type = type,
      alpha = alpha,
      boundary_row = unname(alpha),
      value = 0,
      description = sprintf(
        "Fixed left operator %.3g * X(a) + %.3g * X'(a) = 0.",
        alpha[[1]], alpha[[2]]
      )
    ),
    class = "left_boundary_k1"
  )
}

.classify_left_boundary_type_k1 <- function(alpha) {
  eps <- sqrt(.Machine$double.eps)
  if (abs(alpha[[1]]) <= eps && abs(alpha[[2]]) > eps) {
    return("neumann")
  }
  if (abs(alpha[[2]]) <= eps && abs(alpha[[1]]) > eps) {
    return("dirichlet")
  }
  "robin"
}

#' @rdname fixed_left_boundary_k1
#' @return An object of class `"right_boundary_operator_k1"`.
#' @export
right_boundary_operator_k1 <- function(alpha = c(1, 0)) {
  alpha <- .coerce_boundary_alpha_k1(alpha, what = "alpha")
  names(alpha) <- c("alpha0", "alpha1")
  structure(
    list(
      side = "right",
      type = .classify_left_boundary_type_k1(alpha),
      alpha = alpha,
      boundary_row = unname(alpha),
      description = "Fixed right operator alpha0 * X(b) + alpha1 * X'(b) = rho."
    ),
    class = "right_boundary_operator_k1"
  )
}

.is_right_boundary_operator_k1 <- function(x) {
  inherits(x, "right_boundary_operator_k1")
}

.is_left_boundary_k1 <- function(x) {
  inherits(x, "left_boundary_k1")
}

.coerce_boundary_alpha_k1 <- function(alpha, what) {
  if (is.matrix(alpha)) {
    if (!identical(dim(alpha), c(1L, 2L))) {
      stop("`", what, "` matrix must have dimension 1 x 2 for K.order = 1.")
    }
    alpha <- as.numeric(alpha[1L, ])
  } else {
    alpha <- as.numeric(alpha)
  }
  if (length(alpha) != 2L) {
    stop("`", what, "` must be a numeric vector of length two.")
  }
  if (anyNA(alpha) || any(!is.finite(alpha))) {
    stop("`", what, "` must contain finite numeric values.")
  }
  if (all(abs(alpha) <= sqrt(.Machine$double.eps))) {
    stop("`", what, "` must contain at least one nonzero value.")
  }
  alpha
}

.coerce_left_boundary_k1 <- function(x) {
  if (.is_left_boundary_k1(x)) {
    return(x)
  }
  if (!is.numeric(x) && !is.matrix(x)) {
    stop(
      "`left_boundary` must be a numeric length-two vector or a one-row ",
      "two-column matrix, or a `fixed_left_boundary_k1()` object."
    )
  }
  fixed_left_boundary_k1(x)
}

.coerce_right_boundary_operator_k1 <- function(x) {
  if (.is_right_boundary_operator_k1(x)) {
    return(x)
  }
  if (!is.numeric(x) && !is.matrix(x)) {
    stop(
      "`right_boundary` must be a numeric length-two vector or a one-row ",
      "two-column matrix, or a `right_boundary_operator_k1()` object."
    )
  }
  right_boundary_operator_k1(x)
}

.validate_right_boundary_operator_k1 <- function(x) {
  if (!.is_right_boundary_operator_k1(x)) {
    stop(
      "`right_boundary` must be created by `right_boundary_operator_k1()` ",
      "or use a supported boundary specification."
    )
  }
  invisible(x)
}

.validate_left_boundary_k1 <- function(x) {
  if (!.is_left_boundary_k1(x)) {
    stop(
      "`left_boundary` must be created by `fixed_left_boundary_k1()` ",
      "or use a supported boundary specification."
    )
  }
  invisible(x)
}

.right_boundary_label_k1 <- function(right_boundary) {
  alpha <- right_boundary$alpha
  sprintf(
    "%.3g * X(b) + %.3g * X'(b) = rho",
    alpha[[1]],
    alpha[[2]]
  )
}

.left_boundary_label_k1 <- function(left_boundary) {
  alpha <- left_boundary$alpha
  sprintf(
    "%.3g * X(a) + %.3g * X'(a) = 0",
    alpha[[1]],
    alpha[[2]]
  )
}

.right_boundary_rows_k1 <- function(right_boundary, n_curves = 1L) {
  .validate_right_boundary_operator_k1(right_boundary)
  out <- matrix(right_boundary$boundary_row, nrow = n_curves, ncol = 2L, byrow = TRUE)
  rownames(out) <- NULL
  colnames(out) <- c("X_b", "Xprime_b")
  out
}

#' Compute the affine carrier basis h(t) for the right-boundary value
#'
#' Constructs the fixed-effect basis `h(t)` that carries the curve-specific
#' right-boundary value `rho`. The basis is chosen as an affine function
#' `h(t) = c0 + c1 * t` (the null space of the pure Laplace `L = -d^2/dt^2`)
#' satisfying the homogeneous left boundary condition
#' `alpha_{a,0} * h(a) + alpha_{a,1} * h'(a) = 0` and the unit right boundary
#' value `alpha_{b,0} * h(b) + alpha_{b,1} * h'(b) = 1`. The two conditions
#' assemble into a 2x2 linear system in `(c0, c1)`.
#'
#' Note: for the pure Laplace operator this is the null-space basis. For the
#' shifted operator, whose homogeneous solutions are exponential, the same
#' affine basis is retained as a convenient fixed-effect carrier because it
#' satisfies the required boundary constraints regardless of the fitted level
#' penalty. The system is singular if (and only if) the two boundary operators
#' are linearly dependent on the affine span.
#'
#' @keywords internal
#' @noRd
.boundary_value_basis_k1 <- function(time, right_boundary, left_boundary) {
  .validate_right_boundary_operator_k1(right_boundary)
  .validate_left_boundary_k1(left_boundary)

  time <- as.numeric(time)
  a <- min(time)
  b <- max(time)
  alpha_a <- left_boundary$boundary_row
  alpha_b <- right_boundary$boundary_row

  Mmat <- rbind(
    c(alpha_a[1], alpha_a[1] * a + alpha_a[2]),
    c(alpha_b[1], alpha_b[1] * b + alpha_b[2])
  )
  rhs <- c(0, 1)

  det_M <- Mmat[1, 1] * Mmat[2, 2] - Mmat[1, 2] * Mmat[2, 1]
  scale <- 1 + max(abs(c(alpha_a, alpha_b, a, b)))
  if (!is.finite(det_M) || abs(det_M) <= sqrt(.Machine$double.eps) * scale) {
    stop(
      "The supplied left and right boundary operators are linearly dependent ",
      "on the affine span (singular system): cannot construct a unique ",
      "carrier basis h(t) for the right-boundary value. Choose boundary rows ",
      "that are linearly independent on the span of (1, t)."
    )
  }

  coefs <- as.numeric(solve(Mmat, rhs))
  c0 <- coefs[1]
  c1 <- coefs[2]

  list(
    values = c0 + c1 * time,
    coefficients = c(c0 = c0, c1 = c1),
    denominator = det_M
  )
}

.build_boundary_value_fixed_part_k1 <- function(prepared, right_boundary, left_boundary) {
  basis <- .boundary_value_basis_k1(
    time = prepared$time_grid,
    right_boundary = right_boundary,
    left_boundary = left_boundary
  )

  n_boundary <- ncol(prepared$boundary_X)
  design <- matrix(0, nrow = length(prepared$y), ncol = n_boundary)
  for (ii in seq_len(prepared$n_curves)) {
    idx <- prepared$curve_index[[ii]]
    design[idx, ] <- outer(basis$values, prepared$boundary_X[ii, ], FUN = "*")
  }

  boundary_names <- colnames(prepared$boundary_X)
  if (is.null(boundary_names)) {
    boundary_names <- paste0("b", seq_len(n_boundary))
  }
  colnames(design) <- paste0("rho_", boundary_names)

  list(
    design = design,
    basis = basis$values,
    coefficients = basis$coefficients,
    denominator = basis$denominator
  )
}

.resolve_curve_boundary_values_k1 <- function(boundary_X, boundary_coef) {
  if (!length(boundary_coef)) {
    return(rep(0, nrow(boundary_X)))
  }
  drop(boundary_X %*% boundary_coef)
}
