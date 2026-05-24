if (getRversion() >= "2.15.1") {
  # Variables used inside RTMB objective closures (injected by RTMB::getAll).
  # The Markussen REML engine (engine = "markussen_reml") uses:
  #   - from parms:     log_lambda_d2, log_lambda_level, log_sigma_u2, u
  #   - from data_list: y_data, X_data, Z_data
  # sigma2, beta, gamma, and the serial-coefficient block s are all profiled
  # analytically and are therefore NOT RTMB parameters.
  utils::globalVariables(c(
    "log_lambda_d2", "log_lambda_level", "log_sigma_u2",
    "u",
    "y_data", "X_data", "Z_data"
  ))
}

.augment_boundary_value_design_k1 <- function(prepared, right_boundary, left_boundary) {
  base_X <- prepared$X
  n_base <- ncol(base_X)
  n_boundary <- ncol(prepared$boundary_X)

  if (n_boundary == 0L) {
    prepared$base_X <- base_X
    prepared$boundary_value_design <- matrix(
      numeric(0L),
      nrow = length(prepared$y),
      ncol = 0L
    )
    colnames(prepared$boundary_value_design) <- character(0L)
    prepared$boundary_value_basis <- numeric(0L)
    prepared$boundary_value_denominator <- NA_real_
    prepared$fixed_column_groups <- list(
      beta = if (n_base) seq_len(n_base) else integer(0L),
      boundary = integer(0L)
    )
    return(prepared)
  }

  boundary_part <- .build_boundary_value_fixed_part_k1(
    prepared = prepared,
    right_boundary = right_boundary,
    left_boundary = left_boundary
  )

  prepared$base_X <- base_X
  prepared$boundary_value_design <- boundary_part$design
  prepared$boundary_value_basis <- boundary_part$basis
  prepared$boundary_value_denominator <- boundary_part$denominator
  prepared$X <- cbind(base_X, boundary_part$design)

  prepared$fixed_column_groups <- list(
    beta = if (n_base) seq_len(n_base) else integer(0),
    boundary = if (n_boundary) n_base + seq_len(n_boundary) else integer(0)
  )

  prepared
}

# Internal matrix-based implementation used by fdaLm_rtmb(). Keep this
# unexported so users fit models through the fdaMixed-style formula interface.
.fit_boundary_value_rtmb <- function(
  y,
  time = NULL,
  curve = NULL,
  data = NULL,
  X = NULL,
  Z = NULL,
  boundary_X = NULL,
  operator = operator_k1(),
  left_boundary = c(1, 0),
  right_boundary = c(1, 0),
  tau_init = NULL,
  sigma2_init = NULL,
  sigma_u2_init = NULL,
  lambda_level_init = NULL,
  trace_quad_n = 64L,
  control = list(),
  boundary_formula = NULL
) {
  call <- match.call()
  caller <- parent.frame()
  y_missing <- missing(y)
  time_missing <- missing(time)
  curve_missing <- missing(curve)
  X_missing <- missing(X)
  Z_missing <- missing(Z)
  boundary_X_missing <- missing(boundary_X)
  boundary_formula_missing <- missing(boundary_formula)
  tau_init_missing <- missing(tau_init)
  sigma2_init_missing <- missing(sigma2_init)
  y_expr <- substitute(y)
  time_expr <- substitute(time)
  curve_expr <- substitute(curve)
  X_expr <- substitute(X)
  Z_expr <- substitute(Z)
  boundary_X_expr <- substitute(boundary_X)
  boundary_formula_expr <- substitute(boundary_formula)

  boundary_X_supplied <- !boundary_X_missing && !identical(boundary_X_expr, quote(NULL))
  boundary_formula_supplied <- !boundary_formula_missing &&
    !identical(boundary_formula_expr, quote(NULL))
  if (boundary_X_supplied && boundary_formula_supplied) {
    stop("Supply only one of `boundary_X` and `boundary_formula`.")
  }

  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("`data` must be a data frame.")
    }
    if (y_missing || time_missing || curve_missing) {
      stop("When `data` is supplied, `y`, `time`, and `curve` must be supplied.")
    }
    y <- .resolve_fit_data_vector(y_expr, data = data, env = caller, what = "y")
    time <- .resolve_fit_data_vector(time_expr, data = data, env = caller, what = "time")
    curve <- .resolve_fit_data_vector(curve_expr, data = data, env = caller, what = "curve")
    X <- if (X_missing) NULL else .resolve_fit_design_arg(X_expr, data = data, env = caller, what = "X")
    Z <- if (Z_missing) NULL else .resolve_fit_design_arg(Z_expr, data = data, env = caller, what = "Z")
    boundary_X <- if (boundary_X_missing) {
      NULL
    } else {
      .resolve_fit_design_arg(
        boundary_X_expr,
        data = data,
        env = caller,
        what = "boundary_X",
        allow_formula = FALSE
      )
    }
    if (boundary_formula_supplied) {
      boundary_X <- .resolve_fit_formula_arg(
        boundary_formula_expr,
        data = data,
        env = caller,
        what = "boundary_formula"
      )
    }
  } else if (boundary_formula_supplied) {
    boundary_X <- .resolve_fit_formula_arg(
      boundary_formula_expr,
      data = NULL,
      env = caller,
      what = "boundary_formula"
    )
  }
  if (inherits(boundary_X, "formula")) {
    stop("`boundary_X` must be a design matrix, not a formula. Use `boundary_formula` for formulas.")
  }

  trace_quad_n <- as.integer(trace_quad_n)
  if (length(trace_quad_n) != 1L || is.na(trace_quad_n) || trace_quad_n < 1L) {
    stop("`trace_quad_n` must be a positive integer.")
  }

  operator <- .coerce_operator_k1(operator)
  left_boundary <- .coerce_left_boundary_k1(left_boundary)
  right_boundary <- .coerce_right_boundary_operator_k1(right_boundary)
  .validate_left_boundary_k1(left_boundary)
  .validate_right_boundary_operator_k1(right_boundary)

  prepared <- .prepare_boundary_value_data(
    y = y,
    time = time,
    curve = curve,
    X = X,
    Z = Z,
    boundary_X = boundary_X,
    intercept_default = FALSE
  )
  prepared <- .augment_boundary_value_design_k1(
    prepared = prepared,
    right_boundary = right_boundary,
    left_boundary = left_boundary
  )

  y_var <- stats::var(prepared$y)
  if (!is.finite(y_var) || y_var <= 0) {
    y_var <- 1
  }

  tau_init <- if (is.null(tau_init)) 0.5 * y_var else as.numeric(tau_init)
  sigma2_init <- if (is.null(sigma2_init)) 0.1 * y_var else as.numeric(sigma2_init)
  sigma_u2_init <- if (is.null(sigma_u2_init)) 0.1 * y_var else as.numeric(sigma_u2_init)

  if (tau_init <= 0 || sigma2_init <= 0 || sigma_u2_init < 0) {
    stop("Initial variance values must be positive, with `sigma_u2_init >= 0`.")
  }

  lambda_d2_init <- if (tau_init_missing && sigma2_init_missing) {
    max(operator$lambda_start[["d2"]], 1e-8)
  } else {
    max(sigma2_init / tau_init, 1e-8)
  }
  if (isTRUE(operator$estimate_level)) {
    lambda_level_init <- if (is.null(lambda_level_init)) {
      max(operator$lambda_start[["level"]], lambda_d2_init, 1e-8)
    } else {
      as.numeric(lambda_level_init)
    }
    if (length(lambda_level_init) != 1L || !is.finite(lambda_level_init) || lambda_level_init <= 0) {
      stop("`lambda_level_init` must be positive for shifted operators.")
    }
  } else {
    if (!is.null(lambda_level_init)) {
      lambda_level_init <- as.numeric(lambda_level_init)
      if (length(lambda_level_init) != 1L || !is.finite(lambda_level_init) || lambda_level_init != 0) {
        stop("`lambda_level_init` must be 0 or NULL when the operator has no level term.")
      }
    }
    lambda_level_init <- 0
  }

  .fit_boundary_value_markussen_reml_k1(
    prepared          = prepared,
    operator          = operator,
    left_boundary     = left_boundary,
    right_boundary    = right_boundary,
    lambda_d2_init    = lambda_d2_init,
    sigma2_init       = sigma2_init,
    sigma_u2_init     = sigma_u2_init,
    lambda_level_init = lambda_level_init,
    control           = control,
    trace_quad_n      = trace_quad_n,
    call              = call
  )
}
