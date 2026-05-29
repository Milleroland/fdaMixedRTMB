#' Methods for boundary-value fits
#'
#' S3 methods for fitted boundary-value models, including standard generics
#' from \pkg{stats} and a home-grown \code{fixef} generic for
#' compatibility with \pkg{emmeans} and the broader model ecosystem.
#'
#' @param object A fitted `"boundary_value_fit"` object.
#' @param x A fitted `"boundary_value_fit"` object, or its summary.
#' @param component Which coefficient or variance block to return.
#' @param type Prediction or residual type.
#' @param digits Printing precision.
#' @param ... Unused extras for generic compatibility.
#' @return Method-specific summaries, coefficients, predictions, or the input
#'   object invisibly for print methods.
#' @name boundary_value_fit_methods
NULL


# This generic is also defined in nlme and lme4.  We provide our own so
# that users who do not have those packages installed can still call
# fixef(fit).  If nlme or lme4 is attached *before* fdaMixedRTMB, R will
# mask our generic with theirs; fixef.boundary_value_fit will still
# dispatch correctly in all cases.

#' Extract fixed-effect coefficients
#'
#' Returns the ordinary (non-boundary) fixed-effect coefficient vector
#' for a fitted model object.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A named numeric vector of fixed-effect coefficients.
#' @export
fixef <- function(object, ...) UseMethod("fixef")

#' @rdname boundary_value_fit_methods
#' @export
fixef.boundary_value_fit <- function(object, ...) {
  object$beta
}

.compact_deparse <- function(x) {
  paste(deparse(x, width.cutoff = 500L), collapse = " ")
}

.boundary_model_label <- function(object) {
  has_boundary <- length(object$boundary_coef) > 0L
  if (!has_boundary) {
    return("none (homogeneous right boundary)")
  }

  if (!is.null(object$boundary_formula)) {
    boundary_formula <- stats::as.formula(object$boundary_formula)
    return(paste("rho ~", .compact_deparse(boundary_formula[[2L]])))
  }

  if (!is.null(object$boundary_terms)) {
    labels <- attr(object$boundary_terms, "term.labels")
    has_intercept <- identical(attr(object$boundary_terms, "intercept"), 1L)
    rhs <- c(if (has_intercept) "1", labels)
    if (!length(rhs)) {
      rhs <- "0"
    }
    return(paste("rho ~", paste(rhs, collapse = " + ")))
  }

  boundary_X <- object$data$boundary_X
  if (!is.null(boundary_X) && ncol(boundary_X) > 0L) {
    cols <- colnames(boundary_X)
    if (is.null(cols)) {
      cols <- paste0("b", seq_len(ncol(boundary_X)))
    }
    cols[cols == "(Intercept)"] <- "1"
    return(paste("rho ~", paste(cols, collapse = " + ")))
  }

  "rho ~ 1"
}

.boundary_grid_labels <- function(grid) {
  if (is.null(grid) || ncol(grid) == 0L) {
    return("(all)")
  }

  apply(
    grid,
    1L,
    function(row) {
      paste(paste(names(row), as.character(row), sep = "="), collapse = ", ")
    }
  )
}

.boundary_model_rho <- function(object) {
  if (length(object$boundary_coef) == 0L) {
    return(numeric(0L))
  }

  if (!is.null(object$boundary_terms) && !is.null(object$boundary_model_data)) {
    dat <- object$boundary_model_data
    if (ncol(dat) == 0L) {
      grid <- dat[1L, , drop = FALSE]
      X_grid <- matrix(1, nrow = 1L, ncol = length(object$boundary_coef))
      colnames(X_grid) <- names(object$boundary_coef)
      labels <- "(all)"
    } else {
      grid <- dat[!duplicated(dat), , drop = FALSE]
      mf <- stats::model.frame(
        object$boundary_terms,
        data = grid,
        na.action = stats::na.pass,
        xlev = object$boundary_xlev
      )
      X_grid <- stats::model.matrix(
        object$boundary_terms,
        mf,
        contrasts.arg = object$boundary_contrasts
      )
      labels <- .boundary_grid_labels(grid)
    }
    out <- drop(X_grid %*% unname(object$boundary_coef))
    return(stats::setNames(out, labels))
  }

  boundary_X <- object$data$boundary_X
  if (is.null(boundary_X) || ncol(boundary_X) == 0L) {
    return(numeric(0L))
  }

  keep <- !duplicated(as.data.frame(boundary_X))
  grid <- boundary_X[keep, , drop = FALSE]
  out <- drop(grid %*% unname(object$boundary_coef))
  labels <- .boundary_grid_labels(as.data.frame(grid))
  stats::setNames(out, labels)
}

.print_named_numeric <- function(x, digits, name_col = "name") {
  if (!length(x)) {
    return(invisible(x))
  }

  nm <- names(x)
  if (is.null(nm)) {
    nm <- as.character(seq_along(x))
  }
  out <- stats::setNames(
    data.frame(nm, round(unname(x), digits), check.names = FALSE),
    c(name_col, "estimate")
  )
  print(out, row.names = FALSE, max = nrow(out) * ncol(out))
  invisible(x)
}

#' @rdname boundary_value_fit_methods
#' @export
coef.boundary_value_fit <- function(
  object,
  component = c("beta", "boundary", "variance", "rho", "all"),
  ...
) {
  component <- match.arg(component)

  switch(
    component,
    beta     = object$beta,
    boundary = object$boundary_coef,
    variance = c(tau = object$tau, sigma2 = object$sigma2),
    rho      = object$curve_rho,
    list(
      beta     = object$beta,
      boundary = object$boundary_coef,
      variance = c(tau = object$tau, sigma2 = object$sigma2),
      rho      = object$curve_rho
    )
  )
}

#' @rdname boundary_value_fit_methods
#' @export
predict.boundary_value_fit <- function(
  object,
  type = c("response", "fixed", "boundary", "serial", "right_derivative", "rho"),
  ...
) {
  type <- match.arg(type)

  switch(
    type,
    response = object$fitted,
    fixed = object$fixed_effect,
    boundary = object$boundary_effect,
    serial = object$serial_effect,
    right_derivative = object$right_boundary_derivative,
    rho = object$curve_rho
  )
}

#' @rdname boundary_value_fit_methods
#' @export
summary.boundary_value_fit <- function(object, ...) {
  kappa_hat <- if (!is.null(object$lambda_level) && !is.null(object$lambda_d2) &&
                   is.finite(object$lambda_level) && is.finite(object$lambda_d2) &&
                   object$lambda_d2 > 0) {
    sqrt(object$lambda_level / object$lambda_d2)
  } else {
    NA_real_
  }
  out <- list(
    call = object$call,
    method = object$method,
    engine = object$engine,
    converged = object$converged,
    logLik = object$logLik,
    logLik_type = object$logLik_type,
    operator = c(tau = object$tau, lambda_d2 = object$lambda_d2,
                 lambda_level = object$lambda_level, kappa = kappa_hat),
    operator_label = object$operator$operator_label,
    noise = c(sigma2 = object$sigma2),
    beta = object$beta,
    boundary_coef = object$boundary_coef,
    boundary_model = .boundary_model_label(object),
    boundary_model_rho = .boundary_model_rho(object),
    right_operator = .right_boundary_label_k1(object$right_boundary),
    curve_rho = object$curve_rho,
    right_boundary_derivative = object$right_boundary_derivative,
    right_boundary_serial_operator = object$right_boundary_serial_operator,
    n_curves = object$data$n_curves,
    n_per_curve = object$data$n_per_curve
  )

  class(out) <- "summary.boundary_value_fit"
  out
}

#' @rdname boundary_value_fit_methods
#' @export
print.boundary_value_fit <- function(x, digits = 5L, ...) {
  cat("boundary_value_fit\n")
  cat("  converged:", x$converged, "\n")
  cat("  curves:", x$data$n_curves, "  points/curve:", x$data$n_per_curve, "\n")
  kappa_val <- if (!is.null(x$lambda_level) && !is.null(x$lambda_d2) &&
                   is.finite(x$lambda_level) && is.finite(x$lambda_d2) &&
                   x$lambda_d2 > 0) {
    sqrt(x$lambda_level / x$lambda_d2)
  } else {
    NA_real_
  }
  cat("  right operator:", .right_boundary_label_k1(x$right_boundary), "\n")
  cat("  boundary model: ", .boundary_model_label(x), "\n", sep = "")
  cat("  operator:", x$operator$operator_label, "\n")
  cat("  tau:", formatC(x$tau, digits = digits, format = "f"),
      " lambda_d2:", formatC(x$lambda_d2, digits = digits, format = "f"),
      " lambda_level:", formatC(x$lambda_level, digits = digits, format = "f"),
      " kappa:", formatC(kappa_val, digits = digits, format = "f"), "\n")
  cat("  sigma2:", formatC(x$sigma2, digits = digits, format = "f"), "\n")
  boundary_rho <- .boundary_model_rho(x)
  if (length(boundary_rho)) {
    cat("  rho by boundary model:\n")
    .print_named_numeric(boundary_rho, digits, name_col = "boundary")
  } else {
    cat("  curve-specific rho:\n")
    .print_named_numeric(x$curve_rho, digits)
  }
  invisible(x)
}

#' @rdname boundary_value_fit_methods
#' @export
print.summary.boundary_value_fit <- function(x, digits = max(3L, getOption("digits") - 2L), ...) {
  cat("Boundary-value fit summary\n")
  cat("  converged:", x$converged, "\n")
  logLik_label <- if (identical(x$logLik_type, "restricted")) "restricted logLik:" else "logLik:"
  cat("  ", logLik_label, " ", formatC(x$logLik, digits = digits, format = "f"), "\n", sep = "")
  cat("  right operator:", x$right_operator, "\n")
  cat("  boundary model: ", x$boundary_model, "\n", sep = "")
  cat("  operator:", x$operator_label, "\n")
  cat("    tau =", formatC(x$operator["tau"], digits = digits, format = "f"),
      " lambda_d2 =", formatC(x$operator["lambda_d2"], digits = digits, format = "f"),
      " lambda_level =", formatC(x$operator["lambda_level"], digits = digits, format = "f"),
      " kappa =", formatC(x$operator["kappa"], digits = digits, format = "f"), "\n")
  cat("  noise:", "sigma2 =", formatC(x$noise["sigma2"], digits = digits, format = "f"), "\n")

  if (length(x$beta)) {
    cat("\nOrdinary fixed effects\n")
    print(round(x$beta, digits), quote = FALSE)
  }
  if (length(x$boundary_coef)) {
    cat("\nBoundary-value coefficients\n")
    print(round(x$boundary_coef, digits), quote = FALSE)
  }
  if (length(x$boundary_model_rho)) {
    cat("\nBoundary-model rho\n")
    .print_named_numeric(x$boundary_model_rho, digits, name_col = "boundary")
  } else {
    cat("\nCurve-specific rho\n")
    .print_named_numeric(x$curve_rho, digits)
  }

  cat("\nHomogeneous serial boundary check\n")
  .print_named_numeric(x$right_boundary_serial_operator, digits)
  invisible(x)
}

# vcov

#' Variance-covariance matrix for fixed-effect coefficients
#'
#' Returns the approximate variance-covariance matrix for the fixed-effect
#' coefficients of a \code{boundary_value_fit} model from the RTMB marginal
#' objective Hessian.
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param component Which block to return:
#'   \describe{
#'     \item{\code{"beta"}}{Covariance for ordinary fixed effects (the
#'       \code{X} design matrix columns, excluding boundary-value columns).
#'       This is the default and the block used by \pkg{emmeans}.}
#'     \item{\code{"boundary"}}{Covariance for the boundary-value regression
#'       coefficients.}
#'     \item{\code{"all"}}{Full joint covariance for both blocks.}
#'   }
#' @param ... Unused.
#' @return A named symmetric numeric matrix.
#' @export
vcov.boundary_value_fit <- function(
  object,
  component = c("beta", "boundary", "all"),
  ...
) {
  component <- match.arg(component)

  full_vcov <- if (!is.null(object$fixed_vcov)) {
    object$fixed_vcov
  } else {
    object$sigma2 * object$Cbeta
  }
  cn_full <- colnames(object$data$X)
  if (!is.null(cn_full) && length(cn_full) == nrow(full_vcov)) {
    rownames(full_vcov) <- colnames(full_vcov) <- cn_full
  }

  base_cols     <- object$data$fixed_column_groups$beta
  boundary_cols <- object$data$fixed_column_groups$boundary

  switch(
    component,
    beta = {
      if (length(base_cols) == 0L) {
        return(matrix(numeric(0L), 0L, 0L))
      }
      v  <- full_vcov[base_cols, base_cols, drop = FALSE]
      cn <- colnames(object$data$base_X)
      if (!is.null(cn)) rownames(v) <- colnames(v) <- cn
      v
    },
    boundary = {
      if (length(boundary_cols) == 0L) {
        return(matrix(numeric(0L), 0L, 0L))
      }
      v  <- full_vcov[boundary_cols, boundary_cols, drop = FALSE]
      cn <- colnames(object$data$boundary_value_design)
      if (!is.null(cn)) rownames(v) <- colnames(v) <- cn
      v
    },
    all = full_vcov
  )
}

# model.matrix

#' Model matrix for fixed effects
#'
#' Returns the design matrix used for the ordinary fixed effects
#' (the \code{X} argument, excluding boundary-value columns).
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param component \code{"beta"} (default) for the ordinary fixed-effect
#'   columns only; \code{"all"} for the full augmented matrix including
#'   boundary-value columns.
#' @param ... Unused.
#' @return A numeric matrix.
#' @export
model.matrix.boundary_value_fit <- function(
  object,
  component = c("beta", "all"),
  ...
) {
  component <- match.arg(component)
  if (identical(component, "beta")) object$data$base_X else object$data$X
}

# residuals

#' Residuals from a boundary-value fit
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param type \code{"response"} (default) for raw residuals
#'   \eqn{y - \hat{y}}, or \code{"pearson"} for residuals standardised by
#'   \eqn{\sqrt{\hat{\sigma}^2}}.
#' @param ... Unused.
#' @return A numeric vector of residuals.
#' @export
residuals.boundary_value_fit <- function(
  object,
  type = c("response", "pearson"),
  ...
) {
  type  <- match.arg(type)
  resid <- object$data$y - object$fitted
  if (identical(type, "pearson")) {
    resid <- resid / sqrt(object$sigma2)
  }
  resid
}

# nobs

#' Number of observations
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param ... Unused.
#' @return A single integer.
#' @export
nobs.boundary_value_fit <- function(object, ...) {
  length(object$data$y)
}

# df.residual

#' Residual degrees of freedom
#'
#' Returns \eqn{n - p} where \eqn{n} is the total number of observations and
#' \eqn{p} is the total number of fixed-effect columns (ordinary \emph{and}
#' boundary-value).
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param ... Unused.
#' @return A single integer.
#' @export
df.residual.boundary_value_fit <- function(object, ...) {
  length(object$data$y) - ncol(object$data$X)
}

# logLik

#' Log-likelihood of a boundary-value fit
#'
#' Returns a \code{"logLik"} object with the fitted marginal log-likelihood
#' value stored in \code{object$logLik}.
#'
#' @param object A fitted \code{boundary_value_fit} object.
#' @param ... Unused.
#' @return An object of class \code{"logLik"}.
#' @export
logLik.boundary_value_fit <- function(object, ...) {
  val <- object$logLik
  attr(val, "nobs") <- length(object$data$y)
  attr(val, "df")   <- ncol(object$data$X)
  class(val)        <- "logLik"
  val
}

# terms

#' Terms object for the fixed-effect formula
#'
#' Returns the \code{terms} object for the ordinary fixed-effect part of the
#' model.  This is only available when the model was fit via
#' \code{\link{fdaLm_rtmb}()} (or when \code{fix_terms} was stored manually).
#' Returns \code{NULL} otherwise.
#'
#' @param x A fitted \code{boundary_value_fit} object.
#' @param ... Unused.
#' @return A \code{terms} object, or \code{NULL}.
#' @export
terms.boundary_value_fit <- function(x, ...) {
  x$fix_terms
}

# emmeans support

#' emmeans support for boundary_value_fit objects
#'
#' These methods make \code{boundary_value_fit} models compatible with the
#' \pkg{emmeans} package.  Two parts of the model can be targeted via the
#' \code{mode} argument:
#'
#' \describe{
#'   \item{\code{"fixed"} (default)}{Estimated marginal means for the
#'     \emph{ordinary} fixed effects (\eqn{\beta} from the \code{X} /
#'     \code{fixed_part} design matrix).}
#'   \item{\code{"boundary"}}{Estimated marginal means for the
#'     \emph{boundary-value regression} (\eqn{\rho} as a function of
#'     curve-level predictors from \code{boundary_formula}).}
#' }
#'
#' Both modes require the model to have been fit via \code{\link{fdaLm_rtmb}()}
#' so that the terms objects and reference data are stored automatically.
#'
#' @examples
#' \dontrun{
#' library(emmeans)
#' fit <- fdaLm_rtmb(y | curve ~ group,
#'                   data = dat, design = des,
#'                   boundary_formula = ~ group)
#'
#' # EMMs for ordinary fixed effects
#' emmeans(fit, ~ group)
#' pairs(emmeans(fit, ~ group))
#'
#' # EMMs for the boundary-value regression (rho ~ group)
#' emmeans(fit, ~ group, mode = "boundary")
#' pairs(emmeans(fit, ~ group, mode = "boundary"))
#' }
#'
#' @name emmeans_support
NULL

#' @rdname emmeans_support
#' @param object A fitted \code{boundary_value_fit} object.
#' @param mode \code{"fixed"} (default) or \code{"boundary"}.
#' @param data Ignored; absorbed for compatibility with the \pkg{emmeans}
#'   internal call signature.
#' @param params Ignored; absorbed for compatibility with the \pkg{emmeans}
#'   internal call signature.
#' @param ... Additional arguments (ignored).
#' @importFrom emmeans recover_data
#' @export
recover_data.boundary_value_fit <- function(object,
                                            mode = c("fixed", "boundary"),
                                            data = NULL,
                                            params = NULL,
                                            ...) {
  mode <- match.arg(mode)

  if (identical(mode, "boundary")) {
    trms <- object$boundary_terms
    dat  <- object$boundary_model_data
    if (is.null(trms)) {
      stop(
        "emmeans boundary mode requires 'boundary_terms' in the fitted object.\n",
        "Fit the model via fdaLm_rtmb() with a 'boundary_formula' for\n",
        "automatic support, or assign 'boundary_terms' and\n",
        "'boundary_model_data' manually."
      )
    }
  } else {
    trms <- object$fix_terms
    dat  <- object$model_data
    if (is.null(trms)) {
      stop(
        "emmeans support requires 'fix_terms' in the fitted object.\n",
        "Fit the model via fdaLm_rtmb() for automatic emmeans support, or\n",
        "assign 'fix_terms' and 'model_data' manually."
      )
    }
  }

  # Build the recovery data frame directly from the stored terms and data,
  # without delegating to recover_data.call().  The emmeans generic passes
  # data = NULL (and sometimes params) as named arguments; we absorb them
  # above so they do not cause "formal argument matched by multiple actual
  # arguments" errors.
  result <- stats::model.frame(trms, data = dat, na.action = stats::na.pass)
  attr(result, "terms")      <- trms
  attr(result, "predictors") <- attr(trms, "term.labels")
  # emmeans checks for na.action attribute; supply an empty one (no rows excluded)
  attr(result, "na.action")  <- structure(integer(0), class = "omit")
  result
}

#' @rdname emmeans_support
#' @param trms A \code{terms} object (passed by \pkg{emmeans}).
#' @param xlev A list of factor levels (passed by \pkg{emmeans}).
#' @param grid A data frame of predictor values (passed by \pkg{emmeans}).
#' @importFrom emmeans emm_basis
#' @export
emm_basis.boundary_value_fit <- function(object, trms, xlev, grid,
                                         mode = c("fixed", "boundary"),
                                         ...) {
  mode <- match.arg(mode)

  if (identical(mode, "boundary")) {
    # Boundary formula keeps its intercept — do not suppress it.
    m      <- stats::model.frame(trms, grid,
                                  na.action = stats::na.pass, xlev = xlev)
    X_grid <- stats::model.matrix(trms, m,
                                   contrasts.arg = object$boundary_contrasts)
    bhat   <- object$boundary_coef
    V      <- vcov(object, component = "boundary")
    cn <- colnames(X_grid)
    if (length(bhat) == length(cn)) {
      names(bhat)            <- cn
      rownames(V) <- colnames(V) <- cn
    }
  } else {
    # Fixed-effects RHS was built with drop_intercept = TRUE; enforce here.
    trms_use <- trms
    attr(trms_use, "intercept") <- 0L
    m      <- stats::model.frame(trms_use, grid,
                                  na.action = stats::na.pass, xlev = xlev)
    X_grid <- stats::model.matrix(trms_use, m,
                                   contrasts.arg = object$fix_contrasts)
    bhat   <- object$beta
    V      <- vcov(object, component = "beta")
  }

  # Both modes share sigma2 and residual df from the joint RTMB model.
  df_res <- df.residual(object)
  list(
    X      = X_grid,
    bhat   = bhat,
    nbasis = matrix(NA_real_),
    V      = V,
    dffun  = function(k, dfargs) dfargs$df,
    dfargs = list(df = df_res)
  )
}

# plot

#' @rdname boundary_value_fit_methods
#' @param which Diagnostic panels to draw. \code{"all"} (default) draws every
#'   panel; pass a character vector from \code{"fit"}, \code{"residual"},
#'   \code{"qq"}, \code{"serial"}, \code{"rho"}, \code{"boundary"}, or
#'   numeric indices \code{1:6}.
#' @param ask Logical; if \code{TRUE} pause between panels.
#' @export
plot.boundary_value_fit <- function(x, which = "all", ask = FALSE, ...) {
  panel_names <- c("fit", "residual", "qq", "serial", "rho", "boundary")

  if (is.null(which) || identical(which, "all")) {
    which <- panel_names
  } else if (is.numeric(which)) {
    if (any(which < 1L | which > length(panel_names))) {
      stop("`which` must be between 1 and ", length(panel_names), ".")
    }
    which <- panel_names[which]
  } else {
    which <- match.arg(which, c("all", panel_names), several.ok = TRUE)
    if ("all" %in% which) which <- panel_names
  }

  resid     <- x$data$y - x$fitted
  std_resid <- resid / sqrt(x$sigma2)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  if (isTRUE(ask)) {
    graphics::par(ask = TRUE)
  } else if (length(which) > 1L) {
    n_col <- if (length(which) <= 2L) length(which) else 2L
    graphics::par(mfrow = c(ceiling(length(which) / n_col), n_col))
  }

  for (panel in which) {
    switch(panel,
      fit = {
        graphics::plot(x$fitted, x$data$y,
                       xlab = "Fitted value", ylab = "Observed value",
                       main = "Observed vs fitted", ...)
        graphics::abline(0, 1, col = "gray45", lwd = 1.5)
      },
      residual = {
        graphics::plot(x$fitted, resid,
                       xlab = "Fitted value", ylab = "Residual",
                       main = "Residuals vs fitted", ...)
        graphics::abline(h = 0, col = "gray45", lwd = 1.5)
      },
      qq = {
        stats::qqnorm(std_resid, main = "Normal Q-Q",
                      ylab = "Standardized residual", ...)
        stats::qqline(std_resid, col = "gray45", lwd = 1.5)
      },
      serial = {
        serial_mat <- matrix(x$serial_effect,
                             nrow = x$data$n_per_curve,
                             ncol = x$data$n_curves)
        graphics::matplot(x$data$time_grid, serial_mat,
                          type = "l", lty = 1,
                          xlab = "Time", ylab = "Serial effect",
                          main = "Serial effects", ...)
        graphics::abline(h = 0, col = "gray80")
      },
      rho = {
        xx <- seq_len(x$data$n_curves)
        graphics::plot(xx, x$curve_rho,
                       xaxt = "n", xlab = "Curve",
                       ylab = "Estimated rho", main = "Boundary values", ...)
        graphics::axis(1, at = xx, labels = x$data$curve_levels)
        graphics::abline(h = 0, col = "gray80")
      },
      boundary = {
        xx <- seq_len(x$data$n_curves)
        graphics::plot(xx, x$right_boundary_serial_operator,
                       xaxt = "n", xlab = "Curve",
                       ylab = "Right operator on serial effect",
                       main = "Homogeneous boundary check", ...)
        graphics::axis(1, at = xx, labels = x$data$curve_levels)
        graphics::abline(h = 0, col = "gray45", lwd = 1.5)
      }
    )
  }

  invisible(x)
}
