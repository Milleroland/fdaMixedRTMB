#' Markussen k = 1 operator
#'
#' Creates the `K.order = 1` operator used by \pkg{fdaMixed}. The supported
#' forms match the original package:
#'
#' \itemize{
#'   \item `lambda_start` of length one gives `L f = -lambda_d2 * f''`.
#'   \item `lambda_start` of length two gives
#'     `L f = -lambda_d2 * f'' + lambda_level * f`.
#' }
#'
#' Positive coefficients are used as starting values. A zero level coefficient
#' fixes the level term at zero; a positive level coefficient lets RTMB estimate
#' it from the data.
#'
#' @param lambda_start Numeric starting-value vector of length one or two, using
#'   the same convention as `fdaMixed::fdaLm(K.order = 1, lambda = ...)`.
#' @param lambda Compatibility alias for `lambda_start`.
#' @return An object of class `"operator_k1"`.
#' @export
operator_k1 <- function(lambda_start = 1, lambda = NULL) {
  if (!missing(lambda)) {
    if (!missing(lambda_start)) {
      stop("Supply only one of `lambda_start` and `lambda`.")
    }
    lambda_start <- lambda
  }
  lambda_start <- .normalize_operator_lambda_start_k1(lambda_start)
  estimate_level <- lambda_start[["level"]] > 0
  operator_label <- if (estimate_level) {
    "L f = -lambda_d2 * f'' + lambda_level * f"
  } else {
    "L f = -lambda_d2 * f''"
  }

  structure(
    list(
      K.order = 1L,
      lambda_start = lambda_start,
      estimate_level = estimate_level,
      operator_label = operator_label
    ),
    class = c("operator_k1", "operator_laplace")
  )
}

#' @rdname operator_k1
#' @export
operator_laplace <- function(lambda_start = 1, lambda = NULL) {
  if (!missing(lambda)) {
    if (!missing(lambda_start)) {
      stop("Supply only one of `lambda_start` and `lambda`.")
    }
    return(operator_k1(lambda_start = lambda))
  }
  operator_k1(lambda_start = lambda_start)
}

.normalize_operator_lambda_start_k1 <- function(lambda_start) {
  if (is.null(lambda_start)) {
    lambda_start <- 1
  }
  lambda_start <- as.numeric(lambda_start)
  if (length(lambda_start) < 1L ||
      length(lambda_start) > 2L ||
      any(!is.finite(lambda_start))) {
    stop("`lambda_start` must be a finite numeric vector of length one or two.")
  }
  if (any(lambda_start < 0) || lambda_start[[1L]] <= 0) {
    stop(
      "`lambda_start[1]` must be positive and all lambda_start values ",
      "must be non-negative."
    )
  }
  if (length(lambda_start) == 1L) {
    lambda_start <- c(lambda_start, 0)
  }
  names(lambda_start) <- c("d2", "level")
  lambda_start
}

.coerce_operator_k1 <- function(operator) {
  if (.is_laplace_operator(operator)) {
    return(operator)
  }
  if (is.numeric(operator)) {
    return(operator_k1(lambda_start = operator))
  }
  if (is.character(operator) && length(operator) == 1L) {
    key <- tolower(operator)
    if (key %in% c("laplace", "pure", "pure_laplace", "d2")) {
      return(operator_k1(lambda_start = 1))
    }
    if (key %in% c("shifted", "shifted_laplace", "fdaMixed_shifted")) {
      return(operator_k1(lambda_start = c(1, 1)))
    }
  }

  stop(
    "`operator` must be created by `operator_k1()`/`operator_laplace()`, ",
    "be a lambda_start vector, or be one of 'laplace' or ",
    "'shifted_laplace'."
  )
}

.is_laplace_operator <- function(x) {
  inherits(x, "operator_k1") || inherits(x, "operator_laplace")
}
