#' Markussen k = 1 operator
#'
#' Creates the `K.order = 1` operator used by \pkg{fdaMixed}. The supported
#' forms match the original package:
#'
#' \itemize{
#'   \item `lambda` of length one gives `L f = -lambda[1] * f''`.
#'   \item `lambda` of length two gives
#'     `L f = -lambda[1] * f'' + lambda[2] * f`.
#' }
#'
#' Positive coefficients are used as starting values. A zero level coefficient
#' fixes the level term at zero; a positive level coefficient lets RTMB estimate
#' it from the data.
#'
#' @param lambda Numeric vector of length one or two using the same convention
#'   as `fdaMixed::fdaLm(K.order = 1, lambda = ...)`.
#' @return An object of class `"operator_k1"`.
#' @export
operator_k1 <- function(lambda = 1) {
  lambda <- .normalize_operator_lambda_k1(lambda)
  estimate_level <- lambda[["level"]] > 0
  operator_label <- if (estimate_level) {
    "L f = -lambda_d2 * f'' + lambda_level * f"
  } else {
    "L f = -lambda_d2 * f''"
  }

  structure(
    list(
      K.order = 1L,
      lambda_start = lambda,
      estimate_level = estimate_level,
      operator_label = operator_label
    ),
    class = c("operator_k1", "operator_laplace")
  )
}

#' @rdname operator_k1
#' @export
operator_laplace <- function(lambda = 1) {
  operator_k1(lambda = lambda)
}

.normalize_operator_lambda_k1 <- function(lambda) {
  if (is.null(lambda)) {
    lambda <- 1
  }
  lambda <- as.numeric(lambda)
  if (length(lambda) < 1L || length(lambda) > 2L || any(!is.finite(lambda))) {
    stop("`lambda` must be a finite numeric vector of length one or two.")
  }
  if (any(lambda < 0) || lambda[[1L]] <= 0) {
    stop("`lambda[1]` must be positive and all lambda values must be non-negative.")
  }
  if (length(lambda) == 1L) {
    lambda <- c(lambda, 0)
  }
  names(lambda) <- c("d2", "level")
  lambda
}

.coerce_operator_k1 <- function(operator) {
  if (.is_laplace_operator(operator)) {
    return(operator)
  }
  if (is.numeric(operator)) {
    return(operator_k1(lambda = operator))
  }
  if (is.character(operator) && length(operator) == 1L) {
    key <- tolower(operator)
    if (key %in% c("laplace", "pure", "pure_laplace", "d2")) {
      return(operator_k1(lambda = 1))
    }
    if (key %in% c("shifted", "shifted_laplace", "fdaMixed_shifted")) {
      return(operator_k1(lambda = c(1, 1)))
    }
  }

  stop(
    "`operator` must be created by `operator_k1()`/`operator_laplace()`, ",
    "be a fdaMixed-style lambda vector, or be one of 'laplace' or ",
    "'shifted_laplace'."
  )
}

.is_laplace_operator <- function(x) {
  inherits(x, "operator_k1") || inherits(x, "operator_laplace")
}
