#' fdaMixedRTMB: operator-based functional mixed models
#'
#' A package for fitting operator-based functional mixed models with RTMB.
#' The main interface is \code{\link{fdaLm_rtmb}()}, a formula interface that
#' mirrors \code{fdaLm()} from \pkg{fdaMixed}. The package supports the
#' `K.order = 1` operators from the original \pkg{fdaMixed} package:
#' `L f = -lambda_d2 * f''` and
#' `L f = -lambda_d2 * f'' + lambda_level * f`. The user supplies operator
#' starting values, boundary operators as numeric length-two alpha vectors, and
#' formulas for the fixed, random, and boundary-value regressions.
#'
#' @keywords internal
#' @importFrom stats coef predict
"_PACKAGE"
