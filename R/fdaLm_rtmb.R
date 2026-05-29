#' Formula interface to the RTMB fit
#'
#' `fdaLm_rtmb()` is the main user-facing interface for \pkg{fdaMixedRTMB}.
#' It mimics the `fdaLm()` interface from the original \pkg{fdaMixed} package
#' and uses the \pkg{Formula} package to parse a two-part left-hand-side
#' formula of the form
#' \preformatted{ y | curve ~ fixed_part }
#'
#' Variables on either right-hand side are first looked up in `design`
#' (a curve-level data frame, one row per curve) when available, and then in
#' `data` (the observation-level data frame). This matches the
#' "marginal vs. global ANOVA" behavior of the original `fdaLm()`.
#'
#' When `boundary_value = TRUE`, the right boundary value is estimated through
#' a fixed null-space carrier. The optional `boundary_formula` is a one-sided
#' formula for the curve-level boundary regression. Its variables are looked up
#' in `design` if available; otherwise they are collapsed from the
#' observation-level `data`, in which case they must be constant within each
#' curve. When `boundary_value = FALSE`, no right-boundary value is estimated
#' and the serial effect is fit in the homogeneous right-boundary space.
#'
#' @param formula A `Formula::Formula` of the form
#'   `y | curve ~ fixed_part | random_part`. The left-hand side may have one
#'   or two parts; the right-hand side may have one or two parts.
#' @param data Observation-level data frame containing the response, the
#'   curve identifier (if specified on the LHS), the time column named by
#'   `time_variable`, and any observation-level covariates referenced on the
#'   right-hand side.
#' @param design Optional curve-level data frame with one row per curve, in
#'   the same order as the levels of the curve identifier. Used to resolve
#'   covariates that vary between curves but are constant within them.
#' @param left_boundary,right_boundary Boundary specifications as a numeric
#'   length-two vector `c(alpha0, alpha1)` or a one-row two-column matrix.
#' @param operator An `operator_k1()` object, a fdaMixed-style `lambda` vector,
#'   or `"laplace"`/`"shifted_laplace"`.
#' @param boundary_formula Optional one-sided formula for the curve-level
#'   boundary regression. Only used when `boundary_value = TRUE`.
#' @param boundary_value Logical; if `TRUE` (default), estimate the
#'   curve-specific right-boundary value. If `FALSE`, omit the boundary-value
#'   fixed component and use a homogeneous right boundary for the serial
#'   effect.
#' @param time_variable Name of the column in `data` holding the observation
#'   times. Defaults to `"time"`.
#' @param tau_init Initial serial-effect scale.
#' @param sigma2_init Initial noise variance.
#' @param lambda_level_init Initial zeroth-order penalty for shifted
#'   operators.
#' @param trace_quad_n Number of fixed Gauss-Legendre quadrature nodes for the
#'   Markussen trace-integral REML log-determinant correction (default 64).
#' @param control Optional list passed to `stats::nlminb()`.
#' @return An object of class `"boundary_value_fit"` with the original call
#'   stored in `$call` and the parsed formula stored in `$formula`.
#' @export
fdaLm_rtmb <- function(
  formula,
  data,
  design = NULL,
  left_boundary = c(1, 0),
  right_boundary = c(1, 0),
  operator = operator_k1(),
  boundary_formula = NULL,
  boundary_value = TRUE,
  time_variable = "time",
  tau_init = NULL,
  sigma2_init = NULL,
  lambda_level_init = NULL,
  trace_quad_n = 64L,
  control = list()
) {
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop(
      "`fdaLm_rtmb()` requires the 'Formula' package. ",
      "Install it with install.packages('Formula')."
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("`data` must be supplied as a data frame.")
  }
  if (!is.null(design) && !is.data.frame(design)) {
    stop("`design` must be `NULL` or a data frame.")
  }
  if (!is.character(time_variable) || length(time_variable) != 1L) {
    stop("`time_variable` must be a single column name.")
  }
  if (!is.logical(boundary_value) || length(boundary_value) != 1L || is.na(boundary_value)) {
    stop("`boundary_value` must be a single `TRUE` or `FALSE` value.")
  }
  if (!boundary_value && !is.null(boundary_formula)) {
    stop("`boundary_formula` can only be supplied when `boundary_value = TRUE`.")
  }
  if (!time_variable %in% names(data)) {
    stop(
      "Column `", time_variable, "` not found in `data`. ",
      "Set `time_variable` to the name of the time column."
    )
  }

  call <- match.call()

  fml <- Formula::Formula(formula)
  parts <- length(fml)
  lhs_n <- parts[1L]
  rhs_n <- parts[2L]

  if (lhs_n < 1L) {
    stop("`formula` must have a response on the left-hand side.")
  }
  if (rhs_n < 1L) {
    stop("`formula` must have at least one right-hand side.")
  }

  formula_env <- environment(formula)
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }

  ##response and curve identifier 
  y_form <- stats::formula(fml, lhs = 1L, rhs = 0L)
  y_expr <- y_form[[2L]]
  y <- eval(y_expr, envir = data, enclos = formula_env)

  if (lhs_n >= 2L) {
    curve_form <- stats::formula(fml, lhs = 2L, rhs = 0L)
    curve_expr <- curve_form[[2L]]
    curve <- eval(curve_expr, envir = data, enclos = formula_env)
  } else {
    curve <- factor(rep(1L, length(y)))
  }
  curve <- factor(curve)
  curve_levels <- levels(curve)
  n_curves <- length(curve_levels)
  curve_idx_per_obs <- as.integer(curve)

  ##Time column 
  time <- data[[time_variable]]

  ## helper to resolve a one-sided RHS 
  ## With estimated boundary values, the curve-level boundary regression already
  ## supplies a constant carrier, so an explicit fixed intercept would be
  ## collinear. Homogeneous-boundary fits keep the ordinary intercept.
  resolve_rhs <- function(rhs_index, drop_intercept = FALSE) {
    rhs_form <- stats::formula(fml, lhs = 0L, rhs = rhs_index)
    rhs_terms <- stats::terms(rhs_form)
    if (isTRUE(drop_intercept)) {
      attr(rhs_terms, "intercept") <- 0L
    }
    rhs_vars <- all.vars(rhs_terms)
    in_design <- !is.null(design) &&
      length(rhs_vars) > 0L &&
      all(rhs_vars %in% names(design))
    if (in_design) {
      mf <- stats::model.frame(rhs_terms, data = design)
      mat <- stats::model.matrix(rhs_terms, mf)
      if (nrow(mat) != n_curves) {
        stop(
          "`design` must have one row per curve (expected ",
          n_curves, ", got ", nrow(mat), ")."
        )
      }
      mat <- mat[curve_idx_per_obs, , drop = FALSE]
    } else {
      mf <- stats::model.frame(rhs_terms, data = data)
      mat <- stats::model.matrix(rhs_terms, mf)
      if (nrow(mat) != length(y)) {
        stop(
          "Right-hand side could not be expanded to the observation grid: ",
          "got ", nrow(mat), " rows, expected ", length(y), "."
        )
      }
    }
    if (ncol(mat) == 0L) {
      return(NULL)
    }
    mat
  }

  X_mat <- if (rhs_n >= 1L) {
    resolve_rhs(1L, drop_intercept = boundary_value)
  } else {
    NULL
  }

  ##  emmeans metadata 
  ## Store the fixed-effects terms, factor levels, contrasts, and source data
  ## so that recover_data() / emm_basis() can reconstruct the reference grid.
  fix_terms_store    <- NULL
  fix_xlev_store     <- NULL
  fix_contrasts_store <- NULL
  model_data_store   <- NULL

  if (rhs_n >= 1L) {
    rhs_form_fix <- stats::formula(fml, lhs = 0L, rhs = 1L)
    fix_t        <- stats::terms(rhs_form_fix)
    if (boundary_value) {
      attr(fix_t, "intercept") <- 0L        # match drop_intercept = TRUE
    }
    fix_vars <- all.vars(fix_t)

    if (length(fix_vars) > 0L) {
      in_design_fix <- !is.null(design) &&
        all(fix_vars %in% names(design))
      fix_src <- if (in_design_fix) {
        design[curve_idx_per_obs, , drop = FALSE]
      } else {
        data
      }

      fix_mf_store <- tryCatch(
        stats::model.frame(fix_t, data = fix_src,
                           na.action = stats::na.pass),
        error = function(e) NULL
      )
      if (!is.null(fix_mf_store)) {
        fix_terms_store <- attr(fix_mf_store, "terms")
        fix_xlev_store  <- lapply(
          fix_mf_store,
          function(v) if (is.factor(v)) levels(v) else NULL
        )
        fix_xlev_store <- fix_xlev_store[!vapply(fix_xlev_store, is.null,
                                                  logical(1L))]
        kept <- fix_vars[fix_vars %in% names(fix_src)]
        model_data_store <- fix_src[, kept, drop = FALSE]
      }
    }
    if (!is.null(X_mat)) {
      fix_contrasts_store <- attr(X_mat, "contrasts")
    }
  }

  ## boundary formula
  boundary_X <- NULL
  boundary_terms_store     <- NULL
  boundary_xlev_store      <- NULL
  boundary_contrasts_store <- NULL
  boundary_model_data_store <- NULL

  if (boundary_value && !is.null(boundary_formula)) {
    bf <- stats::as.formula(boundary_formula)
    if (length(bf) != 2L) {
      stop("`boundary_formula` must be a one-sided formula.")
    }
    bterms <- stats::terms(bf)
    bvars <- all.vars(bterms)

    in_design <- !is.null(design) &&
      length(bvars) > 0L &&
      all(bvars %in% names(design))
    if (in_design) {
      bmf <- stats::model.frame(bterms, data = design)
      boundary_X <- stats::model.matrix(bterms, bmf)
      if (nrow(boundary_X) != n_curves) {
        stop(
          "`design` must have one row per curve for `boundary_formula` ",
          "(expected ", n_curves, ", got ", nrow(boundary_X), ")."
        )
      }
    } else {
      bmf <- stats::model.frame(bterms, data = data)
      bmat <- stats::model.matrix(bterms, bmf)
      if (nrow(bmat) != length(y)) {
        stop("`boundary_formula` could not be expanded against `data`.")
      }
      collapsed <- matrix(NA_real_, nrow = n_curves, ncol = ncol(bmat))
      colnames(collapsed) <- colnames(bmat)
      for (i in seq_len(n_curves)) {
        rows <- which(curve_idx_per_obs == i)
        ref <- bmat[rows[1L], , drop = FALSE]
        same <- apply(
          abs(bmat[rows, , drop = FALSE] - ref[rep(1L, length(rows)), , drop = FALSE]) < 1e-10,
          2L,
          all
        )
        if (!all(same)) {
          stop(
            "`boundary_formula` variables must be constant within each curve ",
            "when `design` is not supplied."
          )
        }
        collapsed[i, ] <- ref
      }
      boundary_X <- collapsed
    }

    ##  emmeans metadata for boundary regression 
    ## Use curve-level source data: design if available, otherwise take the
    ## first observation of each curve (boundary vars are constant within curves).
    boundary_src <- if (in_design) {
      design
    } else {
      data[!duplicated(curve_idx_per_obs), , drop = FALSE]
    }
    b_mf_store <- tryCatch(
      stats::model.frame(bterms, data = boundary_src,
                         na.action = stats::na.pass),
      error = function(e) NULL
    )
    if (!is.null(b_mf_store)) {
      boundary_terms_store <- attr(b_mf_store, "terms")
      bxlev_tmp            <- lapply(
        b_mf_store,
        function(v) if (is.factor(v)) levels(v) else NULL
      )
      boundary_xlev_store <- bxlev_tmp[
        !vapply(bxlev_tmp, is.null, logical(1L))
      ]
      kept_b <- bvars[bvars %in% names(boundary_src)]
      boundary_model_data_store <- boundary_src[, kept_b, drop = FALSE]
      bmat_tmp <- tryCatch(
        stats::model.matrix(bterms, b_mf_store),
        error = function(e) NULL
      )
      boundary_contrasts_store <- if (!is.null(bmat_tmp)) attr(bmat_tmp, "contrasts") else NULL
    }
  } else if (!boundary_value) {
    boundary_X <- matrix(numeric(0L), nrow = n_curves, ncol = 0L)
    colnames(boundary_X) <- character(0L)
  }

  result <- .fit_boundary_value_rtmb(
    y = y,
    time = time,
    curve = curve,
    X = X_mat,
    boundary_X = boundary_X,
    operator = operator,
    left_boundary = left_boundary,
    right_boundary = right_boundary,
    tau_init = tau_init,
    sigma2_init = sigma2_init,
    lambda_level_init = lambda_level_init,
    trace_quad_n = trace_quad_n,
    control = control
  )

  result$call             <- call
  result$formula          <- formula
  result$boundary_value   <- boundary_value
  result$boundary_formula <- boundary_formula
  result$fix_terms              <- fix_terms_store
  result$fix_xlev               <- fix_xlev_store
  result$fix_contrasts          <- fix_contrasts_store
  result$model_data             <- model_data_store
  result$boundary_terms         <- boundary_terms_store
  result$boundary_xlev          <- boundary_xlev_store
  result$boundary_contrasts     <- boundary_contrasts_store
  result$boundary_model_data    <- boundary_model_data_store
  result
}
