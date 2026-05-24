.coerce_design_matrix <- function(x, n_rows, what, intercept_default = FALSE) {
  if (is.null(x)) {
    if (!intercept_default) {
      return(matrix(0, nrow = n_rows, ncol = 0L))
    }

    out <- matrix(1, nrow = n_rows, ncol = 1L)
    colnames(out) <- "(Intercept)"
    return(out)
  }

  if (is.vector(x) && !is.list(x)) {
    x <- matrix(x, ncol = 1L)
  } else {
    x <- as.matrix(x)
  }

  if (nrow(x) != n_rows) {
    stop("`", what, "` must have ", n_rows, " rows.")
  }

  storage.mode(x) <- "double"
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`", what, "` must contain only finite numeric values.")
  }

  if (is.null(colnames(x))) {
    colnames(x) <- if (ncol(x)) paste0(what, seq_len(ncol(x))) else character(0L)
  }

  x
}

.eval_fit_data_arg <- function(expr, data, env, what) {
  if (is.character(expr) && length(expr) == 1L && expr %in% names(data)) {
    return(data[[expr]])
  }
  if (is.symbol(expr)) {
    name <- as.character(expr)
    if (name %in% names(data)) {
      return(data[[name]])
    }
  }

  value <- eval(expr, envir = data, enclos = env)
  if (is.character(value) && length(value) == 1L && value %in% names(data)) {
    return(data[[value]])
  }
  value
}

.resolve_fit_data_vector <- function(expr, data, env, what) {
  value <- .eval_fit_data_arg(expr = expr, data = data, env = env, what = what)
  if (is.data.frame(value)) {
    if (ncol(value) != 1L) {
      stop("`", what, "` must resolve to one column.")
    }
    value <- value[[1L]]
  }
  value
}

.resolve_fit_design_arg <- function(expr, data, env, what, allow_formula = TRUE) {
  value <- .eval_fit_data_arg(expr = expr, data = data, env = env, what = what)
  if (is.null(value)) {
    return(NULL)
  }
  if (inherits(value, "formula")) {
    if (!isTRUE(allow_formula)) {
      stop("`", what, "` must be a design matrix, not a formula. Use `boundary_formula` for formulas.")
    }
    model_frame <- stats::model.frame(value, data = data, na.action = stats::na.fail)
    return(stats::model.matrix(value, model_frame))
  }
  if (is.character(value) && all(value %in% names(data))) {
    return(as.matrix(data[value]))
  }
  value
}

.resolve_fit_formula_arg <- function(expr, data = NULL, env, what) {
  value <- if (is.null(data)) eval(expr, envir = env) else eval(expr, envir = data, enclos = env)
  if (is.null(value)) {
    return(NULL)
  }
  if (is.character(value) && length(value) == 1L) {
    value <- stats::as.formula(value, env = env)
  }
  if (!inherits(value, "formula") || length(value) != 2L) {
    stop("`", what, "` must be a one-sided formula.")
  }
  if (is.null(environment(value))) {
    environment(value) <- env
  }

  if (is.null(data) && !length(all.vars(value))) {
    return(NULL)
  }

  model_frame <- if (is.null(data)) stats::model.frame(value, na.action = stats::na.fail) else stats::model.frame(value, data = data, na.action = stats::na.fail)
  stats::model.matrix(value, model_frame)
}

.coerce_boundary_design <- function(boundary_X, curve, n_curves) {
  if (is.null(boundary_X)) {
    out <- matrix(1, nrow = n_curves, ncol = 1L)
    colnames(out) <- "(Intercept)"
    return(out)
  }

  if (is.vector(boundary_X) && !is.list(boundary_X)) {
    boundary_X <- matrix(boundary_X, ncol = 1L)
  } else {
    boundary_X <- as.matrix(boundary_X)
  }

  storage.mode(boundary_X) <- "double"
  if (anyNA(boundary_X) || any(!is.finite(boundary_X))) {
    stop("`boundary_X` must contain only finite numeric values.")
  }

  if (nrow(boundary_X) == n_curves) {
    if (is.null(colnames(boundary_X))) {
      colnames(boundary_X) <- if (ncol(boundary_X)) paste0("b", seq_len(ncol(boundary_X))) else character(0L)
    }
    return(boundary_X)
  }

  if (nrow(boundary_X) == length(curve)) {
    curve_rows <- split(seq_len(length(curve)), curve)
    collapsed <- matrix(NA_real_, nrow = n_curves, ncol = ncol(boundary_X))
    for (i in seq_along(curve_rows)) {
      idx <- curve_rows[[i]]
      ref <- boundary_X[idx[1], , drop = FALSE]
      ref_expanded <- ref[rep(1L, length(idx)), , drop = FALSE]
      same <- apply(abs(boundary_X[idx, , drop = FALSE] - ref_expanded) < 1e-10, 2, all)
      if (!all(same)) {
        stop("Observation-level `boundary_X` must be constant within each curve.")
      }
      collapsed[i, ] <- ref
    }
    colnames(collapsed) <- colnames(boundary_X)
    if (is.null(colnames(collapsed))) {
      colnames(collapsed) <- if (ncol(collapsed)) paste0("b", seq_len(ncol(collapsed))) else character(0L)
    }
    return(collapsed)
  }

  stop("`boundary_X` must have either one row per curve or one row per observation.")
}

.prepare_boundary_value_data <- function(
  y,
  time,
  curve,
  X = NULL,
  Z = NULL,
  boundary_X = NULL,
  intercept_default = TRUE
) {
  y <- as.numeric(y)
  time <- as.numeric(time)

  if (!length(y) || length(time) != length(y) || length(curve) != length(y)) {
    stop("`y`, `time`, and `curve` must have the same positive length.")
  }
  if (anyNA(y) || anyNA(time) || anyNA(curve)) {
    stop("`y`, `time`, and `curve` must not contain missing values.")
  }
  if (any(!is.finite(y)) || any(!is.finite(time))) {
    stop("`y` and `time` must be finite.")
  }

  curve <- factor(curve)
  ord <- order(curve, time)
  y <- y[ord]
  time <- time[ord]
  curve <- factor(curve[ord], levels = levels(curve))

  counts <- table(curve)
  if (length(unique(counts)) != 1L) {
    stop("All curves must have the same number of observations.")
  }

  curve_index <- split(seq_along(y), curve)
  time_by_curve <- lapply(curve_index, function(idx) time[idx])
  time_grid <- unname(time_by_curve[[1]])

  same_grid <- vapply(
    time_by_curve,
    function(x) isTRUE(all.equal(unname(x), time_grid, tolerance = 1e-10)),
    logical(1)
  )
  if (!all(same_grid)) {
    stop("All curves must share the same time grid.")
  }

  diffs <- diff(time_grid)
  if (!all(abs(diffs - diffs[1]) < 1e-10)) {
    stop("The common time grid must be equidistant.")
  }

  X <- .coerce_design_matrix(X, n_rows = length(y), what = "X", intercept_default = intercept_default)
  Z <- .coerce_design_matrix(Z, n_rows = length(y), what = "Z", intercept_default = FALSE)
  if (!is.null(boundary_X)) {
    boundary_X_check <- if (is.vector(boundary_X) && !is.list(boundary_X)) matrix(boundary_X, ncol = 1L) else as.matrix(boundary_X)
    if (nrow(boundary_X_check) == length(y)) {
      boundary_X <- boundary_X_check[ord, , drop = FALSE]
    }
  }
  boundary_X <- .coerce_boundary_design(boundary_X, curve = curve, n_curves = nlevels(curve))

  X <- X[ord, , drop = FALSE]
  Z <- Z[ord, , drop = FALSE]

  list(
    y = y,
    time = time,
    curve = curve,
    curve_levels = levels(curve),
    curve_index = curve_index,
    n_curves = nlevels(curve),
    n_per_curve = as.integer(unique(counts)),
    time_grid = as.numeric(time_grid),
    grid_spacing = diffs[1],
    X = X,
    Z = Z,
    boundary_X = boundary_X,
    order = ord
  )
}
