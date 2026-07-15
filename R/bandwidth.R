#' Select a VCMoE bandwidth by K-fold cross-validation
#'
#' @param formula A formula of the form `y ~ expert_terms | gating_terms`.
#' @param data A data frame.
#' @param u Continuous index column name or numeric vector.
#' @param k Number of mixture components. Values from 2 through 10 are
#'   supported.
#' @param family Model family. `"gaussian"`, `"binomial"`, and
#'   `"negative-binomial"` are supported.
#' @param bandwidth_grid Candidate bandwidth values. If `NULL`, uses multiples
#'   of the default bandwidth.
#' @param folds Number of random cross-validation folds.
#' @param u_grid Grid where coefficient functions are estimated.
#' @param control Named list passed to the selected fitting engine.
#' @param label Label strategy passed to the selected fitting engine.
#' @param u_scale `u` scaling strategy passed to the selected fitting engine.
#' @param parameterization Estimator convention passed to the selected fitting
#'   engine.
#' @param seed Optional random seed for fold assignment and, when
#'   `control$seed` is absent, deterministic CV refits.
#' @param refit Whether to refit the final model on all data using the selected
#'   bandwidth.
#' @param engine Fitting engine used for every cross-validation fold and the
#'   optional final refit. The default `"local_grid_em"` preserves the 0.1.0
#'   behavior; `"joint_path_em"` uses the joint-path engine of `vcmoe_fit()`.
#' @return An object of class `vcmoe_bandwidth_selection`.
#' @export
vcmoe_select_bandwidth <- function(formula, data, u, k = 2L, family = "gaussian",
                                   bandwidth_grid = NULL, folds = 5L,
                                   u_grid = NULL, control = list(),
                                   label = "align",
                                   parameterization = "a1_epanechnikov_scaled",
                                   u_scale = c("unit", "none"),
                                   seed = NULL,
                                   refit = TRUE,
                                   engine = c("local_grid_em", "joint_path_em")) {
  family <- match.arg(family, c("gaussian", "binomial", "negative-binomial"))
  label <- match.arg(label, c("align", "global", "greedy"))
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  engine <- match.arg(engine)
  k <- as.integer(k)
  if (k < 2L || k > 10L) {
    stop("`k` must be between 2 and 10.", call. = FALSE)
  }
  if (!is.logical(refit) || length(refit) != 1L || is.na(refit)) {
    stop("`refit` must be TRUE or FALSE.", call. = FALSE)
  }

  info <- .bandwidth_selection_data_info(formula, data, u, family, u_scale)
  n_complete <- length(info$complete_rows)
  folds <- .validate_cv_folds(folds, n_complete)
  default_bandwidth <- .default_bandwidth(info$u_values)
  candidates <- .validate_bandwidth_grid(bandwidth_grid, default_bandwidth)
  u_grid_values <- .validate_selection_u_grid(u_grid, info$u_values)
  fold_id <- .make_cv_folds(n_complete, folds, seed)
  if (identical(engine, "joint_path_em")) {
    .validate_joint_path_cv_grid(
      info$u_values,
      u_grid_values,
      fold_id,
      .vcmoe_default_control(control),
      warn_unique = !isTRUE(refit)
    )
    warning(
      sprintf(
        "Joint-path bandwidth selection will run %d candidate(s) x %d fold(s), plus the optional final refit; runtime scales with every joint-path fit.",
        length(candidates),
        folds
      ),
      call. = FALSE
    )
  }

  fold_results <- vector("list", length(candidates) * folds)
  result_id <- 1L
  for (bandwidth_id in seq_along(candidates)) {
    for (fold in seq_len(folds)) {
      fold_results[[result_id]] <- .score_bandwidth_fold(
        formula = formula,
        data = data,
        u = u,
        family = family,
        k = k,
        bandwidth = candidates[[bandwidth_id]],
        bandwidth_id = bandwidth_id,
        fold = fold,
        fold_id = fold_id,
        complete_rows = info$complete_rows,
        u_grid = u_grid_values,
        control = control,
        label = label,
        u_scale = u_scale,
        parameterization = parameterization,
        seed = seed,
        engine = engine
      )
      result_id <- result_id + 1L
    }
  }
  cv_details <- do.call(rbind, fold_results)
  cv_summary <- .summarise_bandwidth_cv(cv_details, candidates, default_bandwidth)
  if (!any(is.finite(cv_summary$total_loglik))) {
    stop("All candidate bandwidths failed during cross-validation.", call. = FALSE)
  }
  best_row <- which(cv_summary$selected)[[1L]]
  best_bandwidth <- cv_summary$bandwidth[[best_row]]

  final_fit <- NULL
  if (isTRUE(refit)) {
    final_control <- .selection_fit_control(control, seed, bandwidth_id = 999L, fold = 999L)
    final_fit <- .selection_fit_engine(
      engine = engine,
      formula = formula,
      data = data,
      u = u,
      k = k,
      family = family,
      bandwidth = best_bandwidth,
      u_grid = u_grid_values,
      control = final_control,
      label = label,
      u_scale = u_scale,
      parameterization = parameterization
    )
  }

  out <- list(
    best_bandwidth = best_bandwidth,
    cv_summary = cv_summary,
    cv_details = cv_details,
    cv_folds = data.frame(
      row = info$complete_rows,
      fold = fold_id,
      stringsAsFactors = FALSE
    ),
    fit = final_fit,
    settings = list(
      family = family,
      k = k,
      u_grid = u_grid_values,
      bandwidth_grid = candidates,
      folds = folds,
      seed = seed,
      criterion = "heldout_predictive_loglik",
      default_bandwidth = default_bandwidth,
      refit = refit,
      engine = engine,
      engine_id = engine,
      u_scale = u_scale,
      u_scaling = info$u_scaling,
      parameterization = parameterization
    )
  )
  class(out) <- "vcmoe_bandwidth_selection"
  out
}

.selection_fit_engine <- function(engine, ...) {
  fit_args <- list(...)
  switch(
    engine,
    local_grid_em = do.call(vcmoe_fit, fit_args),
    joint_path_em = {
      fit_args$engine <- "joint_path_em"
      do.call(vcmoe_fit, fit_args)
    },
    stop("Unknown bandwidth-selection engine: ", engine, call. = FALSE)
  )
}

.validate_joint_path_cv_grid <- function(u_values, u_grid, fold_id, control,
                                         warn_unique = TRUE) {
  training_n <- length(u_values) - max(tabulate(fold_id))
  .joint_path_cost_guard(max(training_n, 1L), length(u_grid), control)
  if (isTRUE(warn_unique)) {
    .joint_path_unique_grid_warning(u_values, u_grid)
  }
  invisible(TRUE)
}

.bandwidth_selection_data_info <- function(formula, data, u, family, u_scale = "unit") {
  pieces <- .split_vcmoe_formula(formula)
  u_info <- .extract_u(u, data)
  u_values_all <- as.numeric(u_info$values)
  expert_frame_all <- stats::model.frame(
    pieces$expert_formula,
    data = data,
    na.action = stats::na.pass
  )
  gating_frame_all <- stats::model.frame(
    pieces$gating_formula,
    data = data,
    na.action = stats::na.pass
  )
  complete <- .complete_cases_frame(expert_frame_all, nrow(data)) &
    .complete_cases_frame(gating_frame_all, nrow(data)) &
    is.finite(u_values_all)
  if (!any(complete)) {
    stop("No complete rows remain after checking response, covariates, and `u`.", call. = FALSE)
  }

  expert_frame <- expert_frame_all[complete, , drop = FALSE]
  .parse_vcmoe_response(family, stats::model.response(expert_frame), pieces$response)
  u_scaling <- .vcmoe_make_u_scaling(u_values_all[complete], u_scale)
  list(
    complete_rows = which(complete),
    u_values = .vcmoe_scale_u(u_values_all[complete], u_scaling),
    u_values_original = u_values_all[complete],
    u_scaling = u_scaling
  )
}

.validate_cv_folds <- function(folds, n_complete) {
  if (!is.numeric(folds) || length(folds) != 1L || !is.finite(folds) ||
      abs(folds - round(folds)) > 1e-8) {
    stop("`folds` must be a single finite whole number.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (folds < 2L || folds > n_complete) {
    stop("`folds` must be between 2 and the number of complete rows.", call. = FALSE)
  }
  folds
}

.default_bandwidth_grid <- function(u) {
  .default_bandwidth(u) * c(0.5, 0.75, 1, 1.25, 1.5, 2)
}

.validate_bandwidth_grid <- function(bandwidth_grid, default_bandwidth) {
  if (is.null(bandwidth_grid)) {
    bandwidth_grid <- default_bandwidth * c(0.5, 0.75, 1, 1.25, 1.5, 2)
  }
  if (!is.numeric(bandwidth_grid) || !length(bandwidth_grid) ||
      any(!is.finite(bandwidth_grid)) || any(bandwidth_grid <= 0)) {
    stop("`bandwidth_grid` must contain positive finite numeric values.", call. = FALSE)
  }
  sort(unique(as.numeric(bandwidth_grid)))
}

.validate_selection_u_grid <- function(u_grid, u_values) {
  if (is.null(u_grid)) {
    return(.default_u_grid(u_values))
  }
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (!length(u_grid) || any(!is.finite(u_grid))) {
    stop("`u_grid` must contain finite numeric values.", call. = FALSE)
  }
  u_grid
}

.make_cv_folds <- function(n, folds, seed) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sample(rep(seq_len(folds), length.out = n))
}

.subset_u_for_rows <- function(u, rows) {
  if (is.character(u) && length(u) == 1L) {
    return(u)
  }
  u[rows]
}

.selection_fit_control <- function(control, seed, bandwidth_id, fold) {
  control <- control %||% list()
  if (is.null(control$seed) && !is.null(seed)) {
    control$seed <- as.integer(seed + bandwidth_id * 1000L + fold)
  }
  control
}

.score_bandwidth_fold <- function(formula, data, u, family, k, bandwidth,
                                  bandwidth_id, fold, fold_id, complete_rows,
                                  u_grid, control, label, u_scale, parameterization, seed,
                                  engine) {
  validation_rows <- complete_rows[fold_id == fold]
  training_rows <- complete_rows[fold_id != fold]
  start_time <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    suppressWarnings(.selection_fit_engine(
      engine = engine,
      formula = formula,
      data = data[training_rows, , drop = FALSE],
      u = .subset_u_for_rows(u, training_rows),
      k = k,
      family = family,
      bandwidth = bandwidth,
      u_grid = u_grid,
      control = .selection_fit_control(control, seed, bandwidth_id, fold),
      label = label,
      u_scale = u_scale,
      parameterization = parameterization
    )),
    error = function(e) e
  )
  runtime <- proc.time()[["elapsed"]] - start_time

  if (inherits(fit, "error")) {
    return(data.frame(
      bandwidth = bandwidth,
      bandwidth_id = bandwidth_id,
      fold = fold,
      fit_status = "failed",
      error_message = conditionMessage(fit),
      heldout_loglik = NA_real_,
      validation_n = length(validation_rows),
      converged_grid_points = NA_integer_,
      grid_points = length(u_grid),
      ambiguity_warnings = NA_integer_,
      runtime_seconds = runtime,
      engine = engine,
      engine_id = engine,
      stringsAsFactors = FALSE
    ))
  }

  loglik <- tryCatch(
    .vcmoe_predictive_loglik(
      fit,
      newdata = data[validation_rows, , drop = FALSE],
      u = .subset_u_for_rows(u, validation_rows)
    ),
    error = function(e) e
  )
  if (inherits(loglik, "error") || !is.finite(loglik)) {
    return(data.frame(
      bandwidth = bandwidth,
      bandwidth_id = bandwidth_id,
      fold = fold,
      fit_status = "failed",
      error_message = if (inherits(loglik, "error")) conditionMessage(loglik) else "Non-finite held-out log-likelihood.",
      heldout_loglik = NA_real_,
      validation_n = length(validation_rows),
      converged_grid_points = sum(fit$diagnostics$converged),
      grid_points = length(fit$u_grid),
      ambiguity_warnings = length(fit$diagnostics$warnings),
      runtime_seconds = runtime,
      engine = engine,
      engine_id = fit$engine_id %||% "local_grid_em",
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    bandwidth = bandwidth,
    bandwidth_id = bandwidth_id,
    fold = fold,
    fit_status = "ok",
    error_message = NA_character_,
    heldout_loglik = loglik,
    validation_n = length(validation_rows),
    converged_grid_points = sum(fit$diagnostics$converged),
    grid_points = length(fit$u_grid),
    ambiguity_warnings = length(fit$diagnostics$warnings),
    runtime_seconds = runtime,
    engine = engine,
    engine_id = fit$engine_id %||% "local_grid_em",
    stringsAsFactors = FALSE
  )
}

.vcmoe_predictive_loglik <- function(object, newdata, u = NULL) {
  pieces <- .predict_components(object, newdata = newdata, u = u)
  response <- .prediction_response(object, newdata)
  if (is.null(response)) {
    stop("Validation data must contain the response columns used by the fit.", call. = FALSE)
  }

  n <- nrow(pieces$prior)
  log_terms <- matrix(NA_real_, nrow = n, ncol = object$k)
  for (component in seq_len(object$k)) {
    if (identical(object$family, "binomial")) {
      expert_loglik <- stats::dbinom(
        response$success,
        size = response$trials,
        prob = pmin(pmax(pieces$mean[, component], 1e-12), 1 - 1e-12),
        log = TRUE
      )
    } else if (identical(object$family, "negative-binomial")) {
      expert_loglik <- stats::dnbinom(
        response$y,
        size = pieces$theta[, component],
        mu = pmax(pieces$mean[, component], 1e-12),
        log = TRUE
      )
    } else {
      expert_loglik <- stats::dnorm(
        response$y,
        mean = pieces$mean[, component],
        sd = pmax(pieces$sigma[, component], 1e-12),
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(pieces$prior[, component], 1e-12)) + expert_loglik
  }
  sum(.row_log_sum_exp(log_terms))
}

.summarise_bandwidth_cv <- function(cv_details, candidates, default_bandwidth) {
  rows <- lapply(seq_along(candidates), function(i) {
    detail <- cv_details[cv_details$bandwidth_id == i, , drop = FALSE]
    ok <- identical(any(detail$fit_status == "ok"), TRUE)
    total_loglik <- if (ok) sum(detail$heldout_loglik, na.rm = TRUE) else -Inf
    scored_n <- sum(detail$validation_n[detail$fit_status == "ok"])
    data.frame(
      bandwidth = candidates[[i]],
      bandwidth_id = i,
      total_loglik = total_loglik,
      mean_loglik = if (scored_n > 0) total_loglik / scored_n else NA_real_,
      failed_folds = sum(detail$fit_status != "ok"),
      successful_folds = sum(detail$fit_status == "ok"),
      validation_n = sum(detail$validation_n),
      scored_n = scored_n,
      converged_grid_points = sum(detail$converged_grid_points, na.rm = TRUE),
      grid_points = sum(detail$grid_points, na.rm = TRUE),
      ambiguity_warnings = sum(detail$ambiguity_warnings, na.rm = TRUE),
      runtime_seconds = sum(detail$runtime_seconds),
      distance_to_default = abs(candidates[[i]] - default_bandwidth),
      selected = FALSE,
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  finite_rows <- which(is.finite(summary$total_loglik))
  if (length(finite_rows)) {
    no_failure_rows <- finite_rows[summary$failed_folds[finite_rows] == 0L]
    ranking_rows <- if (length(no_failure_rows)) no_failure_rows else finite_rows
    rank_order <- order(-summary$total_loglik[ranking_rows], summary$distance_to_default[ranking_rows])
    summary$selected[ranking_rows[rank_order[[1L]]]] <- TRUE
  }
  summary
}

#' @export
print.vcmoe_bandwidth_selection <- function(x, ...) {
  cat("VCMoE bandwidth selection\n")
  cat("  engine: ", x$settings$engine %||% "local_grid_em", "\n", sep = "")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  components: ", x$settings$k, "\n", sep = "")
  cat("  parameterization: ", x$settings$parameterization %||% "unknown", "\n", sep = "")
  cat("  criterion: ", x$settings$criterion, "\n", sep = "")
  cat("  candidates: ", length(x$settings$bandwidth_grid), "\n", sep = "")
  cat("  folds: ", x$settings$folds, "\n", sep = "")
  cat("  best bandwidth: ", signif(x$best_bandwidth, 4), "\n", sep = "")
  if (!is.null(x$fit)) {
    cat("  final fit converged grid points: ", sum(x$fit$diagnostics$converged), "/",
        length(x$fit$u_grid), "\n", sep = "")
  }
  invisible(x)
}
