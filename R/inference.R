.validate_bootstrap_coefficient_set <- function(coefficient_set) {
  match.arg(coefficient_set, c("expert", "gating"), several.ok = TRUE)
}

.bootstrap_base_data <- function(fit, data) {
  if (!is.data.frame(data)) {
    stop("`data` must be the original data frame used to fit `fit`.", call. = FALSE)
  }
  if (is.null(fit$rows_used) || !length(fit$rows_used) || max(fit$rows_used) > nrow(data)) {
    stop("`data` does not appear to match the rows used by `fit`.", call. = FALSE)
  }
  data[fit$rows_used, , drop = FALSE]
}

.bootstrap_u_info <- function(fit, data, base_data, u) {
  if (is.null(u)) {
    if (!is.null(fit$u_name) && fit$u_name %in% names(base_data)) {
      return(list(values = as.numeric(base_data[[fit$u_name]]), refit = fit$u_name))
    }
    if (!is.null(fit$fitted$u_original) && length(fit$fitted$u_original) == nrow(base_data)) {
      return(list(values = as.numeric(fit$fitted$u_original), refit = as.numeric(fit$fitted$u_original)))
    }
    if (!is.null(fit$fitted$u) && length(fit$fitted$u) == nrow(base_data)) {
      return(list(values = as.numeric(fit$fitted$u), refit = as.numeric(fit$fitted$u)))
    }
    stop("Provide `u`; the original fit did not store a reusable `u` column.", call. = FALSE)
  }
  if (is.character(u) && length(u) == 1L) {
    if (!u %in% names(base_data)) {
      stop("`u` was provided as a name but is not present in `data`.", call. = FALSE)
    }
    return(list(values = as.numeric(base_data[[u]]), refit = u))
  }
  if (length(u) == nrow(data)) {
    values <- as.numeric(u[fit$rows_used])
  } else if (length(u) == nrow(base_data)) {
    values <- as.numeric(u)
  } else {
    stop("`u` must have one value per original data row or one value per used row.", call. = FALSE)
  }
  list(values = values, refit = values)
}

.vcmoe_refit_label <- function(fit) {
  label <- fit$label %||% fit$diagnostics$alignment_method %||% "align"
  label <- as.character(label[[1L]])
  if (identical(label, "sequential")) {
    return("align")
  }
  if (!label %in% c("align", "global", "greedy")) {
    return("align")
  }
  label
}

.sample_component <- function(prior) {
  apply(prior, 1L, function(p) {
    sample.int(ncol(prior), size = 1L, prob = p)
  })
}

.simulate_bootstrap_response <- function(fit, base_data, u_values) {
  pieces <- .predict_components(fit, newdata = base_data, u = u_values)
  n <- nrow(base_data)
  component <- .sample_component(pieces$prior)
  component_index <- cbind(seq_len(n), component)
  out <- base_data

  if (identical(fit$family, "gaussian")) {
    if (is.null(fit$response) || !fit$response %in% names(out)) {
      stop("Gaussian bootstrap requires the original response column in `data`.", call. = FALSE)
    }
    out[[fit$response]] <- stats::rnorm(
      n,
      mean = pieces$mean[component_index],
      sd = pmax(pieces$sigma[component_index], 1e-8)
    )
    return(out)
  }

  if (identical(fit$family, "negative-binomial")) {
    if (is.null(fit$response) || !fit$response %in% names(out)) {
      stop("Negative-Binomial bootstrap requires the original response column in `data`.", call. = FALSE)
    }
    out[[fit$response]] <- stats::rnbinom(
      n,
      size = pmax(pieces$theta[component_index], 1e-8),
      mu = pmax(pieces$mean[component_index], 1e-12)
    )
    return(out)
  }

  info <- fit$response_info
  prob <- pmin(pmax(pieces$mean[component_index], 1e-12), 1 - 1e-12)
  if (identical(info$type, "bernoulli")) {
    if (is.null(info$response) || !info$response %in% names(out)) {
      stop("Bernoulli bootstrap requires the original response column in `data`.", call. = FALSE)
    }
    out[[info$response]] <- stats::rbinom(n, size = 1L, prob = prob)
    return(out)
  }

  if (!identical(info$type, "grouped") || is.null(info$success) || is.null(info$failure) ||
      !info$success %in% names(out) || !info$failure %in% names(out)) {
    stop("Grouped Binomial bootstrap requires success and failure columns in `data`.", call. = FALSE)
  }
  trials <- as.numeric(out[[info$success]]) + as.numeric(out[[info$failure]])
  success <- stats::rbinom(n, size = trials, prob = prob)
  out[[info$success]] <- success
  out[[info$failure]] <- trials - success
  out
}

.bootstrap_match_cost <- function(reference, candidate, permutation, coefficient_set) {
  cost <- 0
  for (set in coefficient_set) {
    ref <- reference$coefficients[[set]]
    cur <- candidate$coefficients[[set]][, permutation, , drop = FALSE]
    cost <- cost + mean((cur - ref)^2)
    slope_set <- paste0(set, "_slope")
    if (!is.null(reference$coefficients[[slope_set]]) && !is.null(candidate$coefficients[[slope_set]])) {
      ref_slope <- reference$coefficients[[slope_set]]
      cur_slope <- candidate$coefficients[[slope_set]][, permutation, , drop = FALSE]
      cost <- cost + 0.25 * mean((cur_slope - ref_slope)^2)
    }
  }
  if (identical(reference$family, "gaussian") && !is.null(reference$coefficients$sigma) &&
      !is.null(candidate$coefficients$sigma)) {
    cost <- cost + 0.10 * mean((candidate$coefficients$sigma[, permutation, drop = FALSE] -
      reference$coefficients$sigma)^2)
    if (!is.null(reference$coefficients$sigma_slope) && !is.null(candidate$coefficients$sigma_slope)) {
      cost <- cost + 0.05 * mean((candidate$coefficients$sigma_slope[, permutation, drop = FALSE] -
        reference$coefficients$sigma_slope)^2)
    }
  }
  if (identical(reference$family, "negative-binomial") && !is.null(reference$coefficients$theta) &&
      !is.null(candidate$coefficients$theta)) {
    cost <- cost + 0.10 * mean((log(candidate$coefficients$theta[, permutation, drop = FALSE]) -
      log(reference$coefficients$theta))^2)
  }
  cost
}

.bootstrap_component_match_cost <- function(reference, candidate, reference_component,
                                           candidate_component, coefficient_set) {
  cost <- 0
  for (set in coefficient_set) {
    ref <- reference$coefficients[[set]][, reference_component, , drop = FALSE]
    cur <- candidate$coefficients[[set]][, candidate_component, , drop = FALSE]
    cost <- cost + mean((cur - ref)^2)
    slope_set <- paste0(set, "_slope")
    if (!is.null(reference$coefficients[[slope_set]]) && !is.null(candidate$coefficients[[slope_set]])) {
      ref_slope <- reference$coefficients[[slope_set]][, reference_component, , drop = FALSE]
      cur_slope <- candidate$coefficients[[slope_set]][, candidate_component, , drop = FALSE]
      cost <- cost + 0.25 * mean((cur_slope - ref_slope)^2)
    }
  }
  if (identical(reference$family, "gaussian") && !is.null(reference$coefficients$sigma) &&
      !is.null(candidate$coefficients$sigma)) {
    cost <- cost + 0.10 * mean((candidate$coefficients$sigma[, candidate_component, drop = FALSE] -
      reference$coefficients$sigma[, reference_component, drop = FALSE])^2)
    if (!is.null(reference$coefficients$sigma_slope) && !is.null(candidate$coefficients$sigma_slope)) {
      cost <- cost + 0.05 * mean((candidate$coefficients$sigma_slope[, candidate_component, drop = FALSE] -
        reference$coefficients$sigma_slope[, reference_component, drop = FALSE])^2)
    }
  }
  if (identical(reference$family, "negative-binomial") && !is.null(reference$coefficients$theta) &&
      !is.null(candidate$coefficients$theta)) {
    cost <- cost + 0.10 * mean((log(candidate$coefficients$theta[, candidate_component, drop = FALSE]) -
      log(reference$coefficients$theta[, reference_component, drop = FALSE]))^2)
  }
  cost
}

.bootstrap_reference_permutation <- function(reference, candidate, coefficient_set) {
  if (reference$k <= 6L) {
    perms <- .all_permutations(reference$k)
    scores <- apply(perms, 1L, function(perm) {
      .bootstrap_match_cost(reference, candidate, perm, coefficient_set)
    })
    ord <- order(scores)
    best <- perms[ord[[1L]], ]
    best_score <- scores[ord[[1L]]]
    second_score <- if (length(ord) > 1L) scores[ord[[2L]]] else Inf
  } else {
    costs <- matrix(NA_real_, nrow = reference$k, ncol = reference$k)
    for (reference_component in seq_len(reference$k)) {
      for (candidate_component in seq_len(reference$k)) {
        costs[candidate_component, reference_component] <- .bootstrap_component_match_cost(
          reference,
          candidate,
          reference_component,
          candidate_component,
          coefficient_set
        )
      }
    }
    assignment <- .assignment_permutation(costs)
    best <- assignment$permutation
    best_score <- assignment$score
    second_score <- assignment$second_score
  }
  margin <- second_score - best_score
  relative_margin <- margin / (abs(best_score) + 1e-8)
  list(
    permutation = best,
    cost = best_score,
    margin = margin,
    relative_margin = relative_margin,
    ambiguous = is.finite(margin) && (margin < 1e-6 ||
      (is.finite(relative_margin) && relative_margin < 0.01))
  )
}

.bootstrap_aligned_coefficients <- function(fit, permutation, coefficient_set) {
  out <- list()
  for (set in coefficient_set) {
    out[[set]] <- fit$coefficients[[set]][, permutation, , drop = FALSE]
  }
  out
}

.stack_bootstrap_replicates <- function(successful, coefficient_set, reference_fit) {
  out <- list()
  for (set in coefficient_set) {
    template <- reference_fit$coefficients[[set]]
    arr <- array(
      NA_real_,
      dim = c(dim(template), length(successful)),
      dimnames = c(dimnames(template), list(replicate = paste0("replicate", seq_along(successful))))
    )
    for (i in seq_along(successful)) {
      arr[, , , i] <- successful[[i]]$coefficients[[set]]
    }
    out[[set]] <- arr
  }
  out
}

.bootstrap_refit_control <- function(fit_control, control, seed, replicate_id) {
  fit_control <- fit_control %||% list()
  user_control <- control %||% list()
  if (!"seed" %in% names(user_control)) {
    fit_control$seed <- NULL
  }
  out <- utils::modifyList(fit_control, user_control)
  if (is.null(out$seed) && !is.null(seed)) {
    out$seed <- as.integer(seed + 10000L + replicate_id)
  }
  out
}

.vcmoe_bootstrap_refit <- function(reference_fit, fit_args) {
  engine_id <- reference_fit$engine_id %||% "local_grid_em"
  if (identical(engine_id, "joint_path_em")) {
    fit_args$engine <- "joint_path_em"
  }
  do.call(vcmoe_fit, fit_args)
}

#' Parametric bootstrap inference for a VCMoE fit
#'
#' @param fit A `vcmoe` fit with `k = 2:10`.
#' @param data Original data frame used to fit `fit`.
#' @param u Optional original `u` values or column name.
#' @param B Number of bootstrap replicates.
#' @param coefficient_set Coefficient sets to store.
#' @param seed Optional random seed.
#' @param control Control overrides for bootstrap refits.
#' @param min_successful Minimum successful replicates for reliable inference.
#' @param keep_fits Whether to store successful bootstrap fit objects.
#' @param verbose Whether to message progress.
#' @details Bootstrap refits preserve the reference fitting engine. A
#'   joint-path reference is therefore refitted with
#'   `vcmoe_fit(..., engine = "joint_path_em")` rather than silently falling
#'   back to local-grid EM.
#' @return An object of class `vcmoe_bootstrap`.
#' @export
vcmoe_bootstrap <- function(fit, data, u = NULL, B = 200L,
                            coefficient_set = c("expert", "gating"),
                            seed = NULL, control = list(),
                            min_successful = max(20L, ceiling(0.5 * B)),
                            keep_fits = FALSE, verbose = FALSE) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 2L || abs(B - round(B)) > 1e-8) {
    stop("`B` must be a whole number of at least 2.", call. = FALSE)
  }
  B <- as.integer(B)
  if (!is.numeric(min_successful) || length(min_successful) != 1L ||
      !is.finite(min_successful) || min_successful < 2L) {
    stop("`min_successful` must be at least 2.", call. = FALSE)
  }
  min_successful <- as.integer(ceiling(min_successful))
  coefficient_set <- .validate_bootstrap_coefficient_set(coefficient_set)
  base_data <- .bootstrap_base_data(fit, data)
  u_info <- .bootstrap_u_info(fit, data, base_data, u)

  warnings <- character(0L)
  if (any(!fit$diagnostics$converged)) {
    warnings <- c(warnings, "Reference fit has non-converged grid points.")
  }
  if (any(fit$diagnostics$ambiguous %||% FALSE)) {
    warnings <- c(warnings, "Reference fit has ambiguous label alignment.")
  }

  successful <- list()
  stored_fits <- list()
  summary_rows <- vector("list", B)
  alignment_rows <- vector("list", B)

  for (replicate_id in seq_len(B)) {
    if (!is.null(seed)) {
      set.seed(as.integer(seed + replicate_id))
    }
    if (isTRUE(verbose)) {
      message("Bootstrap replicate ", replicate_id, " / ", B)
    }
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      boot_data <- .simulate_bootstrap_response(fit, base_data, u_info$values)
      fit_args <- list(
        formula = fit$formula,
        data = boot_data,
        u = u_info$refit,
        k = fit$k,
        family = fit$family,
        bandwidth = fit$bandwidth,
        u_grid = fit$u_grid,
        control = .bootstrap_refit_control(fit$control, control, seed, replicate_id),
        label = .vcmoe_refit_label(fit),
        u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
        parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id %||% "a1_epanechnikov_scaled"
      )
      boot_fit <- suppressWarnings(.vcmoe_bootstrap_refit(fit, fit_args))
      match <- .bootstrap_reference_permutation(fit, boot_fit, coefficient_set)
      aligned <- .bootstrap_aligned_coefficients(boot_fit, match$permutation, coefficient_set)
      list(fit = boot_fit, coefficients = aligned, match = match)
    }, error = function(e) e)
    runtime <- proc.time()[["elapsed"]] - start_time

    if (inherits(result, "error")) {
      summary_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = "failed",
        successful_index = NA_integer_,
        converged_grid_points = NA_integer_,
        grid_points = length(fit$u_grid),
        ambiguity_warnings = NA_integer_,
        runtime_seconds = runtime,
        error_message = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      alignment_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        permutation = NA_character_,
        cost = NA_real_,
        margin = NA_real_,
        relative_margin = NA_real_,
        ambiguous = NA,
        stringsAsFactors = FALSE
      )
    } else {
      successful_index <- length(successful) + 1L
      successful[[successful_index]] <- list(coefficients = result$coefficients)
      if (isTRUE(keep_fits)) {
        stored_fits[[successful_index]] <- result$fit
      }
      summary_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = "ok",
        successful_index = successful_index,
        converged_grid_points = sum(result$fit$diagnostics$converged),
        grid_points = length(result$fit$u_grid),
        ambiguity_warnings = length(result$fit$diagnostics$warnings),
        runtime_seconds = runtime,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
      alignment_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        permutation = paste(result$match$permutation, collapse = ","),
        cost = result$match$cost,
        margin = result$match$margin,
        relative_margin = result$match$relative_margin,
        ambiguous = result$match$ambiguous,
        stringsAsFactors = FALSE
      )
    }
  }

  replicate_summary <- do.call(rbind, summary_rows)
  alignment_summary <- do.call(rbind, alignment_rows)
  n_successful <- length(successful)
  if (n_successful < 2L) {
    stop("Fewer than 2 bootstrap replicates succeeded; inference cannot be summarized.", call. = FALSE)
  }
  if (n_successful < min_successful) {
    warnings <- c(
      warnings,
      sprintf("Only %d bootstrap replicates succeeded; requested minimum was %d.",
              n_successful, min_successful)
    )
  }
  if (any(replicate_summary$status != "ok")) {
    warnings <- c(warnings, sprintf("%d bootstrap replicate(s) failed.", sum(replicate_summary$status != "ok")))
  }
  if (any(alignment_summary$ambiguous %in% TRUE, na.rm = TRUE)) {
    warnings <- c(warnings, "Some bootstrap-to-reference label matches were ambiguous.")
  }

  out <- list(
    fit = fit,
    replicates = .stack_bootstrap_replicates(successful, coefficient_set, fit),
    replicate_summary = replicate_summary,
    alignment_summary = alignment_summary,
    settings = list(
      B = B,
      seed = seed,
      bandwidth = fit$bandwidth,
      family = fit$family,
      k = fit$k,
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id %||% "a1_epanechnikov_scaled",
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      u_scaling = fit$u_scaling,
      u_grid = fit$u_grid,
      coefficient_set = coefficient_set,
      min_successful = min_successful,
      keep_fits = keep_fits,
      engine_id = fit$engine_id %||% "local_grid_em"
    ),
    warnings = unique(warnings),
    fits = if (isTRUE(keep_fits)) stored_fits else NULL
  )
  class(out) <- "vcmoe_bootstrap"
  if (length(out$warnings)) {
    warning(paste(out$warnings, collapse = " "), call. = FALSE)
  }
  out
}

.bootstrap_confint_one <- function(object, set, level, type) {
  boot <- object$replicates[[set]]
  if (is.null(boot)) {
    stop("Bootstrap object does not contain coefficient set `", set, "`.", call. = FALSE)
  }
  estimate <- object$fit$coefficients[[set]]
  n_successful <- dim(boot)[4L]
  alpha <- 1 - level
  se <- apply(boot, 1:3, stats::sd, na.rm = TRUE)
  lower <- upper <- estimate

  for (grid_id in seq_len(dim(estimate)[1L])) {
    for (component in seq_len(dim(estimate)[2L])) {
      for (term in seq_len(dim(estimate)[3L])) {
        values <- boot[grid_id, component, term, ]
        if (identical(type, "pointwise")) {
          qs <- stats::quantile(values, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
          lower[grid_id, component, term] <- qs[[1L]]
          upper[grid_id, component, term] <- qs[[2L]]
        }
      }
    }
  }

  if (identical(type, "simultaneous")) {
    for (component in seq_len(dim(estimate)[2L])) {
      for (term in seq_len(dim(estimate)[3L])) {
        sd_path <- pmax(se[, component, term], 1e-8)
        max_stats <- vapply(seq_len(n_successful), function(rep_id) {
          max(abs((boot[, component, term, rep_id] - estimate[, component, term]) / sd_path), na.rm = TRUE)
        }, numeric(1L))
        critical <- stats::quantile(max_stats, probs = level, na.rm = TRUE, names = FALSE)
        lower[, component, term] <- estimate[, component, term] - critical * sd_path
        upper[, component, term] <- estimate[, component, term] + critical * sd_path
      }
    }
  }

  grid <- expand.grid(
    u = object$fit$u_grid,
    component = dimnames(estimate)[[2L]],
    term = dimnames(estimate)[[3L]],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$coefficient_set <- set
  grid$estimate <- as.vector(estimate)
  grid$se <- as.vector(se)
  grid$lower <- as.vector(lower)
  grid$upper <- as.vector(upper)
  grid$type <- type
  grid$level <- level
  grid$n_successful <- n_successful
  grid[, c("coefficient_set", "term", "component", "u", "estimate", "se", "lower", "upper",
           "type", "level", "n_successful")]
}

#' Bootstrap confidence intervals for VCMoE coefficients
#'
#' @param object A `vcmoe_bootstrap` object.
#' @param parm Coefficient set to summarize.
#' @param level Confidence level.
#' @param type Interval type.
#' @param ... Unused.
#' @return A tidy data frame of coefficient intervals.
#' @export
confint.vcmoe_bootstrap <- function(object, parm = c("expert", "gating"),
                                    level = 0.95,
                                    type = c("pointwise", "simultaneous"), ...) {
  if (!inherits(object, "vcmoe_bootstrap")) {
    stop("`object` must be a VCMoE bootstrap object.", call. = FALSE)
  }
  parm <- match.arg(parm, c("expert", "gating"), several.ok = TRUE)
  type <- match.arg(type)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }
  do.call(rbind, lapply(parm, function(set) {
    .bootstrap_confint_one(object, set, level, type)
  }))
}

#' Plot bootstrap inference intervals
#'
#' @param object A `vcmoe_bootstrap` object.
#' @param coefficient_set Coefficient set to plot.
#' @param type Interval type.
#' @param level Confidence level.
#' @return A `ggplot` object.
#' @export
plot_inference <- function(object, coefficient_set = "expert",
                           type = c("pointwise", "simultaneous"),
                           level = 0.95) {
  type <- match.arg(type)
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating"))
  df <- confint(object, parm = coefficient_set, level = level, type = type)
  ggplot2::ggplot(df, ggplot2::aes(x = u, y = estimate, color = component, fill = component)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.18, color = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~term, scales = "free_y") +
    ggplot2::labs(x = "u", y = "coefficient", color = "component", fill = "component") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @export
print.vcmoe_bootstrap <- function(x, ...) {
  successful <- sum(x$replicate_summary$status == "ok")
  cat("VCMoE bootstrap inference\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  components: ", x$settings$k, "\n", sep = "")
  cat("  bootstrap replicates: ", successful, "/", x$settings$B, " successful\n", sep = "")
  cat("  coefficient sets: ", paste(x$settings$coefficient_set, collapse = ", "), "\n", sep = "")
  if (length(x$warnings)) {
    cat("  warnings: ", length(x$warnings), "\n", sep = "")
  }
  invisible(x)
}
