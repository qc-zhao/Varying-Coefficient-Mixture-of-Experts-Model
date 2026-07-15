.vcmoe_confband_filter <- function(intervals, coefficient_set) {
  keep <- rep(FALSE, nrow(intervals))
  if ("expert" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "expert"
  }
  if ("gating" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "gating_contrast"
  }
  if ("sigma" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "nuisance" & intervals$term == "log_sigma"
  }
  if ("theta" %in% coefficient_set) {
    keep <- keep | intervals$coefficient_set == "nuisance" & intervals$term == "log_theta"
  }
  intervals[keep, , drop = FALSE]
}

#' Analytic-style confidence bands for a VCMoE fit
#'
#' @param fit A `vcmoe` fit with `k = 2:10`.
#' @param data Optional original data frame. The current implementation uses
#'   the data stored in `fit$fitted`; refit with `keep_data = TRUE` if needed.
#' @param level Confidence level.
#' @param type Interval columns to expose as `lower` and `upper`.
#' @param coefficient_set Coefficient blocks to return.
#' @param strict Whether weak local fits should return blocked intervals.
#' @param control Optional development inference controls. HC0 is the only
#'   active covariance adjustment.
#' @details For `engine = "joint_path_em"`, the covariance follows the JASA
#'   observed local-likelihood asymptotic sandwich plug-in. It does not include
#'   shared-path, label-selection, or finite-grid cross-grid responsibility
#'   uncertainty. Joint-path convergence and the returned score-imbalance
#'   diagnostics should therefore be inspected.
#' @return A `vcmoe_confband` object with interval and diagnostic data frames.
#' @export
vcmoe_confband <- function(fit, data = NULL, level = 0.95,
                           type = c("pointwise", "simultaneous"),
                           coefficient_set = c("expert", "gating", "sigma", "theta"),
                           strict = TRUE,
                           control = list()) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  type <- match.arg(type)
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating", "sigma", "theta"), several.ok = TRUE)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number between 0 and 1.", call. = FALSE)
  }
  joint_path <- identical(fit$engine_id %||% "local_grid_em", "joint_path_em")
  if (joint_path) {
    warning(
      paste0(
        "Joint-path analytic-style bands use the JASA observed local-likelihood ",
        "sandwich plug-in and exclude shared-path, label-selection, and finite-grid ",
        "cross-grid responsibility uncertainty; inspect convergence and ",
        "`score_imbalance_max`."
      ),
      call. = FALSE
    )
  }
  domain_length <- fit$u_scaling$domain_length %||%
    vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  inference_control <- utils::modifyList(
    list(
      simultaneous_method = "analytic_epanechnikov_path",
      covariance_adjustment = "HC0",
      scb_domain_length = domain_length
    ),
    control %||% list()
  )
  result <- .vcmoe_dev_intervals(
    fit = fit,
    data = data,
    level = level,
    strict = strict,
    control = inference_control
  )
  intervals <- .vcmoe_confband_filter(result$intervals, coefficient_set)
  if (identical(type, "pointwise")) {
    intervals$lower <- intervals$pointwise_lower
    intervals$upper <- intervals$pointwise_upper
  } else {
    intervals$lower <- intervals$simultaneous_lower
    intervals$upper <- intervals$simultaneous_upper
  }
  intervals$type <- type
  rownames(intervals) <- NULL
  out <- list(
    fit = fit,
    intervals = intervals,
    diagnostics = result$diagnostics,
    settings = list(
      family = fit$family,
      k = fit$k,
      engine_id = fit$engine_id %||% "local_grid_em",
      estimating_equation = if (joint_path) {
        "jasa_observed_local_likelihood_plugin"
      } else {
        "observed_local_likelihood"
      },
      covariance_scope = if (joint_path) {
        "local_asymptotic_plugin_excludes_finite_grid_cross_grid_coupling"
      } else {
        "local_sandwich"
      },
      covariance_target = "observed_local_marginal_likelihood",
      estimator_covariance_match = if (joint_path) {
        "asymptotic_plugin_not_exact_finite_grid_match"
      } else {
        "matched_when_local_fit_converged_and_penalty_negligible"
      },
      shared_path_uncertainty_accounted = FALSE,
      label_uncertainty_accounted = FALSE,
      coverage_theory = if (joint_path) {
        "local_likelihood_asymptotic_plugin_not_finite_grid_joint_path_theory"
      } else {
        "local_likelihood_asymptotic_with_bias_and_boundary_caveats"
      },
      level = level,
      type = type,
      coefficient_set = coefficient_set,
      covariance_adjustment = "HC0",
      simultaneous_method = "analytic_epanechnikov_path",
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id,
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      domain_length = domain_length
    )
  )
  class(out) <- "vcmoe_confband"
  out
}

#' @export
print.vcmoe_confband <- function(x, ...) {
  ok <- sum(x$intervals$status == "ok", na.rm = TRUE)
  total <- nrow(x$intervals)
  cat("VCMoE analytic-style confidence bands\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  components: ", x$settings$k, "\n", sep = "")
  cat("  engine: ", x$settings$engine_id, "\n", sep = "")
  cat("  type: ", x$settings$type, "\n", sep = "")
  cat("  level: ", x$settings$level, "\n", sep = "")
  cat("  coverage theory: ", x$settings$coverage_theory, "\n", sep = "")
  cat("  interval rows ok: ", ok, "/", total, "\n", sep = "")
  warning_rows <- sum(x$intervals$status == "warning", na.rm = TRUE)
  if (warning_rows) {
    cat("  interval rows warning: ", warning_rows, "/", total, "\n", sep = "")
  }
  invisible(x)
}
