.row_nonzero_count <- function(matrix_value) {
  apply(matrix_value, 1L, function(row) {
    if (all(is.na(row))) {
      return(NA_integer_)
    }
    sum(is.finite(row) & row != 0)
  })
}

.row_na_count <- function(matrix_value) {
  apply(matrix_value, 1L, function(row) {
    if (all(is.na(row))) {
      return(NA_integer_)
    }
    sum(is.na(row))
  })
}

.row_max_finite <- function(matrix_value) {
  apply(matrix_value, 1L, function(row) {
    if (!any(is.finite(row))) {
      return(NA_real_)
    }
    max(row[is.finite(row)])
  })
}

.row_min_finite <- function(matrix_value) {
  apply(matrix_value, 1L, function(row) {
    if (!any(is.finite(row))) {
      return(NA_real_)
    }
    min(row[is.finite(row)])
  })
}

.diagnostic_effective_n <- function(object) {
  weights <- .diagnostic_local_weights(object)
  if (is.null(weights)) {
    return(rep(NA_real_, length(object$u_grid)))
  }
  apply(weights, 2L, function(w) {
    sum(w)^2 / sum(w^2)
  })
}

.diagnostic_local_weights <- function(object) {
  if (is.null(object$fitted) || is.null(object$fitted$u) ||
      is.null(object$bandwidth) || !is.finite(object$bandwidth)) {
    return(NULL)
  }
  estimation_spec <- .vcmoe_estimation_spec(
    object$parameterization_id %||% vcmoe_parameterization(object)$id %||% "a1_epanechnikov_scaled"
  )
  weights <- matrix(NA_real_, nrow = length(object$fitted$u), ncol = length(object$u_grid))
  for (grid_id in seq_along(object$u_grid)) {
    w <- tryCatch(
      .local_weights(object$fitted$u, object$u_grid[[grid_id]], object$bandwidth, estimation_spec),
      error = function(e) NULL
    )
    if (is.null(w)) {
      return(NULL)
    }
    weights[, grid_id] <- w
  }
  weights
}

.diagnostic_component_means <- function(object) {
  weights <- .diagnostic_local_weights(object)
  posterior_mean <- matrix(NA_real_, nrow = length(object$u_grid), ncol = object$k)
  for (grid_id in seq_along(object$u_grid)) {
    posterior_grid <- object$posterior[, , grid_id, drop = FALSE][, , 1L]
    if (is.null(weights)) {
      posterior_mean[grid_id, ] <- colMeans(posterior_grid)
    } else {
      w <- weights[, grid_id]
      posterior_mean[grid_id, ] <- colSums(posterior_grid * w) / sum(w)
    }
  }
  colnames(posterior_mean) <- paste0("posterior_", dimnames(object$posterior)[[2L]])
  posterior_mean
}

.diagnostic_component_effective_n <- function(object) {
  weights <- .diagnostic_local_weights(object)
  component_ess <- matrix(NA_real_, nrow = length(object$u_grid), ncol = object$k)
  if (is.null(weights) || is.null(object$posterior)) {
    colnames(component_ess) <- paste0("component_effective_n_", dimnames(object$posterior)[[2L]] %||% seq_len(object$k))
    return(component_ess)
  }
  for (grid_id in seq_along(object$u_grid)) {
    posterior_grid <- object$posterior[, , grid_id, drop = FALSE][, , 1L]
    w <- weights[, grid_id]
    for (component in seq_len(object$k)) {
      wc <- w * posterior_grid[, component]
      denom <- sum(wc^2)
      component_ess[grid_id, component] <- if (is.finite(denom) && denom > 0) {
        sum(wc)^2 / denom
      } else {
        NA_real_
      }
    }
  }
  colnames(component_ess) <- paste0("component_effective_n_", dimnames(object$posterior)[[2L]])
  component_ess
}

.diagnostic_theta_summary <- function(object) {
  n_grid <- length(object$u_grid)
  if (!identical(object$family, "negative-binomial") || is.null(object$coefficients$theta)) {
    return(data.frame(
      theta_min = rep(NA_real_, n_grid),
      theta_max = rep(NA_real_, n_grid),
      theta_boundary_count = rep(NA_integer_, n_grid),
      theta_boundary = rep(NA, n_grid)
    ))
  }
  control <- .vcmoe_default_control(object$control %||% list())
  theta <- object$coefficients$theta
  boundary <- theta <= control$negbin_theta_min * 1.01 |
    theta >= control$negbin_theta_max / 1.01
  data.frame(
    theta_min = apply(theta, 1L, min, na.rm = TRUE),
    theta_max = apply(theta, 1L, max, na.rm = TRUE),
    theta_boundary_count = rowSums(boundary, na.rm = TRUE),
    theta_boundary = rowSums(boundary, na.rm = TRUE) > 0
  )
}

.diagnostic_negbin_eta_summary <- function(object) {
  n_grid <- length(object$u_grid)
  empty <- data.frame(
    negbin_eta_min_observed = rep(NA_real_, n_grid),
    negbin_eta_max_observed = rep(NA_real_, n_grid),
    negbin_eta_clipping_count = rep(NA_integer_, n_grid),
    negbin_eta_near_clipping_count = rep(NA_integer_, n_grid)
  )
  if (!identical(object$family, "negative-binomial") || is.null(object$fitted)) {
    return(empty)
  }
  control <- .vcmoe_default_control(object$control %||% list())
  estimation_spec <- .vcmoe_estimation_spec(
    object$parameterization_id %||% vcmoe_parameterization(object)$id %||% "a1_epanechnikov_scaled"
  )
  z_design <- object$fitted$z_design
  offset <- object$fitted$expert_offset %||% rep(0, nrow(z_design))
  out <- empty
  for (grid_id in seq_along(object$u_grid)) {
    du <- .local_du(object$fitted$u, object$u_grid[[grid_id]], object$bandwidth, estimation_spec)
    eta_values <- numeric(0L)
    for (component in seq_len(object$k)) {
      eta <- offset +
        as.numeric(z_design %*% object$coefficients$expert[grid_id, component, ]) +
        du * as.numeric(z_design %*% object$coefficients$expert_slope[grid_id, component, ])
      eta_values <- c(eta_values, eta)
    }
    eta_values <- eta_values[is.finite(eta_values)]
    if (length(eta_values)) {
      out$negbin_eta_min_observed[[grid_id]] <- min(eta_values)
      out$negbin_eta_max_observed[[grid_id]] <- max(eta_values)
      out$negbin_eta_clipping_count[[grid_id]] <- sum(
        eta_values <= control$negbin_eta_min | eta_values >= control$negbin_eta_max
      )
      out$negbin_eta_near_clipping_count[[grid_id]] <- sum(
        eta_values <= control$negbin_eta_min + control$negbin_eta_margin |
          eta_values >= control$negbin_eta_max - control$negbin_eta_margin
      )
    }
  }
  out
}

#' Summarize VCMoE fit diagnostics
#'
#' @param object A `vcmoe` object.
#' @return A data frame with one row per coefficient grid point.
#' @export
vcmoe_diagnostics <- function(object) {
  if (!inherits(object, "vcmoe")) {
    stop("`object` must be a VCMoE fit.", call. = FALSE)
  }
  n_grid <- length(object$u_grid)
  diagnostics <- object$diagnostics
  ambiguous <- diagnostics$ambiguous %||% rep(FALSE, n_grid)
  expert_convergence <- diagnostics$expert_optimizer_convergence
  expert_gradient_norm <- diagnostics$expert_gradient_norm
  expert_coef_norm <- diagnostics$expert_coef_norm
  if (is.null(expert_convergence)) {
    expert_convergence <- matrix(NA_integer_, nrow = n_grid, ncol = object$k)
  }
  if (is.null(expert_gradient_norm)) {
    expert_gradient_norm <- matrix(NA_real_, nrow = n_grid, ncol = object$k)
  }
  if (is.null(expert_coef_norm)) {
    expert_coef_norm <- matrix(NA_real_, nrow = n_grid, ncol = object$k)
  }

  posterior_mean <- .diagnostic_component_means(object)
  component_effective_n <- .diagnostic_component_effective_n(object)
  theta_summary <- .diagnostic_theta_summary(object)
  negbin_eta_summary <- .diagnostic_negbin_eta_summary(object)
  out <- data.frame(
    grid_id = seq_len(n_grid),
    u = object$u_grid,
    engine_id = object$engine_id %||% "local_grid_em",
    family = object$family,
    k = object$k,
    bandwidth = object$bandwidth,
    converged = diagnostics$converged,
    iterations = diagnostics$iterations,
    loglik = diagnostics$loglik,
    selected_start = diagnostics$selected_start,
    posterior_entropy = diagnostics$posterior_entropy,
    ambiguous = ambiguous,
    alignment_margin = diagnostics$alignment_margin,
    effective_n = .diagnostic_effective_n(object),
    min_component_effective_n = .row_min_finite(component_effective_n),
    min_posterior_mean = apply(posterior_mean, 1L, min),
    max_posterior_mean = apply(posterior_mean, 1L, max),
    expert_optimizer_nonzero_count = .row_nonzero_count(expert_convergence),
    expert_optimizer_na_count = .row_na_count(expert_convergence),
    expert_gradient_norm_max = .row_max_finite(expert_gradient_norm),
    expert_coef_norm_max = .row_max_finite(expert_coef_norm),
    stringsAsFactors = FALSE
  )
  cbind(
    out,
    theta_summary,
    negbin_eta_summary,
    as.data.frame(posterior_mean, stringsAsFactors = FALSE),
    as.data.frame(component_effective_n, stringsAsFactors = FALSE)
  )
}

.diagnostics_plot_data <- function(diagnostics) {
  component_columns <- grep("^posterior_component", names(diagnostics), value = TRUE)
  base <- data.frame(
    u = rep(diagnostics$u, times = 5L),
    metric = rep(
      c("converged", "posterior entropy", "effective local sample size",
        "minimum component proportion", "ambiguity warning"),
      each = nrow(diagnostics)
    ),
    series = "fit",
    value = c(
      as.numeric(diagnostics$converged),
      diagnostics$posterior_entropy,
      diagnostics$effective_n,
      diagnostics$min_posterior_mean,
      as.numeric(diagnostics$ambiguous)
    ),
    stringsAsFactors = FALSE
  )
  component_df <- data.frame()
  if (length(component_columns)) {
    component_df <- do.call(rbind, lapply(component_columns, function(column) {
      data.frame(
        u = diagnostics$u,
        metric = "component proportion",
        series = sub("^posterior_", "", column),
        value = diagnostics[[column]],
        stringsAsFactors = FALSE
      )
    }))
  }
  out <- rbind(base, component_df)
  out$metric <- factor(
    out$metric,
    levels = c(
      "converged",
      "posterior entropy",
      "component proportion",
      "minimum component proportion",
      "effective local sample size",
      "ambiguity warning"
    )
  )
  out
}

#' Plot VCMoE fit diagnostics
#'
#' @param object A `vcmoe` object.
#' @return A `ggplot` object.
#' @export
plot_diagnostics <- function(object) {
  diagnostics <- vcmoe_diagnostics(object)
  plot_data <- .diagnostics_plot_data(diagnostics)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = u, y = value, color = series)) +
    ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.6, na.rm = TRUE) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::labs(x = "u", y = "diagnostic value", color = "series") +
    ggplot2::theme_minimal(base_size = 12)
}
