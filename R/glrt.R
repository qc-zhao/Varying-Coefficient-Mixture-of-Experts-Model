.vcmoe_glrt_epanechnikov_kernel <- function(t) {
  0.75 * pmax(0, 1 - t^2)
}

.vcmoe_glrt_epanechnikov_convolution <- function(t) {
  f <- function(s) .vcmoe_glrt_epanechnikov_kernel(s) *
    .vcmoe_glrt_epanechnikov_kernel(t - s)
  lower <- max(-1, t - 1)
  upper <- min(1, t + 1)
  if (lower >= upper) {
    return(0)
  }
  stats::integrate(f, lower, upper, rel.tol = 1e-8)$value
}

.vcmoe_glrt_rk_delta_epanechnikov <- function(bandwidth, p_test = 1,
                                              C_const = 2, U_length = 1) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(p_test) || length(p_test) != 1L ||
      !is.finite(p_test) || p_test <= 0) {
    stop("`p_test` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(C_const) || length(C_const) != 1L ||
      !is.finite(C_const) || C_const <= 0) {
    stop("`C_const` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(U_length) || length(U_length) != 1L ||
      !is.finite(U_length) || U_length <= 0) {
    stop("`U_length` must be a positive finite number.", call. = FALSE)
  }
  int_K2 <- stats::integrate(
    function(t) .vcmoe_glrt_epanechnikov_kernel(t)^2,
    -1, 1, rel.tol = 1e-8
  )$value
  K0 <- .vcmoe_glrt_epanechnikov_kernel(0)
  g <- function(t) {
    .vcmoe_glrt_epanechnikov_kernel(t) -
      0.5 * vapply(t, .vcmoe_glrt_epanechnikov_convolution, numeric(1L))
  }
  denom <- stats::integrate(function(t) g(t)^2, -2, 2, rel.tol = 1e-8)$value
  num <- K0 - 0.5 * int_K2
  rK <- num / denom
  delta <- rK * p_test * C_const * U_length * num / bandwidth
  list(
    rK = as.numeric(rK),
    delta = as.numeric(delta),
    p_test = as.numeric(p_test),
    C_const = as.numeric(C_const),
    U_length = as.numeric(U_length),
    int_K2 = as.numeric(int_K2),
    K0 = as.numeric(K0)
  )
}

.vcmoe_glrt_control <- function(control) {
  defaults <- list(
    maxit = 200L,
    reltol = 1e-7,
    lambda_tol = 1e-7,
    strict = TRUE,
    min_converged_fraction = 0.80,
    min_component_proportion = 0.05,
    C_const = 2
  )
  utils::modifyList(defaults, control %||% list())
}

.vcmoe_glrt_normalize_calibration <- function(calibration) {
  calibration <- match.arg(
    calibration,
    c("none", "bootstrap", "analytic_epanechnikov", "both", "parametric_bootstrap")
  )
  if (identical(calibration, "parametric_bootstrap")) {
    return("bootstrap")
  }
  calibration
}

.vcmoe_glrt_public_set <- function(coefficient_set) {
  coefficient_set <- match.arg(coefficient_set, c("expert", "gating", "sigma", "theta"))
  coefficient_set
}

.vcmoe_glrt_select_term <- function(terms, term, label) {
  if (is.null(term)) {
    if (length(terms) == 1L) {
      return(terms[[1L]])
    }
    stop("`term` must be provided for coefficient-specific ", label,
         " tests. Available terms: ", paste(terms, collapse = ", "),
         call. = FALSE)
  }
  if (!is.character(term) || length(term) != 1L || !term %in% terms) {
    stop("Unknown `term` for ", label, " test. Available terms: ",
         paste(terms, collapse = ", "), call. = FALSE)
  }
  term
}

.vcmoe_glrt_select_component <- function(labels, component, label) {
  if (is.null(component)) {
    return(labels[[1L]])
  }
  if (is.character(component) && length(component) == 1L) {
    if (!component %in% labels) {
      stop("Unknown `component` for ", label, " test. Available components: ",
           paste(labels, collapse = ", "), call. = FALSE)
    }
    return(component)
  }
  component <- as.integer(component)
  if (length(component) != 1L || is.na(component) ||
      component < 1L || component > length(labels)) {
    stop("`component` must identify one ", label, " component.", call. = FALSE)
  }
  labels[[component]]
}

.vcmoe_glrt_select_gating_component <- function(fit, gating_labels, component) {
  if (is.null(component)) {
    return(gating_labels[[1L]])
  }
  component_labels <- dimnames(fit$coefficients$gating)[[2L]]
  baseline <- if (fit$k == 2L) component_labels[[2L]] else component_labels[[1L]]
  if (is.character(component) && length(component) == 1L) {
    if (component %in% gating_labels) {
      return(component)
    }
    if (component %in% component_labels) {
      if (identical(component, baseline) && fit$k > 2L) {
        stop("For `k > 2` gating tests, component 1 is the baseline; choose ",
             "a non-baseline component or a contrast label. Available contrasts: ",
             paste(gating_labels, collapse = ", "), call. = FALSE)
      }
      if (fit$k == 2L) {
        return(gating_labels[[1L]])
      }
      contrast <- paste(component, baseline, sep = "_vs_")
      if (contrast %in% gating_labels) {
        return(contrast)
      }
    }
    stop("Unknown `component` for gating test. Available components: ",
         paste(component_labels, collapse = ", "), "; available contrasts: ",
         paste(gating_labels, collapse = ", "), call. = FALSE)
  }
  component_id <- as.integer(component)
  if (length(component_id) != 1L || is.na(component_id) ||
      component_id < 1L || component_id > length(component_labels)) {
    stop("`component` must identify one gating component or contrast.",
         call. = FALSE)
  }
  if (fit$k == 2L) {
    return(gating_labels[[1L]])
  }
  if (component_id == 1L) {
    stop("For `k > 2` gating tests, component 1 is the baseline; choose ",
         "a non-baseline component or a contrast label. Available contrasts: ",
         paste(gating_labels, collapse = ", "), call. = FALSE)
  }
  contrast <- paste(component_labels[[component_id]], baseline, sep = "_vs_")
  if (!contrast %in% gating_labels) {
    stop("Selected gating contrast was not found in the fitted parameter table.",
         call. = FALSE)
  }
  contrast
}

.vcmoe_glrt_coefficient_selection <- function(fit, coefficient_set, component, term) {
  coefficient_set <- .vcmoe_glrt_public_set(coefficient_set)
  spec <- .vcmoe_dev_parameter_spec(fit, 1L)
  component_labels <- dimnames(fit$coefficients$expert)[[2L]]
  expert_terms <- dimnames(fit$coefficients$expert)[[3L]]
  gating_terms <- dimnames(fit$coefficients$gating)[[3L]]
  gating_labels <- unique(spec$parameter_table$component[
    spec$parameter_table$coefficient_set == "gating_contrast"
  ])

  if (identical(coefficient_set, "expert")) {
    return(list(
      public_set = "expert",
      internal_set = "expert",
      component = .vcmoe_glrt_select_component(component_labels, component, "expert"),
      term = .vcmoe_glrt_select_term(expert_terms, term, "expert")
    ))
  }
  if (identical(coefficient_set, "gating")) {
    selected_component <- .vcmoe_glrt_select_gating_component(
      fit,
      gating_labels,
      component
    )
    return(list(
      public_set = "gating",
      internal_set = "gating_contrast",
      component = selected_component,
      term = .vcmoe_glrt_select_term(gating_terms, term, "gating")
    ))
  }
  if (identical(coefficient_set, "sigma")) {
    if (!identical(fit$family, "gaussian")) {
      stop("`coefficient_set = \"sigma\"` is only available for Gaussian fits.",
           call. = FALSE)
    }
    return(list(
      public_set = "sigma",
      internal_set = "nuisance",
      component = .vcmoe_glrt_select_component(component_labels, component, "sigma"),
      term = .vcmoe_glrt_select_term("log_sigma", term, "sigma")
    ))
  }
  if (!identical(fit$family, "negative-binomial")) {
    stop("`coefficient_set = \"theta\"` is only available for Negative-Binomial fits.",
         call. = FALSE)
  }
  list(
    public_set = "theta",
    internal_set = "nuisance",
    component = .vcmoe_glrt_select_component(component_labels, component, "theta"),
    term = .vcmoe_glrt_select_term("log_theta", term, "theta")
  )
}

.vcmoe_glrt_constraint <- function(fit, test, coefficient_set, component, term) {
  test <- match.arg(test, c("coefficient", "constant_block", "constant_all"))
  table <- .vcmoe_dev_parameter_spec(fit, 1L)$parameter_table
  if (identical(test, "constant_all")) {
    shared <- table$parameter[table$block %in% c("intercept", "log_theta")]
    fixed <- table$parameter[table$block == "slope"]
    return(list(
      test = "constant_all",
      coefficient_set = "all",
      component = NA_character_,
      term = NA_character_,
      shared_parameters = shared,
      fixed_zero_parameters = fixed,
      p_test = length(shared),
      label = "all_coefficient_functions"
    ))
  }
  if (identical(test, "constant_block")) {
    if (!is.null(component) || !is.null(term)) {
      stop("`component` and `term` are not used for `test = \"constant_block\"`.",
           call. = FALSE)
    }
    coefficient_set <- .vcmoe_glrt_public_set(coefficient_set)
    if (!coefficient_set %in% c("expert", "gating")) {
      stop("`test = \"constant_block\"` supports `coefficient_set = \"expert\"` or `\"gating\"`.",
           call. = FALSE)
    }
    if (identical(coefficient_set, "gating")) {
      mask <- table$coefficient_set == "gating_contrast"
      internal_set <- "gating_contrast"
      label <- "gating_coefficient_functions"
    } else {
      mask <- table$coefficient_set %in% c("expert", "nuisance")
      internal_set <- "expert_nuisance"
      label <- "expert_and_dispersion_coefficient_functions"
    }
    shared <- table$parameter[mask & table$block %in% c("intercept", "log_theta")]
    fixed <- table$parameter[mask & table$block == "slope"]
    if (!length(shared)) {
      stop("Selected coefficient block has no intercept/log-theta parameter to constrain.",
           call. = FALSE)
    }
    return(list(
      test = "constant_block",
      coefficient_set = coefficient_set,
      internal_set = internal_set,
      component = NA_character_,
      term = NA_character_,
      shared_parameters = shared,
      fixed_zero_parameters = fixed,
      p_test = length(shared),
      label = label
    ))
  }

  selected <- .vcmoe_glrt_coefficient_selection(fit, coefficient_set, component, term)
  mask <- table$coefficient_set == selected$internal_set &
    table$component == selected$component &
    table$term == selected$term
  if (!any(mask)) {
    stop("Selected coefficient path was not found in the fitted parameter table.",
         call. = FALSE)
  }
  shared <- table$parameter[mask & table$block %in% c("intercept", "log_theta")]
  fixed <- table$parameter[mask & table$block == "slope"]
  if (!length(shared)) {
    stop("Selected coefficient path has no intercept/log-theta parameter to constrain.",
         call. = FALSE)
  }
  list(
    test = "coefficient",
    coefficient_set = selected$public_set,
    internal_set = selected$internal_set,
    component = selected$component,
    term = selected$term,
    shared_parameters = shared,
    fixed_zero_parameters = fixed,
    p_test = 1,
    label = paste(selected$public_set, selected$component, selected$term, sep = ":")
  )
}

.vcmoe_glrt_recompute_reduced_cache <- function(fit, constraint = NULL) {
  if (is.null(fit$fitted) || is.null(fit$fitted$u) || is.null(fit$u_grid)) {
    return(fit)
  }
  n_grid <- length(fit$u_grid)
  n <- length(fit$fitted$u)
  component_names <- dimnames(fit$coefficients$expert)[[2L]] %||%
    paste0("component", seq_len(fit$k))
  fit$posterior <- array(
    NA_real_,
    dim = c(n, fit$k, n_grid),
    dimnames = list(
      observation = NULL,
      component = component_names,
      u = signif(fit$u_grid, 8)
    )
  )
  for (grid_id in seq_along(fit$u_grid)) {
    spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
    par <- .vcmoe_dev_parameter_vector(fit, grid_id)
    params <- .vcmoe_dev_unpack_parameter(par, spec)
    log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
    log_norm <- .row_log_sum_exp(log_terms)
    posterior <- exp(log_terms - log_norm)
    posterior <- posterior / rowSums(posterior)
    fit$posterior[, , grid_id] <- posterior
    fit$diagnostics$loglik[[grid_id]] <- sum(spec$weights * log_norm)
    fit$diagnostics$posterior_entropy[[grid_id]] <- .weighted_entropy(
      posterior,
      spec$weights
    )
  }
  if (identical(fit$engine_id %||% "local_grid_em", "joint_path_em")) {
    assignment <- .nearest_grid_assignment(fit$fitted$u, fit$u_grid)
    path_posterior <- matrix(
      NA_real_, nrow = n, ncol = fit$k,
      dimnames = list(observation = NULL, component = component_names)
    )
    for (grid_id in seq_along(fit$u_grid)) {
      rows <- assignment == grid_id
      if (any(rows)) {
        path_posterior[rows, ] <- fit$posterior[rows, , grid_id, drop = FALSE][, , 1L]
      }
    }
    fit$joint_path_posterior <- path_posterior
    fit$diagnostics$joint_path_assignment <- data.frame(
      grid_id = seq_len(n_grid),
      u = fit$u_grid,
      u_original = fit$u_grid_original %||% rep(NA_real_, n_grid),
      n_assigned = tabulate(assignment, nbins = n_grid)
    )
    fit$diagnostics$source_engine <- "joint_path_em"
  }
  null_convergence <- fit$diagnostics$glrt_null_optim_convergence %||% NA_integer_
  fit$diagnostics$reduced_fit <- TRUE
  fit$diagnostics$reduced_constraint <- constraint$label %||% NA_character_
  fit$diagnostics$converged <- rep(
    is.finite(null_convergence) && null_convergence == 0L,
    n_grid
  )
  fit$diagnostics$iterations <- rep(NA_integer_, n_grid)
  fit$diagnostics$selected_start <- rep(NA_integer_, n_grid)
  fit$diagnostics$ambiguous <- rep(FALSE, n_grid)
  fit$diagnostics$alignment_margin <- rep(NA_real_, n_grid)
  for (field in c("permutations", "greedy_permutations", "global_permutations")) {
    if (!is.null(fit$diagnostics[[field]])) {
      fit$diagnostics[[field]][] <- NA_integer_
    }
  }
  for (field in c(
    "expert_optimizer_convergence", "expert_gradient_norm", "expert_coef_norm",
    "global_transition_cost"
  )) {
    if (!is.null(fit$diagnostics[[field]])) {
      fit$diagnostics[[field]][] <- NA
    }
  }
  fit$diagnostics$global_path_cost <- NA_real_
  fit$diagnostics$start_loglik <- vector("list", n_grid)
  fit$diagnostics$warnings <- unique(c(
    fit$diagnostics$warnings %||% character(0L),
    "Reduced fit diagnostics were recomputed from constrained coefficients; local EM iteration/start fields are not applicable."
  ))
  fit
}

.vcmoe_glrt_packing <- function(fit, constraint) {
  grid_pars <- lapply(seq_along(fit$u_grid), function(grid_id) {
    .vcmoe_dev_parameter_vector(fit, grid_id)
  })
  shared <- unique(constraint$shared_parameters)
  fixed <- unique(constraint$fixed_zero_parameters)
  shared_start <- vapply(shared, function(name) {
    mean(vapply(grid_pars, function(par) unname(par[[name]]), numeric(1L)), na.rm = TRUE)
  }, numeric(1L))
  names(shared_start) <- paste0("shared:", shared)

  free_entries <- data.frame()
  free_start <- numeric(0L)
  for (grid_id in seq_along(grid_pars)) {
    par <- grid_pars[[grid_id]]
    free_names <- setdiff(names(par), c(shared, fixed))
    if (length(free_names)) {
      free_entries <- rbind(
        free_entries,
        data.frame(
          grid_id = grid_id,
          parameter = free_names,
          stringsAsFactors = FALSE
        )
      )
      values <- par[free_names]
      names(values) <- paste0("grid", grid_id, ":", free_names)
      free_start <- c(free_start, values)
    }
  }

  list(
    fit = fit,
    constraint = constraint,
    grid_pars = grid_pars,
    shared = shared,
    fixed = fixed,
    free_entries = free_entries,
    start = c(shared_start, free_start)
  )
}

.vcmoe_glrt_unpack_grid_parameters <- function(par, packing) {
  n_shared <- length(packing$shared)
  shared_values <- if (n_shared) par[seq_len(n_shared)] else numeric(0L)
  names(shared_values) <- packing$shared
  free_values <- if (length(par) > n_shared) par[-seq_len(n_shared)] else numeric(0L)

  out <- packing$grid_pars
  for (grid_id in seq_along(out)) {
    grid_par <- out[[grid_id]]
    if (length(shared_values)) {
      grid_par[names(shared_values)] <- shared_values
    }
    if (length(packing$fixed)) {
      grid_par[packing$fixed] <- 0
    }
    out[[grid_id]] <- grid_par
  }
  if (length(free_values)) {
    for (row_id in seq_len(nrow(packing$free_entries))) {
      entry <- packing$free_entries[row_id, , drop = FALSE]
      out[[entry$grid_id]][entry$parameter] <- free_values[[row_id]]
    }
  }
  out
}

.vcmoe_glrt_apply_grid_parameters <- function(fit, grid_parameters) {
  null_fit <- fit
  for (grid_id in seq_along(grid_parameters)) {
    spec <- .vcmoe_dev_parameter_spec(null_fit, grid_id)
    params <- .vcmoe_dev_unpack_parameter(grid_parameters[[grid_id]], spec)
    null_fit$coefficients$expert[grid_id, , ] <- params$expert_coef[, seq_len(spec$q), drop = FALSE]
    null_fit$coefficients$expert_slope[grid_id, , ] <- params$expert_coef[, spec$q + seq_len(spec$q), drop = FALSE]
    null_fit$coefficients$gating[grid_id, , ] <- params$gating_coef[, seq_len(spec$p), drop = FALSE]
    null_fit$coefficients$gating_slope[grid_id, , ] <- params$gating_coef[, spec$p + seq_len(spec$p), drop = FALSE]
    if (!is.null(null_fit$coefficients$sigma)) {
      null_fit$coefficients$sigma[grid_id, ] <- params$sigma
    }
    if (!is.null(null_fit$coefficients$sigma_slope)) {
      null_fit$coefficients$sigma_slope[grid_id, ] <- params$sigma_slope
    }
    if (!is.null(null_fit$coefficients$theta)) {
      null_fit$coefficients$theta[grid_id, ] <- params$theta
    }
  }
  null_fit
}

.vcmoe_glrt_params_to_parameter_vector <- function(params, spec) {
  q <- spec$q
  p <- spec$p
  expert <- params$expert_coef[, seq_len(q), drop = FALSE]
  expert_slope <- params$expert_coef[, q + seq_len(q), drop = FALSE]
  gating <- params$gating_coef[, seq_len(p), drop = FALSE]
  gating_slope <- params$gating_coef[, p + seq_len(p), drop = FALSE]
  contrasts <- .vcmoe_dev_gating_contrasts(gating, gating_slope)
  par <- c(
    as.vector(t(expert)),
    as.vector(t(expert_slope)),
    as.vector(t(contrasts$intercept)),
    as.vector(t(contrasts$slope))
  )
  if (identical(spec$family, "gaussian")) {
    par <- c(
      par,
      log(pmax(params$sigma, spec$control$min_sigma)),
      params$sigma_slope %||% rep(0, spec$k)
    )
  } else if (identical(spec$family, "negative-binomial")) {
    par <- c(par, log(pmax(params$theta, spec$control$negbin_theta_min)))
  }
  names(par) <- spec$parameter_table$parameter
  par
}

.vcmoe_glrt_project_grid_parameters <- function(grid_parameters, constraint,
                                                projection_weights = NULL) {
  n_grid <- length(grid_parameters)
  if (is.null(projection_weights)) {
    projection_weights <- rep(1, n_grid)
    averaging_measure <- "uniform_grid"
  } else {
    projection_weights <- as.numeric(projection_weights)
    averaging_measure <- "nearest_grid_assignment_frequency"
  }
  if (length(projection_weights) != n_grid ||
      any(!is.finite(projection_weights)) ||
      any(projection_weights < 0) ||
      sum(projection_weights) <= 0) {
    stop(
      "`projection_weights` must be finite, nonnegative, match the grid, and have positive sum.",
      call. = FALSE
    )
  }
  shared <- unique(constraint$shared_parameters)
  fixed <- unique(constraint$fixed_zero_parameters)
  if (length(fixed)) {
    for (grid_id in seq_along(grid_parameters)) {
      grid_parameters[[grid_id]][fixed] <- 0
    }
  }
  projection <- data.frame(
    parameter = character(0L),
    projected_value = numeric(0L),
    averaging_measure = character(0L),
    weight_sum = numeric(0L),
    nonzero_grid_points = integer(0L),
    stringsAsFactors = FALSE
  )
  for (parameter in shared) {
    parameter_values <- vapply(
      grid_parameters,
      function(par) unname(par[[parameter]]),
      numeric(1L)
    )
    finite <- is.finite(parameter_values) & projection_weights > 0
    if (!any(finite)) {
      stop("No finite positively weighted values are available for projection.",
           call. = FALSE)
    }
    value <- stats::weighted.mean(
      parameter_values[finite],
      projection_weights[finite]
    )
    for (grid_id in seq_along(grid_parameters)) {
      grid_parameters[[grid_id]][parameter] <- value
    }
    projection <- rbind(
      projection,
      data.frame(
        parameter = parameter,
        projected_value = value,
        averaging_measure = averaging_measure,
        weight_sum = sum(projection_weights[finite]),
        nonzero_grid_points = sum(finite),
        stringsAsFactors = FALSE
      )
    )
  }
  list(parameters = grid_parameters, projection = projection)
}

.vcmoe_glrt_joint_path_initial_posterior <- function(fit) {
  n <- length(fit$fitted$u)
  if (!is.null(fit$joint_path_posterior) &&
      is.matrix(fit$joint_path_posterior) &&
      identical(dim(fit$joint_path_posterior), c(n, fit$k))) {
    posterior <- pmax(fit$joint_path_posterior, 1e-12)
    return(posterior / rowSums(posterior))
  }
  assignment <- .nearest_grid_assignment(fit$fitted$u, fit$u_grid)
  posterior <- matrix(NA_real_, nrow = n, ncol = fit$k)
  for (grid_id in seq_along(fit$u_grid)) {
    rows <- assignment == grid_id
    if (any(rows)) {
      posterior[rows, ] <- fit$posterior[rows, , grid_id, drop = FALSE][, , 1L]
    }
  }
  posterior <- pmax(posterior, 1e-12)
  posterior / rowSums(posterior)
}

.vcmoe_glrt_joint_path_null_fit <- function(fit, constraint, control = list()) {
  .vcmoe_dev_require_fit(fit)
  if (!identical(fit$engine_id %||% "local_grid_em", "joint_path_em")) {
    stop("Joint-path constrained null fitting requires a joint-path fit.",
         call. = FALSE)
  }
  if (is.null(fit$fitted) || is.null(fit$fitted$u)) {
    stop("Joint-path constrained null fitting requires retained fitted data.",
         call. = FALSE)
  }
  glrt_control <- .vcmoe_glrt_control(control)
  spec1 <- .vcmoe_dev_parameter_spec(fit, 1L)
  fit_control <- spec1$control
  tol <- min(fit_control$tol %||% 1e-6, glrt_control$reltol)
  y <- spec1$y
  trials <- spec1$trials
  z_design <- spec1$z_design
  x_design <- spec1$x_design
  offset <- spec1$offset
  u <- spec1$u
  assignment <- .nearest_grid_assignment(u, fit$u_grid)
  assignment_counts <- tabulate(assignment, nbins = length(fit$u_grid))
  path_posterior <- .vcmoe_glrt_joint_path_initial_posterior(fit)
  params_list <- lapply(seq_along(fit$u_grid), function(grid_id) {
    .vcmoe_dev_unpack_parameter(
      .vcmoe_dev_parameter_vector(fit, grid_id),
      .vcmoe_dev_parameter_spec(fit, grid_id)
    )
  })
  posterior_by_grid <- array(
    NA_real_, dim = c(length(u), fit$k, length(fit$u_grid))
  )
  loglik <- rep(NA_real_, length(fit$u_grid))
  trace <- data.frame(
    iteration = integer(0L),
    objective = numeric(0L),
    objective_delta = numeric(0L),
    posterior_delta = numeric(0L),
    parameter_delta = numeric(0L),
    converged = logical(0L)
  )
  projection_summary <- data.frame()
  objective <- -Inf
  converged <- FALSE

  for (iter in seq_len(as.integer(glrt_control$maxit))) {
    previous_posterior <- path_posterior
    previous_parameters <- lapply(seq_along(params_list), function(grid_id) {
      .vcmoe_glrt_params_to_parameter_vector(
        params_list[[grid_id]],
        .vcmoe_dev_parameter_spec(fit, grid_id)
      )
    })
    objective_old <- objective

    for (grid_id in seq_along(fit$u_grid)) {
      weights <- .local_weights(
        u, fit$u_grid[[grid_id]], fit$bandwidth, spec1$estimation_spec
      )
      du <- .local_du(
        u, fit$u_grid[[grid_id]], fit$bandwidth, spec1$estimation_spec
      )
      params_list[[grid_id]] <- .m_step(
        y, trials, fit$family, z_design, x_design, du, offset,
        weights, path_posterior, params_list[[grid_id]], fit_control
      )
    }

    grid_parameters <- lapply(seq_along(params_list), function(grid_id) {
      .vcmoe_glrt_params_to_parameter_vector(
        params_list[[grid_id]],
        .vcmoe_dev_parameter_spec(fit, grid_id)
      )
    })
    projected <- .vcmoe_glrt_project_grid_parameters(
      grid_parameters,
      constraint,
      projection_weights = assignment_counts
    )
    grid_parameters <- projected$parameters
    projection_summary <- projected$projection
    params_list <- lapply(seq_along(grid_parameters), function(grid_id) {
      .vcmoe_dev_unpack_parameter(
        grid_parameters[[grid_id]],
        .vcmoe_dev_parameter_spec(fit, grid_id)
      )
    })

    param_delta_values <- vapply(seq_along(grid_parameters), function(grid_id) {
      delta <- abs(grid_parameters[[grid_id]] - previous_parameters[[grid_id]])
      if (!any(is.finite(delta))) NA_real_ else max(delta[is.finite(delta)])
    }, numeric(1L))
    max_param_delta <- if (any(is.finite(param_delta_values))) {
      max(param_delta_values[is.finite(param_delta_values)])
    } else {
      NA_real_
    }

    for (grid_id in seq_along(fit$u_grid)) {
      weights <- .local_weights(
        u, fit$u_grid[[grid_id]], fit$bandwidth, spec1$estimation_spec
      )
      du <- .local_du(
        u, fit$u_grid[[grid_id]], fit$bandwidth, spec1$estimation_spec
      )
      posterior_by_grid[, , grid_id] <- .e_step(
        y, trials, fit$family, z_design, x_design, du, offset,
        params_list[[grid_id]], fit_control
      )
      loglik[[grid_id]] <- .local_loglik(
        y, trials, fit$family, z_design, x_design, du, offset,
        weights, params_list[[grid_id]], fit_control
      )
    }
    for (grid_id in seq_along(fit$u_grid)) {
      rows <- assignment == grid_id
      if (any(rows)) {
        path_posterior[rows, ] <- posterior_by_grid[rows, , grid_id, drop = FALSE][, , 1L]
      }
    }
    path_posterior <- pmax(path_posterior, 1e-12)
    path_posterior <- path_posterior / rowSums(path_posterior)
    objective <- .joint_path_sample_loglik(
      y,
      trials,
      fit$family,
      z_design,
      x_design,
      offset,
      u,
      fit$u_grid,
      fit$bandwidth,
      assignment,
      params_list,
      fit_control,
      spec1$estimation_spec
    )
    objective_delta <- if (is.finite(objective_old)) objective - objective_old else NA_real_
    posterior_delta <- max(abs(path_posterior - previous_posterior), na.rm = TRUE)
    did_converge <- is.finite(posterior_delta) &&
      posterior_delta <= tol &&
      (!is.finite(max_param_delta) || max_param_delta <= sqrt(tol))
    trace <- rbind(trace, data.frame(
      iteration = iter,
      objective = objective,
      objective_delta = objective_delta,
      posterior_delta = posterior_delta,
      parameter_delta = max_param_delta,
      converged = did_converge
    ))
    if (did_converge) {
      converged <- TRUE
      break
    }
  }

  final_parameters <- lapply(seq_along(params_list), function(grid_id) {
    .vcmoe_glrt_params_to_parameter_vector(
      params_list[[grid_id]],
      .vcmoe_dev_parameter_spec(fit, grid_id)
    )
  })
  final_parameters <- .vcmoe_glrt_project_grid_parameters(
    final_parameters,
    constraint,
    projection_weights = assignment_counts
  )$parameters
  null_fit <- .vcmoe_glrt_apply_grid_parameters(fit, final_parameters)
  null_fit$diagnostics$glrt_null_optim_convergence <- if (converged) 0L else 1L
  null_fit$diagnostics$glrt_null_optim_value <- if (is.finite(objective)) -objective else NA_real_
  null_fit$diagnostics$glrt_null_optim_error <- NA_character_
  null_fit$diagnostics$glrt_null_engine <- "joint_path_em_constrained_null"
  null_fit$diagnostics$glrt_statistic_criterion <- "nearest_grid_sample_loglik"
  null_fit$diagnostics$glrt_null_averaging_measure <-
    "nearest_grid_assignment_frequency"
  null_fit$diagnostics$glrt_null_constrained_mle <- FALSE
  null_fit$diagnostics$glrt_null_paper_match <- "nearest_grid_approximation"
  null_fit$diagnostics$glrt_constraint <- constraint$label
  null_fit <- .vcmoe_glrt_recompute_reduced_cache(null_fit, constraint)
  null_fit$joint_path_posterior <- path_posterior
  null_fit$diagnostics$joint_path_iterations <- if (nrow(trace)) {
    utils::tail(trace$iteration, 1L)
  } else {
    0L
  }
  null_fit$diagnostics$joint_path_converged <- converged
  null_fit$diagnostics$joint_path_trace <- trace
  null_fit$diagnostics$joint_path_selected_start <- NA_integer_
  null_fit$diagnostics$glrt_null_projection <- projection_summary
  null_fit$diagnostics$source_engine <- "joint_path_em"
  null_fit
}

.vcmoe_glrt_null_fit_bfgs <- function(fit, constraint, control = list()) {
  .vcmoe_dev_require_fit(fit)
  control <- .vcmoe_glrt_control(control)
  packing <- .vcmoe_glrt_packing(fit, constraint)
  objective <- function(par) {
    grid_parameters <- .vcmoe_glrt_unpack_grid_parameters(par, packing)
    value <- tryCatch(sum(vapply(seq_along(grid_parameters), function(grid_id) {
      .vcmoe_dev_local_loglik(
        fit,
        grid_id = grid_id,
        par = grid_parameters[[grid_id]],
        penalized = FALSE
      )
    }, numeric(1L))), error = function(e) NA_real_)
    if (!is.finite(value)) {
      return(1e100)
    }
    -value
  }
  opt <- tryCatch(
    stats::optim(
      packing$start,
      objective,
      method = "BFGS",
      control = list(maxit = as.integer(control$maxit), reltol = control$reltol)
    ),
    error = function(e) e
  )
  if (inherits(opt, "error")) {
    stop(
      "Local-grid constrained null optimization failed: ",
      conditionMessage(opt),
      call. = FALSE
    )
  }
  grid_parameters <- .vcmoe_glrt_unpack_grid_parameters(opt$par, packing)
  null_fit <- .vcmoe_glrt_apply_grid_parameters(fit, grid_parameters)
  null_fit$diagnostics$glrt_null_optim_convergence <- opt$convergence
  null_fit$diagnostics$glrt_null_optim_value <- opt$value
  null_fit$diagnostics$glrt_null_optim_error <- NA_character_
  null_fit$diagnostics$glrt_null_engine <- "local_grid_bfgs_constrained_null"
  null_fit$diagnostics$glrt_statistic_criterion <-
    "summed_kernel_weighted_local_pseudologlik"
  null_fit$diagnostics$glrt_null_averaging_measure <- "constrained_bfgs"
  null_fit$diagnostics$glrt_null_constrained_mle <- TRUE
  null_fit$diagnostics$glrt_null_paper_match <- "no"
  null_fit$diagnostics$glrt_constraint <- constraint$label
  .vcmoe_glrt_recompute_reduced_cache(null_fit, constraint)
}

.vcmoe_glrt_null_fit <- function(fit, constraint, control = list()) {
  if (identical(fit$engine_id %||% "local_grid_em", "joint_path_em")) {
    return(.vcmoe_glrt_joint_path_null_fit(fit, constraint, control = control))
  }
  .vcmoe_glrt_null_fit_bfgs(fit, constraint, control = control)
}

.vcmoe_glrt_global_pseudologlik <- function(fit) {
  if (identical(fit$engine_id %||% "local_grid_em", "joint_path_em")) {
    spec1 <- .vcmoe_dev_parameter_spec(fit, 1L)
    assignment <- .nearest_grid_assignment(spec1$u, fit$u_grid)
    params_list <- lapply(seq_along(fit$u_grid), function(grid_id) {
      .vcmoe_dev_unpack_parameter(
        .vcmoe_dev_parameter_vector(fit, grid_id),
        .vcmoe_dev_parameter_spec(fit, grid_id)
      )
    })
    return(.joint_path_sample_loglik(
      spec1$y,
      spec1$trials,
      fit$family,
      spec1$z_design,
      spec1$x_design,
      spec1$offset,
      spec1$u,
      fit$u_grid,
      fit$bandwidth,
      assignment,
      params_list,
      spec1$control,
      spec1$estimation_spec
    ))
  }
  .vcmoe_dev_global_pseudologlik(fit)
}

.vcmoe_glrt_statistic <- function(full_fit, null_fit, lambda_tol = 1e-7) {
  full_loglik <- .vcmoe_glrt_global_pseudologlik(full_fit)
  null_loglik <- .vcmoe_glrt_global_pseudologlik(null_fit)
  lambda <- full_loglik - null_loglik
  block_reason <- NA_character_
  if (!is.finite(lambda)) {
    block_reason <- "nonfinite_lambda"
  } else if (lambda < -lambda_tol) {
    block_reason <- "null_loglik_exceeds_full_loglik"
  }
  lambda_positive <- max(0, lambda)
  list(
    lambda = lambda,
    lambda_positive = lambda_positive,
    lrt_statistic = 2 * lambda_positive,
    full_loglik = full_loglik,
    null_loglik = null_loglik,
    block_reason = block_reason
  )
}

.vcmoe_glrt_fit_warnings <- function(fit, control) {
  diagnostics <- tryCatch(vcmoe_diagnostics(fit), error = function(e) NULL)
  reasons <- character(0L)
  if (mean(fit$diagnostics$converged %in% TRUE) < control$min_converged_fraction) {
    reasons <- c(reasons, "low_converged_fraction")
  }
  if (any(fit$diagnostics$ambiguous %||% FALSE, na.rm = TRUE)) {
    reasons <- c(reasons, "ambiguous_label")
  }
  if (!is.null(diagnostics) &&
      any(diagnostics$min_posterior_mean < control$min_component_proportion, na.rm = TRUE)) {
    reasons <- c(reasons, "component_collapse")
  }
  unique(reasons)
}

.vcmoe_glrt_require_analytic_fit <- function(fit) {
  .vcmoe_dev_require_fit(fit)
  metadata <- vcmoe_parameterization(fit)
  if (!identical(metadata$id, "a1_epanechnikov_scaled")) {
    stop("Analytic Epanechnikov GLRT requires `parameterization = \"a1_epanechnikov_scaled\"`.",
         call. = FALSE)
  }
  if (!identical(metadata$kernel$name, "epanechnikov") ||
      !identical(metadata$kernel$weight_normalization, "density_over_bandwidth") ||
      !identical(metadata$local_linear_basis$slope_storage, "scaled")) {
    stop("Analytic Epanechnikov GLRT requires Epanechnikov density weights and scaled local-linear slopes.",
         call. = FALSE)
  }
  method <- fit$u_scale %||% fit$u_scaling$method %||% "unit"
  if (!identical(method, "unit")) {
    stop("Analytic Epanechnikov GLRT requires `u_scale = \"unit\"`.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_glrt_analytic_calibration <- function(fit, constraint, observed, control) {
  .vcmoe_glrt_require_analytic_fit(fit)
  domain_length <- fit$u_scaling$domain_length %||%
    vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
  constants <- .vcmoe_glrt_rk_delta_epanechnikov(
    bandwidth = fit$bandwidth,
    p_test = constraint$p_test,
    C_const = control$C_const,
    U_length = domain_length
  )
  analytic_statistic <- constants$rK * observed$lambda_positive
  list(
    analytic_statistic = analytic_statistic,
    analytic_p_value = stats::pchisq(analytic_statistic, df = constants$delta, lower.tail = FALSE),
    rK = constants$rK,
    delta = constants$delta,
    p_test = constants$p_test,
    C_const = constants$C_const,
    U_length = constants$U_length,
    int_K2 = constants$int_K2,
    K0 = constants$K0
  )
}

.vcmoe_glrt_bootstrap <- function(fit, data, null_fit, constraint, B, seed,
                                  control, refit_control, verbose) {
  base_data <- .bootstrap_base_data(fit, data)
  u_info <- .bootstrap_u_info(fit, data, base_data, u = NULL)
  replicate_rows <- vector("list", B)
  successful_stats <- numeric(0L)

  for (replicate_id in seq_len(B)) {
    if (!is.null(seed)) {
      set.seed(as.integer(seed + replicate_id))
    }
    if (isTRUE(verbose)) {
      message("GLRT bootstrap replicate ", replicate_id, " / ", B)
    }
    start_time <- proc.time()[["elapsed"]]
    result <- tryCatch({
      boot_data <- .simulate_bootstrap_response(null_fit, base_data, u_info$values)
      fit_args <- list(
        formula = fit$formula,
        data = boot_data,
        u = u_info$refit,
        k = fit$k,
        family = fit$family,
        bandwidth = fit$bandwidth,
        u_grid = fit$u_grid,
        control = .bootstrap_refit_control(fit$control, refit_control, seed, replicate_id),
        label = .vcmoe_refit_label(fit),
        u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
        parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id %||% "a1_epanechnikov_scaled"
      )
      boot_fit <- suppressWarnings(.vcmoe_bootstrap_refit(fit, fit_args))
      boot_null <- .vcmoe_glrt_null_fit(boot_fit, constraint, control = control)
      statistic <- .vcmoe_glrt_statistic(
        boot_fit,
        boot_null,
        .vcmoe_glrt_control(control)$lambda_tol
      )
      statistic$engine_id <- boot_fit$engine_id %||% "local_grid_em"
      statistic$null_engine_id <- boot_null$diagnostics$glrt_null_engine %||%
        "local_grid_bfgs_constrained_null"
      statistic$null_convergence <- boot_null$diagnostics$glrt_null_optim_convergence %||%
        NA_integer_
      statistic$null_error <- boot_null$diagnostics$glrt_null_optim_error %||%
        NA_character_
      statistic
    }, error = function(e) e)
    runtime <- proc.time()[["elapsed"]] - start_time
    if (inherits(result, "error")) {
      replicate_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = "failed",
        engine_id = NA_character_,
        null_engine_id = NA_character_,
        null_convergence = NA_integer_,
        lambda = NA_real_,
        lrt_statistic = NA_real_,
        block_reason = NA_character_,
        runtime_seconds = runtime,
        error_message = conditionMessage(result),
        stringsAsFactors = FALSE
      )
    } else {
      null_ok <- is.finite(result$null_convergence) && result$null_convergence == 0L
      statistic_ok <- is.na(result$block_reason) && is.finite(result$lambda_positive)
      replicate_status <- if (null_ok && statistic_ok) "ok" else "blocked"
      if (identical(replicate_status, "ok")) {
        successful_stats <- c(successful_stats, result$lambda_positive)
      }
      replicate_rows[[replicate_id]] <- data.frame(
        replicate = replicate_id,
        status = replicate_status,
        engine_id = result$engine_id,
        null_engine_id = result$null_engine_id,
        null_convergence = result$null_convergence,
        lambda = result$lambda,
        lrt_statistic = result$lrt_statistic,
        block_reason = result$block_reason,
        runtime_seconds = runtime,
        error_message = if (null_ok) NA_character_ else result$null_error,
        stringsAsFactors = FALSE
      )
    }
  }
  list(
    replicate_summary = do.call(rbind, replicate_rows),
    successful_lambdas = successful_stats
  )
}

#' Fit a block-constant reduced VCMoE model
#'
#' Refits a VCMoE object under a block-constant coefficient constraint. A
#' local-grid reference uses the established constrained BFGS optimizer; a
#' joint-path reference uses a paper-inspired sample-weighted grid projection and
#' applies the selected constraint after every M-step. This projection is not a
#' generic constrained optimizer, and its diagnostic likelihood need not be
#' monotone.
#'
#' @param fit A `vcmoe` fit.
#' @param constrain Constraint to impose. `"gating_constant"` freezes the
#'   gating contrast functions in `u`; `"expert_constant"` freezes expert mean
#'   functions and family dispersion paths; `"all_constant"` freezes all
#'   fitted coefficient functions.
#' @param control Controls for constrained null fitting. With the default
#'   `strict = TRUE`, a nonconverged reduced fit is rejected.
#' @return A reduced object of class `vcmoe` with recomputed posterior and
#'   likelihood caches.
#' @export
vcmoe_fit_reduced <- function(fit,
                              constrain = c(
                                "gating_constant",
                                "expert_constant",
                                "all_constant"
                              ),
                              control = list()) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  constrain <- match.arg(constrain)
  spec <- switch(
    constrain,
    gating_constant = list(test = "constant_block", coefficient_set = "gating"),
    expert_constant = list(test = "constant_block", coefficient_set = "expert"),
    all_constant = list(test = "constant_all", coefficient_set = "expert")
  )
  constraint <- .vcmoe_glrt_constraint(
    fit,
    test = spec$test,
    coefficient_set = spec$coefficient_set,
    component = NULL,
    term = NULL
  )
  reduced <- .vcmoe_glrt_null_fit(fit, constraint, control = control)
  null_convergence <- reduced$diagnostics$glrt_null_optim_convergence %||% NA_integer_
  if (!is.finite(null_convergence) || null_convergence != 0L) {
    message <- paste0(
      "Constrained reduced-model fitting did not converge; inspect ",
      "`fit$diagnostics$joint_path_trace` or the null optimizer diagnostics."
    )
    if (isTRUE(.vcmoe_glrt_control(control)$strict)) {
      stop(message, call. = FALSE)
    }
    warning(message, call. = FALSE)
  }
  reduced$reduced_model <- list(
    constrain = constrain,
    constraint = constraint,
    source_class = class(fit),
    source_engine = fit$engine_id %||% "local_grid_em",
    null_engine = reduced$diagnostics$glrt_null_engine %||% NA_character_,
    null_averaging_measure = reduced$diagnostics$glrt_null_averaging_measure %||%
      NA_character_,
    constrained_mle = isTRUE(reduced$diagnostics$glrt_null_constrained_mle),
    paper_criterion_match = reduced$diagnostics$glrt_null_paper_match %||%
      NA_character_,
    null_convergence = reduced$diagnostics$glrt_null_optim_convergence,
    null_error = reduced$diagnostics$glrt_null_optim_error
  )
  reduced
}

#' Generalized likelihood-ratio test for VCMoE coefficient variation
#'
#' @param fit A `vcmoe` fit.
#' @param data Original data frame used to fit `fit`.
#' @param test Test type. `"coefficient"` tests one coefficient function;
#'   `"constant_block"` tests all expert or gating functions jointly;
#'   `"constant_all"` tests all fitted coefficient functions jointly.
#' @param coefficient_set Coefficient block for coefficient-specific or
#'   block-constant tests.
#' @param component Component label or index for coefficient-specific tests.
#' @param term Term name for coefficient-specific tests.
#' @param calibration Calibration method. The default `"none"` returns the
#'   statistic without attaching a reference distribution. `"analytic_epanechnikov"` uses the
#'   Epanechnikov modified chi-square calibration; `"bootstrap"` uses
#'   parametric bootstrap calibration; `"both"` reports both. The analytic
#'   calibration is retained as an explicitly requested approximation because
#'   the implemented statistic is not identical to the manuscript criterion.
#' @param B Number of bootstrap calibration replicates.
#' @param seed Optional random seed.
#' @param control Controls for constrained null optimization and diagnostics.
#' @param refit_control Controls overriding bootstrap full-model refits.
#' @param verbose Whether to message bootstrap progress.
#' @details Local-grid fits retain the 0.1.0 constrained BFGS null optimizer.
#'   Joint-path fits use a paper-inspired sample-weighted grid-projected null:
#'   after every M-step, each constrained coefficient path is replaced by its
#'   mean weighted by the number of observations assigned to each nearest grid
#'   point, and constrained local slopes are set to zero.
#'   Its statistic compares sample-level likelihood contributions evaluated at
#'   each observation's nearest grid point. The projected update is not a
#'   generic constrained optimizer and its diagnostic likelihood trace need not
#'   be monotone. Bootstrap calibration preserves both the full-fit engine and
#'   its matching null engine.
#' @return A `vcmoe_glrt` object.
#' @export
vcmoe_glrt <- function(fit, data,
                       test = c("coefficient", "constant_block", "constant_all"),
                       coefficient_set = c("expert", "gating", "sigma", "theta"),
                       component = NULL, term = NULL,
                       calibration = c("none", "bootstrap", "analytic_epanechnikov", "both", "parametric_bootstrap"),
                       B = 200L, seed = NULL, control = list(),
                       refit_control = list(), verbose = FALSE) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  test <- match.arg(test)
  calibration <- .vcmoe_glrt_normalize_calibration(calibration)
  glrt_control <- .vcmoe_glrt_control(control)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 1L || abs(B - round(B)) > 1e-8) {
    stop("`B` must be a whole number of at least 1.", call. = FALSE)
  }
  B <- as.integer(B)
  if (calibration %in% c("analytic_epanechnikov", "both")) {
    .vcmoe_glrt_require_analytic_fit(fit)
  }

  constraint <- .vcmoe_glrt_constraint(fit, test, coefficient_set[[1L]], component, term)
  null_fit <- .vcmoe_glrt_null_fit(fit, constraint, control = control)
  observed <- .vcmoe_glrt_statistic(fit, null_fit, glrt_control$lambda_tol)

  warnings <- .vcmoe_glrt_fit_warnings(fit, glrt_control)
  if (any(c(.vcmoe_default_control(fit$control %||% list())$ridge,
            .vcmoe_default_control(fit$control %||% list())$binomial_ridge,
            .vcmoe_default_control(fit$control %||% list())$negbin_ridge,
            .vcmoe_default_control(fit$control %||% list())$negbin_theta_ridge) > 0)) {
    warnings <- c(warnings, "ridge_used_for_fitting_excluded_from_glrt_statistic")
  }
  null_convergence <- null_fit$diagnostics$glrt_null_optim_convergence
  if (!is.finite(null_convergence) || null_convergence != 0L) {
    warnings <- c(warnings, "null_optimizer_not_fully_converged")
  }

  block_reasons <- character(0L)
  if (!is.na(observed$block_reason)) {
    block_reasons <- c(block_reasons, observed$block_reason)
  }
  if (!is.finite(null_convergence) || null_convergence != 0L) {
    block_reasons <- c(block_reasons, "null_optimizer_not_fully_converged")
  }
  null_error <- null_fit$diagnostics$glrt_null_optim_error %||% NA_character_
  if (!is.na(null_error) && nzchar(null_error)) {
    block_reasons <- c(block_reasons, "null_optimizer_error")
  }
  if (isTRUE(glrt_control$strict)) {
    block_reasons <- c(block_reasons, setdiff(warnings, "ridge_used_for_fitting_excluded_from_glrt_statistic"))
  }
  substantive_warnings <- setdiff(
    unique(warnings),
    "ridge_used_for_fitting_excluded_from_glrt_statistic"
  )
  status <- if (length(block_reasons)) {
    "blocked"
  } else if (length(substantive_warnings)) {
    "warning"
  } else {
    "ok"
  }
  block_reason <- if (length(block_reasons)) paste(unique(block_reasons), collapse = ";") else NA_character_
  inference_available <- status %in% c("ok", "warning")

  analytic <- list(
    analytic_statistic = NA_real_,
    analytic_p_value = NA_real_,
    rK = NA_real_,
    delta = NA_real_,
    p_test = constraint$p_test,
    C_const = glrt_control$C_const,
    U_length = fit$u_scaling$domain_length %||% 1,
    int_K2 = NA_real_,
    K0 = NA_real_
  )
  if (calibration %in% c("analytic_epanechnikov", "both") && inference_available) {
    analytic <- .vcmoe_glrt_analytic_calibration(fit, constraint, observed, glrt_control)
  }

  bootstrap_p_value <- NA_real_
  replicate_summary <- data.frame()
  if (calibration %in% c("bootstrap", "both") && inference_available) {
    boot <- .vcmoe_glrt_bootstrap(
      fit, data, null_fit, constraint, B, seed,
      control = control,
      refit_control = refit_control,
      verbose = verbose
    )
    replicate_summary <- boot$replicate_summary
    if (length(boot$successful_lambdas)) {
      bootstrap_p_value <- (1 + sum(boot$successful_lambdas >= observed$lambda_positive)) /
        (1 + length(boot$successful_lambdas))
    }
  }

  p_value <- if (identical(calibration, "analytic_epanechnikov")) {
    analytic$analytic_p_value
  } else if (identical(calibration, "bootstrap")) {
    bootstrap_p_value
  } else if (identical(calibration, "both")) {
    analytic$analytic_p_value
  } else {
    NA_real_
  }

  fit_defaults <- .vcmoe_default_control(fit$control %||% list())
  joint_path <- identical(fit$engine_id %||% "local_grid_em", "joint_path_em")
  statistic_criterion <- if (joint_path) {
    "nearest_grid_sample_loglik"
  } else {
    "summed_kernel_weighted_local_pseudologlik"
  }
  paper_criterion_match <- if (joint_path) "nearest_grid_approximation" else "no"
  null_averaging_measure <- null_fit$diagnostics$glrt_null_averaging_measure %||%
    if (joint_path) "nearest_grid_assignment_frequency" else "constrained_bfgs"
  constrained_mle <- isTRUE(null_fit$diagnostics$glrt_null_constrained_mle)
  theory_status <- if (joint_path) {
    "nearest_grid_approximation_not_exact_manuscript_theory"
  } else {
    "analytic_calibration_not_established_for_implementation_criterion"
  }
  calibration_status <- if (identical(status, "blocked")) {
    "blocked"
  } else {
    switch(
      calibration,
      none = "not_requested",
      bootstrap = "empirical_parametric_bootstrap",
      analytic_epanechnikov = "approximate_not_theorem_matched",
      both = "analytic_approximation_plus_empirical_bootstrap"
    )
  }
  out <- list(
    fit = fit,
    null_fit = null_fit,
    status = status,
    block_reason = block_reason,
    warnings = unique(warnings),
    statistic = analytic$analytic_statistic,
    lambda = observed$lambda,
    lrt_statistic = observed$lrt_statistic,
    analytic_statistic = analytic$analytic_statistic,
    full_loglik = observed$full_loglik,
    null_loglik = observed$null_loglik,
    p_value = p_value,
    analytic_p_value = analytic$analytic_p_value,
    bootstrap_p_value = bootstrap_p_value,
    rK = analytic$rK,
    delta = analytic$delta,
    p_test = analytic$p_test,
    C_const = analytic$C_const,
    U_length = analytic$U_length,
    int_K2 = analytic$int_K2,
    K0 = analytic$K0,
    null_convergence = null_convergence,
    replicate_summary = replicate_summary,
    settings = list(
      test = test,
      coefficient_set = constraint$coefficient_set,
      component = constraint$component,
      term = constraint$term,
      constraint = constraint$label,
      calibration = calibration,
      B = B,
      seed = seed,
      family = fit$family,
      k = fit$k,
      bandwidth = fit$bandwidth,
      parameterization = fit$parameterization_id %||% vcmoe_parameterization(fit)$id,
      u_scale = fit$u_scale %||% fit$u_scaling$method %||% "unit",
      ridge = fit_defaults$ridge,
      binomial_ridge = fit_defaults$binomial_ridge,
      negbin_ridge = fit_defaults$negbin_ridge,
      negbin_theta_ridge = fit_defaults$negbin_theta_ridge,
      negbin_theta_target = fit_defaults$negbin_theta_target,
      engine_id = fit$engine_id %||% "local_grid_em",
      null_engine_id = null_fit$diagnostics$glrt_null_engine %||%
        "local_grid_bfgs_constrained_null",
      statistic_criterion = statistic_criterion,
      paper_criterion_match = paper_criterion_match,
      null_averaging_measure = null_averaging_measure,
      constrained_mle = constrained_mle,
      calibration_status = calibration_status,
      theory_status = theory_status,
      statistic_note = if (joint_path) {
        paste0(
          "Joint-path calibration uses sample-level nearest-grid log-likelihood ",
          "contributions and a paper-inspired sample-weighted grid-projected null; ",
          "ridge penalties are excluded."
        )
      } else {
        "Analytic calibration uses rK * (full_loglik - null_loglik); ridge penalties are excluded."
      }
    )
  )
  class(out) <- "vcmoe_glrt"
  out
}

#' @export
print.vcmoe_glrt <- function(x, ...) {
  cat("VCMoE GLRT\n")
  cat("  family: ", x$settings$family, "\n", sep = "")
  cat("  test: ", x$settings$test, "\n", sep = "")
  if (identical(x$settings$test, "coefficient")) {
    cat("  coefficient: ", x$settings$coefficient_set, " / ",
        x$settings$component, " / ", x$settings$term, "\n", sep = "")
  }
  cat("  status: ", x$status, "\n", sep = "")
  if (identical(x$status, "blocked")) {
    cat("  block reason: ", x$block_reason, "\n", sep = "")
  } else if (identical(x$status, "warning") && length(x$warnings)) {
    cat("  warnings: ", paste(x$warnings, collapse = ";"), "\n", sep = "")
  }
  cat("  lambda: ", signif(x$lambda, 5), "\n", sep = "")
  cat("  analytic statistic: ", signif(x$analytic_statistic, 5), "\n", sep = "")
  cat("  calibration: ", x$settings$calibration, "\n", sep = "")
  cat("  calibration status: ", x$settings$calibration_status, "\n", sep = "")
  cat("  theory status: ", x$settings$theory_status, "\n", sep = "")
  if (is.finite(x$analytic_p_value)) {
    cat("  analytic p-value: ", signif(x$analytic_p_value, 4), "\n", sep = "")
  }
  if (is.finite(x$bootstrap_p_value)) {
    cat("  bootstrap p-value: ", signif(x$bootstrap_p_value, 4), "\n", sep = "")
  }
  invisible(x)
}
