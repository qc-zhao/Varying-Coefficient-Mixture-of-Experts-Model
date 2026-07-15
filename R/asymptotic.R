.vcmoe_dev_simultaneous_method <- function(method) {
  method <- method %||% "multiplier_grid"
  method <- match.arg(
    method,
    c("multiplier_grid", "multiplier", "analytic_epanechnikov_path",
      "analytic", "analytic_path")
  )
  if (identical(method, "multiplier")) {
    return("multiplier_grid")
  }
  if (method %in% c("analytic", "analytic_path")) {
    return("analytic_epanechnikov_path")
  }
  method
}

.vcmoe_dev_covariance_adjustment <- function(method) {
  method <- method %||% "HC0"
  if (!is.character(method) || length(method) != 1L || is.na(method) || !nzchar(method)) {
    stop("`covariance_adjustment` must be \"HC0\". The alias \"none\" is accepted as HC0.",
         call. = FALSE)
  }
  key <- tolower(method)
  if (key %in% c("hc0", "none")) {
    return("HC0")
  }
  stop("Analytic confidence bands use HC0 only. Set `covariance_adjustment = \"HC0\"`.",
       call. = FALSE)
}

.vcmoe_dev_control <- function(control) {
  defaults <- list(
    level = 0.90,
    finite_diff_eps = 1e-5,
    hessian_eps = 1e-4,
    max_hessian_dim = 45L,
    max_condition = 1e8,
    min_eigenvalue = 1e-8,
    min_converged_fraction = 0.80,
    min_component_proportion = 0.05,
    min_effective_n = 10,
    multiplier_B = 199L,
    simultaneous_method = "multiplier_grid",
    covariance_adjustment = "HC0",
    scb_domain_length = NULL,
    scb_kernel_constant = NULL,
    seed = NULL,
    min_alignment_margin = .vcmoe_default_control(list())$align_margin_tol,
    negbin_eta_margin = 0.25,
    strict = TRUE
  )
  out <- utils::modifyList(defaults, control %||% list())
  out$simultaneous_method <- .vcmoe_dev_simultaneous_method(out$simultaneous_method)
  out$covariance_adjustment <- .vcmoe_dev_covariance_adjustment(out$covariance_adjustment)
  out
}

.vcmoe_dev_require_fit <- function(fit) {
  if (!inherits(fit, "vcmoe")) {
    stop("`fit` must be a VCMoE fit.", call. = FALSE)
  }
  if (is.null(fit$fitted)) {
    stop("Development asymptotic helpers require `fit$fitted`; refit with `control = list(keep_data = TRUE)`.", call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_dev_family_response <- function(fit) {
  response <- fit$fitted$response
  if (identical(fit$family, "gaussian")) {
    return(list(y = response$y, trials = NULL))
  }
  if (identical(fit$family, "binomial")) {
    return(list(y = response$y, trials = response$trials))
  }
  if (identical(fit$family, "negative-binomial")) {
    return(list(y = response$y, trials = NULL))
  }
  stop("Unsupported family.", call. = FALSE)
}

.vcmoe_dev_gating_contrasts <- function(gating, gating_slope) {
  k <- nrow(gating)
  if (k == 2L) {
    return(list(
      intercept = matrix(gating[1L, ] - gating[2L, ], nrow = 1L),
      slope = matrix(gating_slope[1L, ] - gating_slope[2L, ], nrow = 1L),
      labels = "component1_vs_component2"
    ))
  }
  labels <- paste0("component", seq.int(2L, k), "_vs_component1")
  list(
    intercept = gating[-1L, , drop = FALSE] -
      matrix(gating[1L, ], nrow = k - 1L, ncol = ncol(gating), byrow = TRUE),
    slope = gating_slope[-1L, , drop = FALSE] -
      matrix(gating_slope[1L, ], nrow = k - 1L, ncol = ncol(gating_slope), byrow = TRUE),
    labels = labels
  )
}

.vcmoe_dev_gating_from_contrasts <- function(intercept, slope, k) {
  if (k == 2L) {
    return(list(
      intercept = rbind(intercept[1L, ] / 2, -intercept[1L, ] / 2),
      slope = rbind(slope[1L, ] / 2, -slope[1L, ] / 2)
    ))
  }
  beta_intercept <- rbind(rep(0, ncol(intercept)), intercept)
  beta_slope <- rbind(rep(0, ncol(slope)), slope)
  list(
    intercept = .center_logits(beta_intercept),
    slope = .center_logits(beta_slope)
  )
}

.vcmoe_dev_parameter_spec <- function(fit, grid_id) {
  .vcmoe_dev_require_fit(fit)
  if (grid_id < 1L || grid_id > length(fit$u_grid)) {
    stop("`grid_id` is outside the fitted `u_grid`.", call. = FALSE)
  }
  k <- fit$k
  q <- dim(fit$coefficients$expert)[3L]
  p <- dim(fit$coefficients$gating)[3L]
  expert_terms <- dimnames(fit$coefficients$expert)[[3L]]
  gating_terms <- dimnames(fit$coefficients$gating)[[3L]]
  component_labels <- dimnames(fit$coefficients$expert)[[2L]]
  gating_labels <- if (k == 2L) "component1_vs_component2" else paste0("component", seq.int(2L, k), "_vs_component1")

  parameter_table <- data.frame()
  add_rows <- function(block, coefficient_set, components, terms) {
    grid <- expand.grid(
      term = terms,
      component = components,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    data.frame(
      block = block,
      coefficient_set = coefficient_set,
      component = grid$component,
      term = grid$term,
      stringsAsFactors = FALSE
    )
  }
  parameter_table <- rbind(
    parameter_table,
    add_rows("intercept", "expert", component_labels, expert_terms),
    add_rows("slope", "expert", component_labels, expert_terms),
    add_rows("intercept", "gating_contrast", gating_labels, gating_terms),
    add_rows("slope", "gating_contrast", gating_labels, gating_terms)
  )
  if (identical(fit$family, "gaussian")) {
    parameter_table <- rbind(
      parameter_table,
      add_rows("intercept", "nuisance", component_labels, "log_sigma"),
      add_rows("slope", "nuisance", component_labels, "log_sigma")
    )
  } else if (identical(fit$family, "negative-binomial")) {
    parameter_table <- rbind(
      parameter_table,
      add_rows("log_theta", "nuisance", component_labels, "log_theta")
    )
  }
  parameter_table$parameter <- paste(
    parameter_table$coefficient_set,
    parameter_table$component,
    parameter_table$term,
    parameter_table$block,
    sep = ":"
  )
  parameter_table$index <- seq_len(nrow(parameter_table))

  defaults <- .vcmoe_default_control(fit$control %||% list())
  response <- .vcmoe_dev_family_response(fit)
  estimation_spec <- .vcmoe_estimation_spec(
    fit$parameterization_id %||% vcmoe_parameterization(fit)$id %||% "a1_epanechnikov_scaled"
  )
  weights <- .local_weights(fit$fitted$u, fit$u_grid[[grid_id]], fit$bandwidth, estimation_spec)
  du <- .local_du(fit$fitted$u, fit$u_grid[[grid_id]], fit$bandwidth, estimation_spec)
  ridge <- numeric(nrow(parameter_table))
  ridge_target <- numeric(nrow(parameter_table))
  ridge[parameter_table$coefficient_set == "expert"] <- if (identical(fit$family, "binomial")) {
    defaults$binomial_ridge
  } else if (identical(fit$family, "negative-binomial")) {
    defaults$negbin_ridge
  } else {
    defaults$ridge
  }
  ridge[parameter_table$coefficient_set == "gating_contrast"] <- defaults$ridge
  if (identical(fit$family, "negative-binomial")) {
    theta_rows <- parameter_table$coefficient_set == "nuisance" &
      parameter_table$term == "log_theta"
    ridge[theta_rows] <- defaults$negbin_theta_ridge * sum(weights)
    theta_target <- defaults$negbin_theta_target
    if (is.null(theta_target)) {
      theta_target <- fit$coefficients$theta[grid_id, ]
    } else if (length(theta_target) == 1L) {
      theta_target <- rep(theta_target, k)
    }
    theta_target <- pmin(pmax(as.numeric(theta_target), defaults$negbin_theta_min), defaults$negbin_theta_max)
    ridge_target[theta_rows] <- log(theta_target)
  }
  list(
    fit = fit,
    grid_id = grid_id,
    u0 = fit$u_grid[[grid_id]],
    family = fit$family,
    k = k,
    q = q,
    p = p,
    parameter_table = parameter_table,
    ridge = ridge,
    ridge_target = ridge_target,
    y = response$y,
    trials = response$trials,
    z_design = fit$fitted$z_design,
    x_design = fit$fitted$x_design,
    offset = fit$fitted$expert_offset %||% rep(0, length(response$y)),
    u = fit$fitted$u,
    bandwidth = fit$bandwidth,
    estimation_spec = estimation_spec,
    weights = weights,
    du = du,
    control = defaults
  )
}

.vcmoe_dev_parameter_vector <- function(fit, grid_id) {
  spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
  expert <- matrix(fit$coefficients$expert[grid_id, , ], nrow = spec$k, ncol = spec$q)
  expert_slope <- matrix(fit$coefficients$expert_slope[grid_id, , ], nrow = spec$k, ncol = spec$q)
  gating <- matrix(fit$coefficients$gating[grid_id, , ], nrow = spec$k, ncol = spec$p)
  gating_slope <- matrix(fit$coefficients$gating_slope[grid_id, , ], nrow = spec$k, ncol = spec$p)
  contrasts <- .vcmoe_dev_gating_contrasts(gating, gating_slope)
  par <- c(
    as.vector(t(expert)),
    as.vector(t(expert_slope)),
    as.vector(t(contrasts$intercept)),
    as.vector(t(contrasts$slope))
  )
  if (identical(fit$family, "gaussian")) {
    sigma_slope <- if (!is.null(fit$coefficients$sigma_slope)) {
      as.numeric(fit$coefficients$sigma_slope[grid_id, ])
    } else {
      rep(0, spec$k)
    }
    par <- c(
      par,
      log(pmax(fit$coefficients$sigma[grid_id, ], spec$control$min_sigma)),
      sigma_slope
    )
  } else if (identical(fit$family, "negative-binomial")) {
    par <- c(par, log(pmax(fit$coefficients$theta[grid_id, ], spec$control$negbin_theta_min)))
  }
  names(par) <- spec$parameter_table$parameter
  par
}

.vcmoe_dev_unpack_parameter <- function(par, spec) {
  k <- spec$k
  q <- spec$q
  p <- spec$p
  pos <- 1L
  take_matrix <- function(nrow, ncol) {
    n <- nrow * ncol
    out <- matrix(par[pos:(pos + n - 1L)], nrow = nrow, ncol = ncol, byrow = TRUE)
    pos <<- pos + n
    out
  }
  expert_intercept <- take_matrix(k, q)
  expert_slope <- take_matrix(k, q)
  gating_rows <- if (k == 2L) 1L else k - 1L
  gating_intercept_contrast <- take_matrix(gating_rows, p)
  gating_slope_contrast <- take_matrix(gating_rows, p)
  gating <- .vcmoe_dev_gating_from_contrasts(gating_intercept_contrast, gating_slope_contrast, k)
  sigma <- NULL
  sigma_slope <- NULL
  theta <- NULL
  if (identical(spec$family, "gaussian")) {
    sigma <- exp(par[pos:(pos + k - 1L)])
    pos <- pos + k
    sigma_slope <- par[pos:(pos + k - 1L)]
    pos <- pos + k
  } else if (identical(spec$family, "negative-binomial")) {
    theta <- exp(par[pos:(pos + k - 1L)])
    pos <- pos + k
  }
  list(
    expert_coef = cbind(expert_intercept, expert_slope),
    gating_coef = cbind(gating$intercept, gating$slope),
    sigma = sigma,
    sigma_slope = sigma_slope,
    theta = theta
  )
}

.vcmoe_dev_family_loglik_matrix <- function(params, spec) {
  expert_design <- .make_expert_design(spec$z_design, spec$du)
  gating_design <- .make_expert_design(spec$x_design, spec$du)
  prior <- .gating_prob(gating_design, params$gating_coef)
  log_terms <- matrix(NA_real_, nrow = length(spec$y), ncol = spec$k)
  eta_clipped <- FALSE
  for (component in seq_len(spec$k)) {
    eta <- as.numeric(expert_design %*% params$expert_coef[component, ])
    if (identical(spec$family, "gaussian")) {
      sigma <- exp(log(pmax(params$sigma[component], spec$control$min_sigma)) +
        params$sigma_slope[component] * spec$du)
      sigma <- pmin(pmax(sigma, spec$control$min_sigma), spec$control$max_sigma)
      expert_loglik <- stats::dnorm(spec$y, mean = eta,
                                    sd = sigma,
                                    log = TRUE)
    } else if (identical(spec$family, "binomial")) {
      expert_loglik <- stats::dbinom(spec$y, size = spec$trials,
                                     prob = stats::plogis(eta), log = TRUE)
    } else {
      eta_nb <- spec$offset + eta
      eta_clipped <- eta_clipped || any(eta_nb <= spec$control$negbin_eta_min |
        eta_nb >= spec$control$negbin_eta_max)
      mu <- .bounded_exp(eta_nb, spec$control)
      expert_loglik <- stats::dnbinom(
        spec$y,
        size = pmin(pmax(params$theta[component], spec$control$negbin_theta_min),
                    spec$control$negbin_theta_max),
        mu = mu,
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(prior[, component], 1e-12)) + expert_loglik
  }
  attr(log_terms, "eta_clipped") <- eta_clipped
  log_terms
}

.vcmoe_dev_local_loglik <- function(fit, data = NULL, grid_id, par,
                                    per_observation = FALSE, penalized = FALSE) {
  spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
  params <- .vcmoe_dev_unpack_parameter(par, spec)
  log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
  values <- spec$weights * .row_log_sum_exp(log_terms)
  if (isTRUE(per_observation)) {
    return(values)
  }
  out <- sum(values)
  if (isTRUE(penalized)) {
    out <- out - 0.5 * sum(spec$ridge * (par - spec$ridge_target)^2)
  }
  out
}

.vcmoe_dev_score_matrix <- function(par, spec, eps = 1e-5) {
  base_values <- {
    params <- .vcmoe_dev_unpack_parameter(par, spec)
    log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
    spec$weights * .row_log_sum_exp(log_terms)
  }
  scores <- matrix(NA_real_, nrow = length(base_values), ncol = length(par))
  for (j in seq_along(par)) {
    step <- eps * pmax(1, abs(par[[j]]))
    plus <- par
    minus <- par
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step
    params_plus <- .vcmoe_dev_unpack_parameter(plus, spec)
    params_minus <- .vcmoe_dev_unpack_parameter(minus, spec)
    values_plus <- spec$weights * .row_log_sum_exp(.vcmoe_dev_family_loglik_matrix(params_plus, spec))
    values_minus <- spec$weights * .row_log_sum_exp(.vcmoe_dev_family_loglik_matrix(params_minus, spec))
    scores[, j] <- (values_plus - values_minus) / (2 * step)
  }
  colnames(scores) <- names(par)
  scores
}

.vcmoe_dev_score_balance <- function(score) {
  score_sum <- colSums(score)
  score_scale <- sqrt(colSums(score^2))
  standardized <- abs(score_sum) / pmax(score_scale, sqrt(.Machine$double.eps))
  list(
    score_sum_max = max(abs(score_sum)),
    score_imbalance_max = max(standardized)
  )
}

.vcmoe_dev_hessian <- function(par, spec, eps = 1e-4) {
  d <- length(par)
  objective <- function(x) {
    -.vcmoe_dev_local_loglik(spec$fit, grid_id = spec$grid_id, par = x, penalized = TRUE)
  }
  hessian <- matrix(NA_real_, nrow = d, ncol = d, dimnames = list(names(par), names(par)))
  f0 <- objective(par)
  for (i in seq_len(d)) {
    step_i <- eps * pmax(1, abs(par[[i]]))
    ei <- rep(0, d)
    ei[[i]] <- step_i
    f_plus <- objective(par + ei)
    f_minus <- objective(par - ei)
    hessian[i, i] <- (f_plus - 2 * f0 + f_minus) / (step_i^2)
    if (i < d) {
      for (j in seq.int(i + 1L, d)) {
        step_j <- eps * pmax(1, abs(par[[j]]))
        ej <- rep(0, d)
        ej[[j]] <- step_j
        f_pp <- objective(par + ei + ej)
        f_pm <- objective(par + ei - ej)
        f_mp <- objective(par - ei + ej)
        f_mm <- objective(par - ei - ej)
        value <- (f_pp - f_pm - f_mp + f_mm) / (4 * step_i * step_j)
        hessian[i, j] <- value
        hessian[j, i] <- value
      }
    }
  }
  0.5 * (hessian + t(hessian))
}

.vcmoe_dev_matrix_condition <- function(matrix_value) {
  values <- tryCatch(eigen(0.5 * (matrix_value + t(matrix_value)), symmetric = TRUE, only.values = TRUE)$values,
                     error = function(e) rep(NA_real_, nrow(matrix_value)))
  if (!all(is.finite(values))) {
    return(list(min_eigenvalue = NA_real_, condition = Inf))
  }
  abs_values <- abs(values)
  list(
    min_eigenvalue = min(values),
    condition = max(abs_values) / max(min(abs_values), 1e-12)
  )
}

.vcmoe_dev_effective_n <- function(weights) {
  weights <- as.numeric(weights)
  weights <- weights[is.finite(weights) & weights > 0]
  if (!length(weights)) {
    return(NA_real_)
  }
  sum(weights)^2 / sum(weights^2)
}

.vcmoe_dev_adjust_meat <- function(score, inv_bread, spec, control) {
  method <- .vcmoe_dev_covariance_adjustment(control$covariance_adjustment)
  effective_n <- .vcmoe_dev_effective_n(spec$weights)
  meat <- crossprod(score)
  list(
    meat = meat,
    covariance_adjustment = method,
    covariance_factor = 1,
    local_effective_n = effective_n,
    leverage_mean = NA_real_,
    leverage_max = NA_real_
  )
}

.vcmoe_dev_score_hessian <- function(fit, data = NULL, grid_id, control = list()) {
  control <- .vcmoe_dev_control(control)
  spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
  par <- .vcmoe_dev_parameter_vector(fit, grid_id)
  if (length(par) > control$max_hessian_dim) {
    return(list(
      status = "blocked",
      block_reason = "parameter_dimension_exceeds_hessian_limit",
      parameter = par,
      spec = spec,
      score = NULL,
      bread = NULL,
      meat = NULL,
      covariance = NULL,
      condition = Inf,
      min_eigenvalue = NA_real_
    ))
  }
  score <- tryCatch(.vcmoe_dev_score_matrix(par, spec, control$finite_diff_eps), error = function(e) e)
  if (inherits(score, "error") || any(!is.finite(score))) {
    return(list(status = "blocked", block_reason = "nonfinite_score", parameter = par, spec = spec))
  }
  score_balance <- .vcmoe_dev_score_balance(score)
  bread <- tryCatch(.vcmoe_dev_hessian(par, spec, control$hessian_eps), error = function(e) e)
  if (inherits(bread, "error") || any(!is.finite(bread))) {
    return(list(status = "blocked", block_reason = "nonfinite_bread", parameter = par, spec = spec,
                score = score, meat = crossprod(score)))
  }
  condition <- .vcmoe_dev_matrix_condition(bread)
  if (!is.finite(condition$condition) || condition$condition > control$max_condition ||
      !is.finite(condition$min_eigenvalue) || condition$min_eigenvalue <= control$min_eigenvalue) {
    return(list(
      status = "blocked",
      block_reason = "ill_conditioned_bread",
      parameter = par,
      spec = spec,
      score = score,
      bread = bread,
      meat = crossprod(score),
      covariance = NULL,
      condition = condition$condition,
      min_eigenvalue = condition$min_eigenvalue
    ))
  }
  inv_bread <- tryCatch(solve(bread), error = function(e) e)
  if (inherits(inv_bread, "error") || any(!is.finite(inv_bread))) {
    return(list(status = "blocked", block_reason = "bread_inversion_failed",
                parameter = par, spec = spec, score = score, bread = bread, meat = crossprod(score),
                condition = condition$condition, min_eigenvalue = condition$min_eigenvalue))
  }
  adjusted <- .vcmoe_dev_adjust_meat(score, inv_bread, spec, control)
  meat <- adjusted$meat
  if (any(!is.finite(meat))) {
    return(list(status = "blocked", block_reason = "invalid_covariance_adjustment",
                parameter = par, spec = spec, score = score, bread = bread, meat = meat,
                condition = condition$condition, min_eigenvalue = condition$min_eigenvalue,
                covariance_adjustment = adjusted$covariance_adjustment,
                covariance_factor = adjusted$covariance_factor,
                local_effective_n = adjusted$local_effective_n,
                leverage_mean = adjusted$leverage_mean,
                leverage_max = adjusted$leverage_max))
  }
  covariance <- inv_bread %*% meat %*% inv_bread
  covariance <- 0.5 * (covariance + t(covariance))
  if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
    return(list(status = "blocked", block_reason = "invalid_covariance",
                parameter = par, spec = spec, score = score, bread = bread, meat = meat,
                covariance = covariance, condition = condition$condition,
                min_eigenvalue = condition$min_eigenvalue,
                covariance_adjustment = adjusted$covariance_adjustment,
                covariance_factor = adjusted$covariance_factor,
                local_effective_n = adjusted$local_effective_n,
                leverage_mean = adjusted$leverage_mean,
                leverage_max = adjusted$leverage_max))
  }
  list(
    status = "ok",
    block_reason = NA_character_,
    parameter = par,
    spec = spec,
    score = score,
    bread = bread,
    meat = meat,
    covariance = covariance,
    condition = condition$condition,
    min_eigenvalue = condition$min_eigenvalue,
    covariance_adjustment = adjusted$covariance_adjustment,
    covariance_factor = adjusted$covariance_factor,
    local_effective_n = adjusted$local_effective_n,
    leverage_mean = adjusted$leverage_mean,
    leverage_max = adjusted$leverage_max,
    score_sum_max = score_balance$score_sum_max,
    score_imbalance_max = score_balance$score_imbalance_max
  )
}

.vcmoe_dev_multiplier_critical <- function(covariance, level = 0.90, B = 199L, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  se <- sqrt(pmax(diag(covariance), 0))
  if (any(!is.finite(se)) || any(se <= 0)) {
    return(NA_real_)
  }
  corr <- covariance / outer(se, se)
  corr <- 0.5 * (corr + t(corr))
  eig <- tryCatch(eigen(corr, symmetric = TRUE), error = function(e) NULL)
  if (is.null(eig) || any(!is.finite(eig$values))) {
    return(NA_real_)
  }
  values <- pmax(eig$values, 0)
  transform <- eig$vectors %*% diag(sqrt(values), nrow = length(values))
  draws <- matrix(stats::rnorm(B * nrow(corr)), nrow = B)
  z <- draws %*% t(transform)
  as.numeric(stats::quantile(apply(abs(z), 1L, max), probs = level, names = FALSE, na.rm = TRUE))
}

.vcmoe_dev_epanechnikov_scb_kernel_constant <- function(kernel_constant = NULL) {
  if (is.null(kernel_constant)) {
    return(list(value = 5 / (8 * pi), source = "epanechnikov"))
  }
  if (!is.numeric(kernel_constant) || length(kernel_constant) != 1L ||
      !is.finite(kernel_constant) || kernel_constant <= 0) {
    stop("`scb_kernel_constant` must be a positive finite number.", call. = FALSE)
  }
  list(value = as.numeric(kernel_constant), source = "user_supplied")
}

.vcmoe_dev_epanechnikov_critical <- function(bandwidth, level = 0.90,
                                             domain_length = 1,
                                             kernel_constant = NULL) {
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L ||
      !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(domain_length) || length(domain_length) != 1L ||
      !is.finite(domain_length) || domain_length <= 0) {
    stop("`scb_domain_length` must be a positive finite number.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be between 0 and 1.", call. = FALSE)
  }
  h_eff <- bandwidth / domain_length
  if (!is.finite(h_eff) || h_eff <= 0 || h_eff >= 1) {
    stop("Analytic Epanechnikov SCB requires `bandwidth / scb_domain_length` in (0, 1).",
         call. = FALSE)
  }
  constant <- .vcmoe_dev_epanechnikov_scb_kernel_constant(kernel_constant)
  s <- sqrt(-2 * log(h_eff))
  dvn <- s + (log(constant$value) - 0.5 * log(log(1 / h_eff))) / s
  critical <- dvn + (log(2) - log(-log(level))) / s
  list(
    critical = as.numeric(critical),
    dvn = as.numeric(dvn),
    s = as.numeric(s),
    h_eff = as.numeric(h_eff),
    kernel_constant = constant$value,
    kernel_constant_source = constant$source,
    domain_length = as.numeric(domain_length)
  )
}

.vcmoe_dev_require_analytic_epanechnikov_fit <- function(fit) {
  .vcmoe_dev_require_fit(fit)
  metadata <- vcmoe_parameterization(fit)
  if (!identical(metadata$id, "a1_epanechnikov_scaled")) {
    stop("Analytic Epanechnikov path SCB requires `parameterization = \"a1_epanechnikov_scaled\"`.",
         call. = FALSE)
  }
  if (!identical(metadata$kernel$name, "epanechnikov") ||
      !identical(metadata$kernel$weight_normalization, "density_over_bandwidth") ||
      !identical(metadata$local_linear_basis$slope_storage, "scaled")) {
    stop("Analytic Epanechnikov path SCB requires Epanechnikov density weights and scaled local-linear slopes.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.vcmoe_dev_grid_block_reason <- function(fit, grid_id, control) {
  diagnostics <- vcmoe_diagnostics(fit)
  fit_control <- .vcmoe_default_control(fit$control %||% list())
  reasons <- character(0L)
  if (mean(fit$diagnostics$converged %in% TRUE) < control$min_converged_fraction) {
    reasons <- c(reasons, "low_converged_fraction")
  }
  if (isTRUE(fit$diagnostics$ambiguous[[grid_id]])) {
    reasons <- c(reasons, "ambiguous_label")
  }
  if (identical(fit$diagnostics$alignment_method, "sequential") &&
      is.finite(diagnostics$alignment_margin[[grid_id]]) &&
      diagnostics$alignment_margin[[grid_id]] < control$min_alignment_margin) {
    reasons <- c(reasons, "weak_sequential_alignment_margin")
  }
  if (is.finite(diagnostics$min_posterior_mean[[grid_id]]) &&
      diagnostics$min_posterior_mean[[grid_id]] < control$min_component_proportion) {
    reasons <- c(reasons, "component_collapse")
  }
  if (is.finite(diagnostics$effective_n[[grid_id]]) &&
      diagnostics$effective_n[[grid_id]] < control$min_effective_n) {
    reasons <- c(reasons, "low_effective_n")
  }
  if (identical(fit$family, "gaussian") && !is.null(fit$coefficients$sigma)) {
    if (any(fit$coefficients$sigma[grid_id, ] <= fit_control$min_sigma * 1.01)) {
      reasons <- c(reasons, "sigma_boundary")
    }
  }
  if (identical(fit$family, "negative-binomial") && !is.null(fit$coefficients$theta)) {
    theta <- fit$coefficients$theta[grid_id, ]
    if (any(theta <= fit_control$negbin_theta_min * 1.01 |
            theta >= fit_control$negbin_theta_max / 1.01)) {
      reasons <- c(reasons, "theta_boundary")
    }
    spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
    params <- .vcmoe_dev_unpack_parameter(.vcmoe_dev_parameter_vector(fit, grid_id), spec)
    log_terms <- .vcmoe_dev_family_loglik_matrix(params, spec)
    if (isTRUE(attr(log_terms, "eta_clipped"))) {
      reasons <- c(reasons, "negbin_eta_clipping")
    }
    expert_design <- .make_expert_design(spec$z_design, spec$du)
    eta_values <- vapply(seq_len(spec$k), function(component) {
      spec$offset + as.numeric(expert_design %*% params$expert_coef[component, ])
    }, numeric(length(spec$y)))
    if (any(eta_values <= fit_control$negbin_eta_min + control$negbin_eta_margin |
            eta_values >= fit_control$negbin_eta_max - control$negbin_eta_margin)) {
      reasons <- c(reasons, "negbin_eta_near_clipping")
    }
  }
  paste(unique(reasons), collapse = ";")
}

.vcmoe_dev_interval_metadata <- function(spec, control, simultaneous_critical,
                                         analytic_meta = NULL) {
  is_analytic <- identical(control$simultaneous_method, "analytic_epanechnikov_path")
  engine_id <- spec$fit$engine_id %||% "local_grid_em"
  joint_path <- identical(engine_id, "joint_path_em")
  data.frame(
    engine_id = engine_id,
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
    parameterization = spec$estimation_spec$id,
    kernel = spec$estimation_spec$kernel,
    weight_normalization = if (identical(spec$estimation_spec$weight_normalization, "density")) {
      "density_over_bandwidth"
    } else {
      "mean_one_per_grid"
    },
    local_basis = spec$estimation_spec$local_basis,
    slope_storage = spec$estimation_spec$slope_storage,
    simultaneous_method = control$simultaneous_method,
    covariance_adjustment = control$covariance_adjustment,
    simultaneous_scope = if (is_analytic) {
      "finite_grid_coefficient_function_path"
    } else {
      "finite_grid_parameter_vector"
    },
    simultaneous_critical = simultaneous_critical,
    simultaneous_dvn = if (is_analytic) analytic_meta$dvn else NA_real_,
    scb_h_eff = if (is_analytic) analytic_meta$h_eff else NA_real_,
    scb_kernel_constant = if (is_analytic) analytic_meta$kernel_constant else NA_real_,
    scb_kernel_constant_source = if (is_analytic) {
      analytic_meta$kernel_constant_source
    } else {
      NA_character_
    },
    scb_domain_length = if (is_analytic) analytic_meta$domain_length else NA_real_,
    bias_correction = "none",
    boundary_correction = "none",
    stringsAsFactors = FALSE
  )
}

.vcmoe_dev_add_interval_metadata <- function(rows, spec, control,
                                             simultaneous_critical,
                                             analytic_meta = NULL) {
  meta <- .vcmoe_dev_interval_metadata(spec, control, simultaneous_critical, analytic_meta)
  cbind(rows, meta[rep(1L, nrow(rows)), , drop = FALSE])
}

.vcmoe_dev_blocked_interval_rows <- function(par, spec, level, reason, diagnostic,
                                             control, analytic_meta = NULL) {
  table <- spec$parameter_table
  rows <- data.frame(
    grid_id = spec$grid_id,
    u = spec$u0,
    parameter = table$parameter,
    coefficient_set = table$coefficient_set,
    component = table$component,
    term = table$term,
    block = table$block,
    estimate = as.numeric(par),
    se = NA_real_,
    pointwise_lower = NA_real_,
    pointwise_upper = NA_real_,
    simultaneous_lower = NA_real_,
    simultaneous_upper = NA_real_,
    level = level,
    status = "blocked",
    block_reason = reason,
    hessian_condition = diagnostic$condition %||% NA_real_,
    min_eigenvalue = diagnostic$min_eigenvalue %||% NA_real_,
    covariance_factor = diagnostic$covariance_factor %||% NA_real_,
    local_effective_n = diagnostic$local_effective_n %||% NA_real_,
    leverage_mean = diagnostic$leverage_mean %||% NA_real_,
    leverage_max = diagnostic$leverage_max %||% NA_real_,
    score_sum_max = diagnostic$score_sum_max %||% NA_real_,
    score_imbalance_max = diagnostic$score_imbalance_max %||% NA_real_,
    stringsAsFactors = FALSE
  )
  .vcmoe_dev_add_interval_metadata(rows, spec, control, NA_real_, analytic_meta)
}

.vcmoe_dev_intervals <- function(fit, data = NULL, level = 0.90, strict = TRUE, control = list()) {
  .vcmoe_dev_require_fit(fit)
  control <- .vcmoe_dev_control(utils::modifyList(control %||% list(), list(level = level, strict = strict)))
  analytic_meta <- NULL
  if (identical(control$simultaneous_method, "analytic_epanechnikov_path")) {
    .vcmoe_dev_require_analytic_epanechnikov_fit(fit)
    domain_length <- control$scb_domain_length %||% fit$u_scaling$domain_length %||%
      vcmoe_parameterization(fit)$u_scaling$domain_length %||% 1
    analytic_meta <- .vcmoe_dev_epanechnikov_critical(
      bandwidth = fit$bandwidth,
      level = level,
      domain_length = domain_length,
      kernel_constant = control$scb_kernel_constant
    )
  }
  interval_rows <- list()
  diagnostic_rows <- list()
  for (grid_id in seq_along(fit$u_grid)) {
    spec <- .vcmoe_dev_parameter_spec(fit, grid_id)
    par <- .vcmoe_dev_parameter_vector(fit, grid_id)
    pre_reason <- .vcmoe_dev_grid_block_reason(fit, grid_id, control)
    if (nzchar(pre_reason) && isTRUE(strict)) {
      diagnostic <- list(condition = NA_real_, min_eigenvalue = NA_real_)
      interval_rows[[grid_id]] <- .vcmoe_dev_blocked_interval_rows(
        par, spec, level, pre_reason, diagnostic, control, analytic_meta
      )
      diagnostic_rows[[grid_id]] <- data.frame(
        grid_id = grid_id, u = fit$u_grid[[grid_id]], status = "blocked",
        block_reason = pre_reason, parameter_count = length(par),
        hessian_condition = NA_real_, min_eigenvalue = NA_real_,
        covariance_factor = NA_real_, local_effective_n = NA_real_,
        leverage_mean = NA_real_, leverage_max = NA_real_,
        score_sum_max = NA_real_, score_imbalance_max = NA_real_,
        simultaneous_method = control$simultaneous_method,
        covariance_adjustment = control$covariance_adjustment,
        simultaneous_critical = if (is.null(analytic_meta)) NA_real_ else analytic_meta$critical,
        stringsAsFactors = FALSE
      )
      next
    }
    sh <- .vcmoe_dev_score_hessian(fit, data, grid_id, control)
    if (!identical(sh$status, "ok")) {
      reason <- paste(c(pre_reason, sh$block_reason), collapse = ";")
      reason <- sub("^;", "", reason)
      interval_rows[[grid_id]] <- .vcmoe_dev_blocked_interval_rows(
        par, spec, level, reason, sh, control, analytic_meta
      )
      diagnostic_rows[[grid_id]] <- data.frame(
        grid_id = grid_id, u = fit$u_grid[[grid_id]], status = "blocked",
        block_reason = reason, parameter_count = length(par),
        hessian_condition = sh$condition %||% NA_real_,
        min_eigenvalue = sh$min_eigenvalue %||% NA_real_,
        covariance_factor = sh$covariance_factor %||% NA_real_,
        local_effective_n = sh$local_effective_n %||% NA_real_,
        leverage_mean = sh$leverage_mean %||% NA_real_,
        leverage_max = sh$leverage_max %||% NA_real_,
        score_sum_max = sh$score_sum_max %||% NA_real_,
        score_imbalance_max = sh$score_imbalance_max %||% NA_real_,
        simultaneous_method = control$simultaneous_method,
        covariance_adjustment = control$covariance_adjustment,
        simultaneous_critical = if (is.null(analytic_meta)) NA_real_ else analytic_meta$critical,
        stringsAsFactors = FALSE
      )
      next
    }
    se <- sqrt(diag(sh$covariance))
    z <- stats::qnorm(1 - (1 - level) / 2)
    multiplier <- if (identical(control$simultaneous_method, "analytic_epanechnikov_path")) {
      analytic_meta$critical
    } else {
      .vcmoe_dev_multiplier_critical(
        sh$covariance,
        level = level,
        B = control$multiplier_B,
        seed = if (is.null(control$seed)) NULL else control$seed + grid_id
      )
    }
    if (!is.finite(multiplier)) {
      multiplier <- NA_real_
    }
    table <- spec$parameter_table
    rows <- data.frame(
      grid_id = grid_id,
      u = spec$u0,
      parameter = table$parameter,
      coefficient_set = table$coefficient_set,
      component = table$component,
      term = table$term,
      block = table$block,
      estimate = as.numeric(sh$parameter),
      se = se,
      pointwise_lower = as.numeric(sh$parameter) - z * se,
      pointwise_upper = as.numeric(sh$parameter) + z * se,
      simultaneous_lower = if (is.finite(multiplier)) as.numeric(sh$parameter) - multiplier * se else NA_real_,
      simultaneous_upper = if (is.finite(multiplier)) as.numeric(sh$parameter) + multiplier * se else NA_real_,
      level = level,
      status = if (nzchar(pre_reason)) "warning" else "ok",
      block_reason = if (nzchar(pre_reason)) pre_reason else NA_character_,
      hessian_condition = sh$condition,
      min_eigenvalue = sh$min_eigenvalue,
      covariance_factor = sh$covariance_factor %||% NA_real_,
      local_effective_n = sh$local_effective_n %||% NA_real_,
      leverage_mean = sh$leverage_mean %||% NA_real_,
      leverage_max = sh$leverage_max %||% NA_real_,
      score_sum_max = sh$score_sum_max %||% NA_real_,
      score_imbalance_max = sh$score_imbalance_max %||% NA_real_,
      stringsAsFactors = FALSE
    )
    interval_rows[[grid_id]] <- .vcmoe_dev_add_interval_metadata(
      rows, spec, control, multiplier, analytic_meta
    )
    diagnostic_rows[[grid_id]] <- data.frame(
      grid_id = grid_id, u = fit$u_grid[[grid_id]],
      status = if (nzchar(pre_reason)) "warning" else "ok",
      block_reason = if (nzchar(pre_reason)) pre_reason else NA_character_,
      parameter_count = length(par),
      hessian_condition = sh$condition,
      min_eigenvalue = sh$min_eigenvalue,
      covariance_factor = sh$covariance_factor %||% NA_real_,
      local_effective_n = sh$local_effective_n %||% NA_real_,
      leverage_mean = sh$leverage_mean %||% NA_real_,
      leverage_max = sh$leverage_max %||% NA_real_,
      score_sum_max = sh$score_sum_max %||% NA_real_,
      score_imbalance_max = sh$score_imbalance_max %||% NA_real_,
      simultaneous_method = control$simultaneous_method,
      covariance_adjustment = control$covariance_adjustment,
      simultaneous_critical = multiplier,
      stringsAsFactors = FALSE
    )
  }
  list(
    intervals = do.call(rbind, interval_rows),
    diagnostics = do.call(rbind, diagnostic_rows)
  )
}

.vcmoe_dev_null_constant_fit <- function(fit) {
  null_fit <- fit
  expert_mean <- apply(fit$coefficients$expert, c(2L, 3L), mean)
  expert_slope_mean <- 0 * apply(fit$coefficients$expert_slope, c(2L, 3L), mean)
  gating_mean <- apply(fit$coefficients$gating, c(2L, 3L), mean)
  gating_slope_mean <- 0 * apply(fit$coefficients$gating_slope, c(2L, 3L), mean)
  for (grid_id in seq_along(fit$u_grid)) {
    null_fit$coefficients$expert[grid_id, , ] <- expert_mean
    null_fit$coefficients$expert_slope[grid_id, , ] <- expert_slope_mean
    null_fit$coefficients$gating[grid_id, , ] <- gating_mean
    null_fit$coefficients$gating_slope[grid_id, , ] <- gating_slope_mean
  }
  if (!is.null(fit$coefficients$sigma)) {
    log_sigma <- colMeans(log(pmax(fit$coefficients$sigma, .vcmoe_default_control(list())$min_sigma)))
    for (grid_id in seq_along(fit$u_grid)) {
      null_fit$coefficients$sigma[grid_id, ] <- exp(log_sigma)
    }
    if (!is.null(fit$coefficients$sigma_slope)) {
      null_fit$coefficients$sigma_slope[,] <- 0
    }
  }
  if (!is.null(fit$coefficients$theta)) {
    log_theta <- colMeans(log(pmax(fit$coefficients$theta, .vcmoe_default_control(list())$negbin_theta_min)))
    for (grid_id in seq_along(fit$u_grid)) {
      null_fit$coefficients$theta[grid_id, ] <- exp(log_theta)
    }
  }
  null_fit
}

.vcmoe_dev_global_pseudologlik <- function(fit) {
  sum(vapply(seq_along(fit$u_grid), function(grid_id) {
    par <- .vcmoe_dev_parameter_vector(fit, grid_id)
    .vcmoe_dev_local_loglik(fit, grid_id = grid_id, par = par, penalized = FALSE)
  }, numeric(1L)))
}

.vcmoe_dev_glrt_like <- function(full_fit, null_fit = NULL, data = NULL, calibration = "none", control = list()) {
  .vcmoe_dev_require_fit(full_fit)
  control <- .vcmoe_dev_control(control)
  block_reasons <- vapply(seq_along(full_fit$u_grid), function(grid_id) {
    .vcmoe_dev_grid_block_reason(full_fit, grid_id, control)
  }, character(1L))
  block_reason <- paste(unique(unlist(strsplit(block_reasons[nzchar(block_reasons)], ";", fixed = TRUE))),
                        collapse = ";")
  if (nzchar(block_reason)) {
    return(data.frame(
      statistic = NA_real_,
      full_loglik = NA_real_,
      null_loglik = NA_real_,
      calibration = calibration,
      calibrated_p = NA_real_,
      calibration_status = "blocked",
      status = "blocked",
      block_reason = block_reason,
      stringsAsFactors = FALSE
    ))
  }
  if (is.null(null_fit)) {
    null_fit <- .vcmoe_dev_null_constant_fit(full_fit)
  }
  full_loglik <- .vcmoe_dev_global_pseudologlik(full_fit)
  null_loglik <- .vcmoe_dev_global_pseudologlik(null_fit)
  statistic <- max(0, 2 * (full_loglik - null_loglik))
  data.frame(
    statistic = statistic,
    full_loglik = full_loglik,
    null_loglik = null_loglik,
    calibration = calibration,
    calibrated_p = NA_real_,
    calibration_status = if (identical(calibration, "none")) "not_run" else "not_implemented_dev_only",
    status = "ok",
    block_reason = NA_character_,
    stringsAsFactors = FALSE
  )
}
