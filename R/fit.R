.make_expert_design <- function(z_design, du) {
  cbind(z_design, z_design * as.numeric(du))
}

.fit_binomial_logistic <- function(success, trials, design, weights, start, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= 1e-10) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start <- as.numeric(start)
  ridge <- control$binomial_ridge %||% control$ridge
  objective <- function(par) {
    eta <- as.numeric(design %*% par)
    prob <- stats::plogis(eta)
    prob <- pmin(pmax(prob, 1e-10), 1 - 1e-10)
    value <- -sum(weights * (
      success * log(prob) + (trials - success) * log1p(-prob)
    ))
    value + 0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    eta <- as.numeric(design %*% par)
    prob <- stats::plogis(eta)
    as.numeric(-crossprod(design, weights * (success - trials * prob)) + ridge * par)
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$gating_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par)
  list(
    coef = coef,
    convergence = opt$convergence,
    gradient_norm = sqrt(sum(gradient(coef)^2)),
    value = opt$value
  )
}

.bounded_exp <- function(eta, control) {
  exp(pmin(pmax(eta, control$negbin_eta_min), control$negbin_eta_max))
}

.fit_negbin_loglink <- function(y, design, offset, weights, start, theta, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start <- as.numeric(start)
  theta <- pmin(pmax(theta, control$negbin_theta_min), control$negbin_theta_max)
  ridge <- control$negbin_ridge %||% control$ridge
  objective <- function(par) {
    eta <- offset + as.numeric(design %*% par)
    mu <- .bounded_exp(eta, control)
    value <- -sum(weights * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
    value + 0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    eta_raw <- offset + as.numeric(design %*% par)
    eta <- pmin(pmax(eta_raw, control$negbin_eta_min), control$negbin_eta_max)
    mu <- exp(eta)
    score_eta <- y - mu * (y + theta) / (mu + theta)
    score_eta <- score_eta * as.numeric(eta_raw > control$negbin_eta_min & eta_raw < control$negbin_eta_max)
    as.numeric(-crossprod(design, weights * score_eta) + ridge * par)
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$negbin_mstep_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par)
  list(
    coef = coef,
    convergence = opt$convergence,
    gradient_norm = sqrt(sum(gradient(coef)^2)),
    value = opt$value
  )
}

.fit_negbin_theta <- function(y, mu, weights, previous_theta, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(pmin(pmax(previous_theta, control$negbin_theta_min), control$negbin_theta_max))
  }
  lower <- log(control$negbin_theta_min)
  upper <- log(control$negbin_theta_max)
  theta_target <- control$negbin_theta_target %||% previous_theta
  theta_target <- as.numeric(theta_target[[1L]])
  if (!is.finite(theta_target) || theta_target <= 0) {
    theta_target <- previous_theta
  }
  theta_target <- pmin(pmax(theta_target, control$negbin_theta_min), control$negbin_theta_max)
  log_theta_target <- log(theta_target)
  theta_ridge <- as.numeric(control$negbin_theta_ridge %||% 0)
  if (!is.finite(theta_ridge) || theta_ridge < 0) {
    theta_ridge <- 0
  }
  penalty_scale <- sum(weights)
  objective <- function(log_theta) {
    theta <- exp(log_theta)
    value <- -sum(weights * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
    value + 0.5 * theta_ridge * penalty_scale * (log_theta - log_theta_target)^2
  }
  opt <- tryCatch(
    stats::optimize(objective, interval = c(lower, upper)),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$minimum)) {
    return(pmin(pmax(previous_theta, control$negbin_theta_min), control$negbin_theta_max))
  }
  exp(opt$minimum)
}

.gaussian_sigma_from_eta <- function(log_sigma_eta, control) {
  exp(pmin(pmax(log_sigma_eta, log(control$min_sigma)), log(control$max_sigma)))
}

.gaussian_component_sigma <- function(params, component, du, control) {
  slope <- if (!is.null(params$sigma_slope)) params$sigma_slope[[component]] else 0
  .gaussian_sigma_from_eta(log(pmax(params$sigma[[component]], control$min_sigma)) + slope * du, control)
}

.fit_gaussian_locallinear <- function(y, design, du, weights, start_coef,
                                      start_sigma, start_sigma_slope, control) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= control$min_component_weight) {
    return(list(
      coef = start_coef,
      sigma = pmin(pmax(start_sigma, control$min_sigma), control$max_sigma),
      sigma_slope = start_sigma_slope,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  start_coef <- as.numeric(start_coef)
  start_sigma <- pmin(pmax(start_sigma, control$min_sigma), control$max_sigma)
  start_sigma_slope <- as.numeric(start_sigma_slope %||% 0)
  if (!is.finite(start_sigma_slope)) {
    start_sigma_slope <- 0
  }
  start <- c(start_coef, log(start_sigma), start_sigma_slope)
  coef_len <- length(start_coef)
  ridge <- control$ridge
  sigma_ridge <- control$gaussian_sigma_ridge %||% control$ridge
  objective <- function(par) {
    if (any(!is.finite(par))) {
      return(1e100)
    }
    beta <- par[seq_len(coef_len)]
    log_sigma_eta <- par[[coef_len + 1L]] + par[[coef_len + 2L]] * du
    sigma <- .gaussian_sigma_from_eta(log_sigma_eta, control)
    mu <- as.numeric(design %*% beta)
    value <- -sum(weights * stats::dnorm(y, mean = mu, sd = sigma, log = TRUE))
    value + 0.5 * ridge * sum(beta^2) + 0.5 * sigma_ridge * par[[coef_len + 2L]]^2
  }
  opt <- tryCatch(
    stats::optim(
      start,
      objective,
      method = "BFGS",
      control = list(maxit = control$gaussian_mstep_maxit)
    ),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(list(
      coef = start_coef,
      sigma = start_sigma,
      sigma_slope = start_sigma_slope,
      convergence = NA_integer_,
      gradient_norm = NA_real_,
      value = NA_real_
    ))
  }
  coef <- as.numeric(opt$par[seq_len(coef_len)])
  sigma <- .gaussian_sigma_from_eta(opt$par[[coef_len + 1L]], control)
  sigma_slope <- as.numeric(opt$par[[coef_len + 2L]])
  list(
    coef = coef,
    sigma = sigma,
    sigma_slope = sigma_slope,
    convergence = opt$convergence,
    gradient_norm = NA_real_,
    value = opt$value
  )
}

.mstep_gating <- function(gating_design, posterior, weights, beta_full, control) {
  k <- ncol(posterior)
  n_terms <- ncol(gating_design)
  weights <- pmax(weights, 0)
  ridge <- control$ridge

  if (k == 2L) {
    start <- beta_full[1L, ] - beta_full[2L, ]
    objective <- function(par) {
      eta <- as.numeric(gating_design %*% par)
      p1 <- stats::plogis(eta)
      p1 <- pmin(pmax(p1, 1e-10), 1 - 1e-10)
      value <- -sum(weights * (
        posterior[, 1L] * log(p1) + posterior[, 2L] * log1p(-p1)
      ))
      value + 0.5 * ridge * sum(par^2)
    }
    gradient <- function(par) {
      eta <- as.numeric(gating_design %*% par)
      p1 <- stats::plogis(eta)
      as.numeric(-crossprod(gating_design, weights * (posterior[, 1L] - p1)) + ridge * par)
    }
    opt <- tryCatch(
      stats::optim(start, objective, gradient, method = "BFGS",
                   control = list(maxit = control$gating_maxit)),
      error = function(e) NULL
    )
    if (is.null(opt) || any(!is.finite(opt$par))) {
      return(beta_full)
    }
    beta <- opt$par
    return(rbind(beta / 2, -beta / 2))
  }

  start_matrix <- beta_full[-1L, , drop = FALSE] -
    matrix(beta_full[1L, ], nrow = k - 1L, ncol = n_terms, byrow = TRUE)
  start <- as.vector(t(start_matrix))

  objective <- function(par) {
    theta <- matrix(par, nrow = k - 1L, byrow = TRUE)
    eta <- cbind(0, gating_design %*% t(theta))
    prob <- .row_softmax(eta)
    -sum(weights * rowSums(posterior * log(pmax(prob, 1e-12)))) +
      0.5 * ridge * sum(par^2)
  }
  gradient <- function(par) {
    theta <- matrix(par, nrow = k - 1L, byrow = TRUE)
    eta <- cbind(0, gating_design %*% t(theta))
    prob <- .row_softmax(eta)
    grad <- matrix(0, nrow = k - 1L, ncol = n_terms)
    for (component in seq_len(k - 1L)) {
      grad[component, ] <- -crossprod(
        gating_design,
        weights * (posterior[, component + 1L] - prob[, component + 1L])
      )
    }
    as.vector(t(grad)) + ridge * par
  }
  opt <- tryCatch(
    stats::optim(start, objective, gradient, method = "BFGS",
                 control = list(maxit = control$gating_maxit)),
    error = function(e) NULL
  )
  if (is.null(opt) || any(!is.finite(opt$par))) {
    return(beta_full)
  }

  theta <- matrix(opt$par, nrow = k - 1L, byrow = TRUE)
  beta_baseline <- rbind(rep(0, n_terms), theta)
  .center_logits(beta_baseline)
}

.m_step <- function(y, trials, family, z_design, x_design, du, offset, weights, posterior, previous, control) {
  k <- ncol(posterior)
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)

  expert_coef <- matrix(0, nrow = k, ncol = 2L * q)
  sigma <- if (identical(family, "gaussian")) numeric(k) else NULL
  sigma_slope <- if (identical(family, "gaussian")) numeric(k) else NULL
  theta <- if (identical(family, "negative-binomial")) numeric(k) else NULL
  expert_convergence <- rep(NA_integer_, k)
  expert_gradient_norm <- rep(NA_real_, k)
  expert_coef_norm <- rep(NA_real_, k)
  for (component in seq_len(k)) {
    wc <- weights * posterior[, component]
    if (identical(family, "gaussian")) {
      gaussian_fit <- .fit_gaussian_locallinear(
        y,
        expert_design,
        du,
        wc,
        previous$expert_coef[component, ],
        previous$sigma[component],
        previous$sigma_slope[component] %||% 0,
        control
      )
      coef_c <- gaussian_fit$coef
      sigma[component] <- gaussian_fit$sigma
      sigma_slope[component] <- gaussian_fit$sigma_slope
      expert_convergence[component] <- gaussian_fit$convergence
      expert_gradient_norm[component] <- gaussian_fit$gradient_norm
    } else if (identical(family, "binomial")) {
      expert_fit <- .fit_binomial_logistic(
        y,
        trials,
        expert_design,
        wc,
        previous$expert_coef[component, ],
        control
      )
      coef_c <- expert_fit$coef
      expert_convergence[component] <- expert_fit$convergence
      expert_gradient_norm[component] <- expert_fit$gradient_norm
    } else {
      expert_fit <- .fit_negbin_loglink(
        y,
        expert_design,
        offset,
        wc,
        previous$expert_coef[component, ],
        previous$theta[component],
        control
      )
      coef_c <- expert_fit$coef
      eta_c <- offset + as.numeric(expert_design %*% coef_c)
      mu_c <- .bounded_exp(eta_c, control)
      theta[component] <- .fit_negbin_theta(
        y,
        mu_c,
        wc,
        previous$theta[component],
        control
      )
      expert_convergence[component] <- expert_fit$convergence
      expert_gradient_norm[component] <- expert_fit$gradient_norm
    }
    expert_coef[component, ] <- coef_c
    expert_coef_norm[component] <- sqrt(sum(coef_c^2))
  }

  beta_full <- .mstep_gating(gating_design, posterior, weights, previous$gating_coef, control)
  list(
    expert_coef = expert_coef,
    gating_coef = beta_full,
    sigma = sigma,
    sigma_slope = sigma_slope,
    theta = theta,
    expert_convergence = expert_convergence,
    expert_gradient_norm = expert_gradient_norm,
    expert_coef_norm = expert_coef_norm
  )
}

.local_loglik <- function(y, trials, family, z_design, x_design, du, offset, weights, params, control) {
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)
  pi <- .gating_prob(gating_design, params$gating_coef)
  k <- nrow(params$expert_coef)
  log_terms <- matrix(0, nrow = length(y), ncol = k)
  for (component in seq_len(k)) {
    eta <- offset + as.numeric(expert_design %*% params$expert_coef[component, ])
    if (identical(family, "gaussian")) {
      expert_loglik <- stats::dnorm(
        y,
        mean = eta,
        sd = .gaussian_component_sigma(params, component, du, control),
        log = TRUE
      )
    } else if (identical(family, "binomial")) {
      expert_loglik <- stats::dbinom(y, size = trials, prob = stats::plogis(eta), log = TRUE)
    } else {
      expert_loglik <- stats::dnbinom(
        y,
        size = params$theta[component],
        mu = .bounded_exp(eta, control),
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(pi[, component], 1e-12)) + expert_loglik
  }
  sum(weights * .row_log_sum_exp(log_terms))
}

.e_step <- function(y, trials, family, z_design, x_design, du, offset, params, control) {
  expert_design <- .make_expert_design(z_design, du)
  gating_design <- .make_expert_design(x_design, du)
  pi <- .gating_prob(gating_design, params$gating_coef)
  k <- nrow(params$expert_coef)
  log_terms <- matrix(0, nrow = length(y), ncol = k)
  for (component in seq_len(k)) {
    eta <- offset + as.numeric(expert_design %*% params$expert_coef[component, ])
    if (identical(family, "gaussian")) {
      expert_loglik <- stats::dnorm(
        y,
        mean = eta,
        sd = .gaussian_component_sigma(params, component, du, control),
        log = TRUE
      )
    } else if (identical(family, "binomial")) {
      expert_loglik <- stats::dbinom(y, size = trials, prob = stats::plogis(eta), log = TRUE)
    } else {
      expert_loglik <- stats::dnbinom(
        y,
        size = params$theta[component],
        mu = .bounded_exp(eta, control),
        log = TRUE
      )
    }
    log_terms[, component] <- log(pmax(pi[, component], 1e-12)) + expert_loglik
  }
  log_norm <- .row_log_sum_exp(log_terms)
  posterior <- exp(log_terms - log_norm)
  posterior / rowSums(posterior)
}

.params_to_local <- function(params, posterior, z_design, x_design, du, loglik, iterations, converged) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  gating_design <- .make_expert_design(x_design, du)
  list(
    expert_intercept = params$expert_coef[, seq_len(q), drop = FALSE],
    expert_slope = params$expert_coef[, q + seq_len(q), drop = FALSE],
    gating_intercept = params$gating_coef[, seq_len(p), drop = FALSE],
    gating_slope = params$gating_coef[, p + seq_len(p), drop = FALSE],
    sigma = params$sigma,
    sigma_slope = params$sigma_slope,
    theta = params$theta,
    posterior = posterior,
    pi = .gating_prob(gating_design, params$gating_coef),
    loglik = loglik,
    iterations = iterations,
    converged = converged,
    expert_convergence = params$expert_convergence,
    expert_gradient_norm = params$expert_gradient_norm,
    expert_coef_norm = params$expert_coef_norm
  )
}

.local_to_params <- function(local) {
  list(
    expert_coef = cbind(local$expert_intercept, local$expert_slope),
    gating_coef = cbind(local$gating_intercept, local$gating_slope),
    sigma = local$sigma,
    sigma_slope = local$sigma_slope,
    theta = local$theta
  )
}

.initial_tau <- function(initial_y, k, random = FALSE) {
  n <- length(initial_y)
  if (random) {
    raw <- matrix(stats::runif(n * k), nrow = n, ncol = k)
    return(raw / rowSums(raw))
  }

  groups <- tryCatch(stats::kmeans(initial_y, centers = k, nstart = 10)$cluster,
                     error = function(e) NULL)
  if (is.null(groups) || any(is.na(groups))) {
    groups <- cut(rank(initial_y, ties.method = "first"),
                  breaks = k, labels = FALSE, include.lowest = TRUE)
  }
  tau <- matrix(0.05 / max(k - 1L, 1L), nrow = n, ncol = k)
  for (component in seq_len(k)) {
    tau[groups == component, component] <- 0.95
  }
  tau / rowSums(tau)
}

.initial_params_from_tau <- function(y, trials, family, z_design, x_design, du, offset, weights, posterior, control) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  k <- ncol(posterior)
  previous <- list(
    expert_coef = matrix(0, nrow = k, ncol = 2L * q),
    gating_coef = matrix(0, nrow = k, ncol = 2L * p),
    sigma = if (identical(family, "gaussian")) rep(max(stats::sd(y), control$min_sigma), k) else NULL,
    sigma_slope = if (identical(family, "gaussian")) rep(0, k) else NULL,
    theta = if (identical(family, "negative-binomial")) rep(control$negbin_init_theta, k) else NULL
  )
  .m_step(y, trials, family, z_design, x_design, du, offset, weights, posterior, previous, control)
}

.initial_binomial_structured_params <- function(y, trials, z_design, x_design, du, weights, k, control, start_id) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  base_fit <- .fit_binomial_logistic(
    y,
    trials,
    expert_design,
    weights,
    rep(0, 2L * q),
    control
  )
  base_coef <- base_fit$coef
  centers <- seq(-(k - 1) / 2, (k - 1) / 2, length.out = k)
  direction <- if (start_id %% 2L == 0L) 1 else -1
  scale <- control$binomial_start_separation * (1 + 0.20 * max(start_id - 1L, 0L))
  expert_coef <- matrix(base_coef, nrow = k, ncol = length(base_coef), byrow = TRUE)
  jitter <- control$binomial_start_jitter
  for (component in seq_len(k)) {
    component_shift <- direction * centers[component] * scale
    expert_coef[component, 1L] <- expert_coef[component, 1L] + component_shift
    if (q >= 2L) {
      expert_coef[component, 2L] <- expert_coef[component, 2L] - 0.50 * component_shift
    }
    if (2L * q >= q + 1L) {
      expert_coef[component, q + 1L] <- expert_coef[component, q + 1L] + 0.20 * component_shift
    }
    if (jitter > 0 && start_id > 1L) {
      expert_coef[component, ] <- expert_coef[component, ] +
        stats::rnorm(ncol(expert_coef), sd = jitter / sqrt(start_id))
    }
  }

  gating_coef <- matrix(0, nrow = k, ncol = 2L * p)
  for (component in seq_len(k)) {
    gating_coef[component, 1L] <- 0.20 * direction * centers[component]
    if (p >= 2L) {
      gating_coef[component, 2L] <- -0.10 * direction * centers[component]
    }
  }
  gating_coef <- .center_logits(gating_coef)

  list(
    expert_coef = expert_coef,
    gating_coef = gating_coef,
    sigma = NULL,
    sigma_slope = NULL,
    expert_convergence = rep(NA_integer_, k),
    expert_gradient_norm = rep(NA_real_, k),
    expert_coef_norm = apply(expert_coef, 1L, function(x) sqrt(sum(x^2)))
  )
}

.initial_negbin_structured_params <- function(y, z_design, x_design, du, offset, weights, k, control, start_id) {
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert_design <- .make_expert_design(z_design, du)
  base_y <- log(pmax(y, 0) + 0.5) - offset
  base_coef <- .fit_wls(base_y, expert_design, weights, control$negbin_ridge %||% control$ridge)
  centers <- seq(-(k - 1) / 2, (k - 1) / 2, length.out = k)
  direction <- if (start_id %% 2L == 0L) 1 else -1
  scale <- control$negbin_start_separation * (1 + 0.15 * max(start_id - 1L, 0L))
  expert_coef <- matrix(base_coef, nrow = k, ncol = length(base_coef), byrow = TRUE)
  jitter <- control$negbin_start_jitter
  for (component in seq_len(k)) {
    component_shift <- direction * centers[component] * scale
    expert_coef[component, 1L] <- expert_coef[component, 1L] + component_shift
    if (q >= 2L) {
      expert_coef[component, 2L] <- expert_coef[component, 2L] - 0.25 * component_shift
    }
    if (jitter > 0 && start_id > 1L) {
      expert_coef[component, ] <- expert_coef[component, ] +
        stats::rnorm(ncol(expert_coef), sd = jitter / sqrt(start_id))
    }
  }

  gating_coef <- matrix(0, nrow = k, ncol = 2L * p)
  for (component in seq_len(k)) {
    gating_coef[component, 1L] <- 0.15 * direction * centers[component]
    if (p >= 2L) {
      gating_coef[component, 2L] <- -0.08 * direction * centers[component]
    }
  }
  gating_coef <- .center_logits(gating_coef)

  list(
    expert_coef = expert_coef,
    gating_coef = gating_coef,
    sigma = NULL,
    sigma_slope = NULL,
    theta = rep(control$negbin_init_theta, k),
    expert_convergence = rep(NA_integer_, k),
    expert_gradient_norm = rep(NA_real_, k),
    expert_coef_norm = apply(expert_coef, 1L, function(x) sqrt(sum(x^2)))
  )
}

.local_em <- function(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                      estimation_spec, init_local = NULL, init_tau = NULL, init_params = NULL) {
  weights <- .local_weights(u, u0, bandwidth, estimation_spec)
  du <- .local_du(u, u0, bandwidth, estimation_spec)
  if (!is.null(init_params)) {
    params <- init_params
  } else if (!is.null(init_local)) {
    params <- .local_to_params(init_local)
  } else {
    posterior0 <- init_tau %||% .initial_tau(initial_y, k)
    params <- .initial_params_from_tau(y, trials, family, z_design, x_design, du, offset, weights, posterior0, control)
  }

  loglik <- -Inf
  converged <- FALSE
  posterior <- NULL
  iter <- 0L
  for (iter in seq_len(control$maxit)) {
    posterior <- .e_step(y, trials, family, z_design, x_design, du, offset, params, control)
    new_params <- .m_step(y, trials, family, z_design, x_design, du, offset, weights, posterior, params, control)
    new_loglik <- .local_loglik(y, trials, family, z_design, x_design, du, offset, weights, new_params, control)
    if (is.finite(loglik) && abs(new_loglik - loglik) <= control$tol * (abs(loglik) + 1)) {
      params <- new_params
      loglik <- new_loglik
      converged <- TRUE
      break
    }
    params <- new_params
    loglik <- new_loglik
  }

  posterior <- .e_step(y, trials, family, z_design, x_design, du, offset, params, control)
  .params_to_local(params, posterior, z_design, x_design, du, loglik, iter, converged)
}

.fit_local_best <- function(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                            previous, estimation_spec) {
  starts <- max(1L, as.integer(control$n_starts))
  candidates <- vector("list", starts)
  loglik <- rep(-Inf, starts)

  for (start_id in seq_len(starts)) {
    init_local <- NULL
    init_tau <- NULL
    init_params <- NULL
    if (start_id == 1L && !is.null(previous)) {
      init_local <- previous
    } else if (identical(family, "binomial") && isTRUE(control$binomial_structured_starts)) {
      weights <- .local_weights(u, u0, bandwidth, estimation_spec)
      du <- .local_du(u, u0, bandwidth, estimation_spec)
      init_params <- .initial_binomial_structured_params(
        y,
        trials,
        z_design,
        x_design,
        du,
        weights,
        k,
        control,
        start_id
      )
    } else if (identical(family, "negative-binomial") && isTRUE(control$negbin_structured_starts)) {
      weights <- .local_weights(u, u0, bandwidth, estimation_spec)
      du <- .local_du(u, u0, bandwidth, estimation_spec)
      init_params <- .initial_negbin_structured_params(
        y,
        z_design,
        x_design,
        du,
        offset,
        weights,
        k,
        control,
        start_id
      )
    } else {
      init_tau <- .initial_tau(initial_y, k, random = start_id > 1L)
    }
    fit <- tryCatch(
      .local_em(y, trials, family, initial_y, z_design, x_design, offset, u, u0, bandwidth, k, control,
                estimation_spec = estimation_spec,
                init_local = init_local, init_tau = init_tau, init_params = init_params),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      candidates[[start_id]] <- fit
      loglik[start_id] <- fit$loglik
    }
  }

  if (!any(is.finite(loglik))) {
    stop("All EM starts failed at u = ", signif(u0, 4), ".", call. = FALSE)
  }
  best <- which.max(loglik)
  candidates[[best]]$start_loglik <- loglik
  candidates[[best]]$selected_start <- best
  candidates[[best]]
}

.nearest_grid_assignment <- function(u, u_grid) {
  u <- as.numeric(u)
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (length(u_grid) == 1L) {
    return(rep(1L, length(u)))
  }
  mids <- (u_grid[-1L] + u_grid[-length(u_grid)]) / 2
  findInterval(u, mids) + 1L
}

.joint_path_cost_guard <- function(n, n_grid, control) {
  n <- as.integer(n)
  n_grid <- as.integer(n_grid)
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("`n` must be a positive finite observation count.", call. = FALSE)
  }
  if (length(n_grid) != 1L || !is.finite(n_grid) || n_grid < 1L) {
    stop("`n_grid` must be a positive finite grid-point count.", call. = FALSE)
  }
  dense_limit <- suppressWarnings(as.integer(
    control$joint_path_dense_grid_limit %||% 100L
  ))
  if (length(dense_limit) != 1L || !is.finite(dense_limit) || dense_limit < 1L) {
    dense_limit <- 100L
  }
  dense_ratio <- suppressWarnings(as.numeric(
    control$joint_path_dense_grid_ratio %||% 0.2
  ))
  if (length(dense_ratio) != 1L || !is.finite(dense_ratio) || dense_ratio <= 0) {
    dense_ratio <- 0.2
  }

  too_dense <- n_grid > dense_limit || n_grid / max(n, 1L) > dense_ratio
  if (too_dense && !isTRUE(control$allow_dense_u_grid)) {
    stop(
      sprintf(
        paste0(
          "`u_grid` has %d point(s) for %d observation(s), exceeding the ",
          "joint-path dense-grid guard. Set `control$allow_dense_u_grid = TRUE` ",
          "only when this computational cost is intended."
        ),
        n_grid,
        n
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.joint_path_unique_grid_warning <- function(u_values, u_grid) {
  u_grid <- sort(unique(as.numeric(u_grid)))
  unique_u <- sort(unique(as.numeric(u_values)))
  if (length(unique_u) == length(u_grid) && length(u_grid) > 20L &&
      max(abs(unique_u - sort(u_grid))) <= sqrt(.Machine$double.eps)) {
    warning(
      "`u_grid` is effectively `unique(u)`; this fits almost one local model per observation.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.joint_path_dense_grid_guard <- function(u_values, u_grid, control) {
  .joint_path_cost_guard(length(u_values), length(u_grid), control)
  .joint_path_unique_grid_warning(u_values, u_grid)
  invisible(TRUE)
}

.joint_path_param_vector <- function(params) {
  as.numeric(unlist(
    list(
      params$expert_coef,
      params$gating_coef,
      params$sigma,
      params$sigma_slope,
      params$theta
    ),
    use.names = FALSE
  ))
}

.joint_path_param_delta <- function(new_params, old_params) {
  if (is.null(old_params)) {
    return(NA_real_)
  }
  new_vec <- .joint_path_param_vector(new_params)
  old_vec <- .joint_path_param_vector(old_params)
  if (length(new_vec) != length(old_vec)) {
    return(NA_real_)
  }
  delta <- abs(new_vec - old_vec)
  if (!any(is.finite(delta))) {
    return(NA_real_)
  }
  max(delta[is.finite(delta)])
}

.joint_path_sample_loglik <- function(y, trials, family, z_design, x_design,
                                      offset, u, u_grid, bandwidth, assignment,
                                      params_list, control, estimation_spec) {
  total <- 0
  for (grid_id in seq_along(u_grid)) {
    rows <- assignment == grid_id
    if (!any(rows)) {
      next
    }
    du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    total <- total + .local_loglik(
      y,
      trials,
      family,
      z_design,
      x_design,
      du,
      offset,
      as.numeric(rows),
      params_list[[grid_id]],
      control
    )
  }
  as.numeric(total)
}

.joint_path_em_one_start <- function(y, trials, family, initial_y, z_design, x_design,
                                     offset, u, u_grid, bandwidth, k, control,
                                     estimation_spec, assignment, progress = NULL,
                                     start_id = 1L, n_starts = 1L) {
  n <- length(y)
  n_grid <- length(u_grid)
  path_posterior <- .initial_tau(initial_y, k, random = start_id > 1L)
  params_list <- vector("list", n_grid)
  loglik <- rep(NA_real_, n_grid)
  trace <- data.frame(
    iteration = integer(0L),
    objective = numeric(0L),
    objective_delta = numeric(0L),
    posterior_delta = numeric(0L),
    parameter_delta = numeric(0L),
    converged = logical(0L)
  )
  objective <- -Inf
  converged <- FALSE
  posterior_by_grid <- array(NA_real_, dim = c(n, k, n_grid))
  start_start_time <- proc.time()[["elapsed"]]

  for (iter in seq_len(control$maxit)) {
    previous_posterior <- path_posterior
    previous_params <- params_list
    max_param_delta <- NA_real_
    objective_old <- objective

    for (grid_id in seq_along(u_grid)) {
      weights <- .local_weights(u, u_grid[[grid_id]], bandwidth, estimation_spec)
      du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
      previous <- params_list[[grid_id]]
      if (is.null(previous)) {
        previous <- .initial_params_from_tau(
          y,
          trials,
          family,
          z_design,
          x_design,
          du,
          offset,
          weights,
          path_posterior,
          control
        )
      }
      params <- .m_step(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        path_posterior,
        previous,
        control
      )
      params_list[[grid_id]] <- params
      posterior_by_grid[, , grid_id] <- .e_step(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        params,
        control
      )
      loglik[[grid_id]] <- .local_loglik(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        params,
        control
      )
      grid_delta <- .joint_path_param_delta(params, previous_params[[grid_id]])
      if (is.finite(grid_delta)) {
        max_param_delta <- if (is.finite(max_param_delta)) {
          max(max_param_delta, grid_delta)
        } else {
          grid_delta
        }
      }
    }

    for (grid_id in seq_along(u_grid)) {
      rows <- assignment == grid_id
      if (any(rows)) {
        path_posterior[rows, ] <- matrix(
          posterior_by_grid[rows, , grid_id],
          ncol = k,
          byrow = FALSE
        )
      }
    }
    path_posterior <- pmax(path_posterior, 1e-12)
    path_posterior <- path_posterior / rowSums(path_posterior)

    objective <- .joint_path_sample_loglik(
      y,
      trials,
      family,
      z_design,
      x_design,
      offset,
      u,
      u_grid,
      bandwidth,
      assignment,
      params_list,
      control,
      estimation_spec
    )
    objective_delta <- if (is.finite(objective_old)) objective - objective_old else NA_real_
    posterior_delta <- max(abs(path_posterior - previous_posterior), na.rm = TRUE)
    did_converge <- is.finite(posterior_delta) &&
      posterior_delta <= control$tol &&
      (!is.finite(max_param_delta) || max_param_delta <= sqrt(control$tol))
    trace <- rbind(
      trace,
      data.frame(
        iteration = iter,
        objective = objective,
        objective_delta = objective_delta,
        posterior_delta = posterior_delta,
        parameter_delta = max_param_delta,
        converged = did_converge
      )
    )
    if (.vcmoe_progress_due(control, iter, did_converge)) {
      .vcmoe_progress_log(
        progress,
        event = "joint_path_iteration",
        status = if (did_converge) "converged" else "running",
        family = family,
        k = k,
        n = n,
        bandwidth = bandwidth,
        n_grid = n_grid,
        start_id = start_id,
        n_starts = n_starts,
        iteration = iter,
        maxit = control$maxit,
        start_start_time = start_start_time,
        loglik = objective,
        loglik_delta = objective_delta,
        converged = did_converge
      )
    }
    if (did_converge) {
      converged <- TRUE
      break
    }
  }

  iterations <- if (nrow(trace)) utils::tail(trace$iteration, 1L) else 0L
  locals <- vector("list", n_grid)
  for (grid_id in seq_along(u_grid)) {
    du <- .local_du(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    weights <- .local_weights(u, u_grid[[grid_id]], bandwidth, estimation_spec)
    posterior <- .e_step(
      y,
      trials,
      family,
      z_design,
      x_design,
      du,
      offset,
      params_list[[grid_id]],
      control
    )
    locals[[grid_id]] <- .params_to_local(
      params_list[[grid_id]],
      posterior,
      z_design,
      x_design,
      du,
      .local_loglik(
        y,
        trials,
        family,
        z_design,
        x_design,
        du,
        offset,
        weights,
        params_list[[grid_id]],
        control
      ),
      iterations,
      converged
    )
    locals[[grid_id]]$selected_start <- start_id
    locals[[grid_id]]$start_loglik <- objective
  }

  list(
    locals = locals,
    path_posterior = path_posterior,
    assignment = assignment,
    loglik = objective,
    loglik_by_grid = loglik,
    trace = trace,
    converged = converged,
    iterations = iterations,
    selected_start = start_id
  )
}

.joint_path_fit_best <- function(y, trials, family, initial_y, z_design, x_design,
                                 offset, u, u_grid, bandwidth, k, control,
                                 estimation_spec, assignment, progress = NULL) {
  starts <- max(1L, as.integer(control$n_starts))
  candidates <- vector("list", starts)
  loglik <- rep(-Inf, starts)
  errors <- rep("", starts)

  for (start_id in seq_len(starts)) {
    start_start_time <- proc.time()[["elapsed"]]
    error_message <- NULL
    fit <- tryCatch(
      .joint_path_em_one_start(
        y,
        trials,
        family,
        initial_y,
        z_design,
        x_design,
        offset,
        u,
        u_grid,
        bandwidth,
        k,
        control,
        estimation_spec,
        assignment,
        progress = progress,
        start_id = start_id,
        n_starts = starts
      ),
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    if (!is.null(fit) && is.finite(fit$loglik)) {
      candidates[[start_id]] <- fit
      loglik[[start_id]] <- fit$loglik
      .vcmoe_progress_log(
        progress,
        event = "joint_path_start_finished",
        status = "completed",
        family = family,
        k = k,
        n = length(y),
        bandwidth = bandwidth,
        n_grid = length(u_grid),
        start_id = start_id,
        n_starts = starts,
        maxit = control$maxit,
        start_start_time = start_start_time,
        loglik = fit$loglik,
        converged = fit$converged
      )
    } else {
      errors[[start_id]] <- error_message %||% "non-finite joint-path objective"
      .vcmoe_progress_log(
        progress,
        event = "joint_path_start_failed",
        status = "failed",
        family = family,
        k = k,
        n = length(y),
        bandwidth = bandwidth,
        n_grid = length(u_grid),
        start_id = start_id,
        n_starts = starts,
        maxit = control$maxit,
        start_start_time = start_start_time,
        error_message = errors[[start_id]]
      )
    }
  }

  if (!any(is.finite(loglik))) {
    detail <- unique(errors[nzchar(errors)])
    suffix <- if (length(detail)) paste0(" ", paste(detail, collapse = "; ")) else ""
    stop("All joint-path EM starts failed.", suffix, call. = FALSE)
  }
  selected <- which.max(loglik)
  out <- candidates[[selected]]
  out$selected_start <- selected
  out$start_loglik <- loglik
  out
}

.vcmoe_fit_joint_path <- function(formula, data, u, k, family, bandwidth, u_grid,
                                  control, label, parameterization, u_scale,
                                  progress, fit_call) {
  if (is.list(formula)) {
    stop("Joint-path EM supports single-response formulas only.", call. = FALSE)
  }
  family <- match.arg(family, c("gaussian", "binomial", "negative-binomial"))
  label <- match.arg(label, c("align", "global", "greedy"))
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  estimation_spec <- .vcmoe_estimation_spec(parameterization)
  k <- as.integer(k)
  if (k < 2L || k > 10L) {
    stop("`k` must be between 2 and 10.", call. = FALSE)
  }
  alignment_method <- if (identical(label, "greedy")) {
    "greedy"
  } else if (identical(label, "global") && k <= 6L) {
    "global"
  } else if (k <= 6L) {
    "global"
  } else {
    "sequential"
  }

  user_control <- control %||% list()
  ridge_user_supplied <- "ridge" %in% names(user_control)
  control <- .vcmoe_joint_path_control(user_control)
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }

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
  expert_offset_all <- .extract_model_offset(expert_frame_all)
  gating_offset_all <- .extract_model_offset(gating_frame_all)
  if (any(abs(gating_offset_all) > 0)) {
    stop("Gating-side `offset()` terms are not supported in v0.", call. = FALSE)
  }
  if (!identical(family, "negative-binomial") && any(abs(expert_offset_all) > 0)) {
    stop("Expert-side `offset()` terms are currently supported only for `family = \"negative-binomial\"`.", call. = FALSE)
  }
  complete_rows <- .complete_cases_frame(expert_frame_all, nrow(data)) &
    .complete_cases_frame(gating_frame_all, nrow(data)) &
    is.finite(u_values_all) &
    is.finite(expert_offset_all)
  if (!any(complete_rows)) {
    stop("No complete rows remain after checking response, covariates, and `u`.", call. = FALSE)
  }
  if (!all(complete_rows)) {
    warning(
      sprintf("Removed %d row(s) with missing or non-finite response, covariates, or `u`.",
              sum(!complete_rows)),
      call. = FALSE
    )
  }

  expert_frame <- expert_frame_all[complete_rows, , drop = FALSE]
  gating_frame <- gating_frame_all[complete_rows, , drop = FALSE]
  response <- .parse_vcmoe_response(family, stats::model.response(expert_frame), pieces$response)
  if (identical(family, "binomial") &&
      identical(response$type, "bernoulli") &&
      !ridge_user_supplied) {
    control$ridge <- 1
  }
  y <- response$y
  trials <- response$trials
  initial_y <- response$observed
  z_design <- stats::model.matrix(stats::terms(pieces$expert_formula), expert_frame)
  x_design <- stats::model.matrix(stats::terms(pieces$gating_formula), data = gating_frame)
  expert_offset <- expert_offset_all[complete_rows]
  u_original <- u_values_all[complete_rows]
  u_scaling <- .vcmoe_make_u_scaling(u_original, u_scale)
  u_values <- .vcmoe_scale_u(u_original, u_scaling)

  if (is.null(bandwidth)) {
    bandwidth <- .default_bandwidth(u_values)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (is.null(u_grid)) {
    u_grid <- .default_u_grid(u_values)
  }
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (!length(u_grid)) {
    stop("`u_grid` must contain at least one value.", call. = FALSE)
  }
  if (any(!is.finite(u_grid))) {
    stop("`u_grid` must contain finite values on the analysis scale.", call. = FALSE)
  }
  if (identical(u_scale, "unit") && (min(u_grid) < -1e-8 || max(u_grid) > 1 + 1e-8)) {
    warning(
      "`u_scale = \"unit\"` interprets `u_grid` on the scaled [0, 1] analysis scale.",
      call. = FALSE
    )
  }
  .joint_path_dense_grid_guard(u_values, u_grid, control)
  u_grid_original <- .vcmoe_unscale_u(u_grid, u_scaling)

  n_grid <- length(u_grid)
  q <- ncol(z_design)
  p <- ncol(x_design)
  assignment <- .nearest_grid_assignment(u_values, u_grid)
  assignment_counts <- tabulate(assignment, nbins = n_grid)
  progress_state <- .vcmoe_progress_setup(progress)
  warning(
    paste0(
      "Joint-path EM fits one local model per grid point in each EM iteration; ",
      "runtime scales approximately with `control$n_starts * control$maxit * ",
      "length(u_grid) * n`."
    ),
    call. = FALSE
  )
  .vcmoe_progress_log(
    progress_state,
    event = "fit_started",
    status = "running",
    family = family,
    k = k,
    n = length(y),
    bandwidth = bandwidth,
    n_grid = n_grid,
    n_starts = control$n_starts,
    maxit = control$maxit
  )
  joint <- .joint_path_fit_best(
    y,
    trials,
    family,
    initial_y,
    z_design,
    x_design,
    expert_offset,
    u_values,
    u_grid,
    bandwidth,
    k,
    control,
    estimation_spec,
    assignment,
    progress = progress_state
  )
  if (!isTRUE(joint$converged)) {
    warning(
      sprintf(
        paste0(
          "Joint-path EM did not converge within %d outer iteration(s); ",
          "inspect `fit$diagnostics$joint_path_trace` before interpretation."
        ),
        control$maxit
      ),
      call. = FALSE
    )
  }
  locals <- joint$locals

  previous <- NULL
  for (grid_id in seq_along(u_grid)) {
    locals[[grid_id]] <- .align_local_fit(
      locals[[grid_id]],
      previous,
      z_design,
      x_design,
      initial_y,
      control
    )
    previous <- locals[[grid_id]]
  }
  diagnostics <- list(
    engine = .vcmoe_joint_path_engine_id,
    alignment_method = alignment_method,
    loglik = numeric(n_grid),
    converged = logical(n_grid),
    iterations = integer(n_grid),
    selected_start = integer(n_grid),
    permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    greedy_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    global_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_optimizer_convergence = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_gradient_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    expert_coef_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    global_path_cost = NA_real_,
    global_transition_cost = rep(NA_real_, n_grid),
    alignment_margin = numeric(n_grid),
    ambiguous = logical(n_grid),
    posterior_entropy = numeric(n_grid),
    warnings = character(0L),
    start_loglik = vector("list", n_grid),
    joint_path_iterations = joint$iterations,
    joint_path_converged = joint$converged,
    joint_path_objective = "nearest_grid_sample_loglik",
    joint_path_paper_criterion_match = "nearest_grid_approximation",
    joint_path_trace = joint$trace,
    joint_path_assignment = data.frame(
      grid_id = seq_len(n_grid),
      u = u_grid,
      u_original = u_grid_original,
      n_assigned = assignment_counts
    ),
    joint_path_selected_start = joint$selected_start,
    joint_path_start_loglik = joint$start_loglik
  )

  if (identical(alignment_method, "global")) {
    global_alignment <- .global_align_local_fits(
      locals,
      u_grid,
      z_design,
      x_design,
      initial_y,
      control,
      bandwidth = bandwidth,
      estimation_spec = estimation_spec
    )
    locals <- global_alignment$locals
    diagnostics$global_path_cost <- global_alignment$path_cost
    diagnostics$global_transition_cost <- global_alignment$transition_cost
  }

  expert <- array(NA_real_, dim = c(n_grid, k, q),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(z_design)))
  expert_slope <- expert
  gating <- array(NA_real_, dim = c(n_grid, k, p),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(x_design)))
  gating_slope <- gating
  sigma <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  sigma_slope <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  theta <- if (identical(family, "negative-binomial")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  posterior <- array(NA_real_, dim = c(length(y), k, n_grid),
                     dimnames = list(observation = NULL, component = paste0("component", seq_len(k)),
                                     u = signif(u_grid, 8)))
  path_posterior <- matrix(
    NA_real_,
    nrow = length(y),
    ncol = k,
    dimnames = list(observation = NULL, component = paste0("component", seq_len(k)))
  )

  for (grid_id in seq_along(u_grid)) {
    local <- locals[[grid_id]]
    diagnostics$ambiguous[grid_id] <- isTRUE(local$ambiguous)
    if (isTRUE(local$ambiguous)) {
      diagnostics$warnings <- c(
        diagnostics$warnings,
        sprintf("Ambiguous label alignment near u = %.4f.", u_grid[grid_id])
      )
    }

    expert[grid_id, , ] <- local$expert_intercept
    expert_slope[grid_id, , ] <- local$expert_slope
    gating[grid_id, , ] <- local$gating_intercept
    gating_slope[grid_id, , ] <- local$gating_slope
    if (!is.null(sigma)) {
      sigma[grid_id, ] <- local$sigma
    }
    if (!is.null(sigma_slope)) {
      sigma_slope[grid_id, ] <- local$sigma_slope
    }
    if (!is.null(theta)) {
      theta[grid_id, ] <- local$theta
    }
    posterior[, , grid_id] <- local$posterior
    rows <- assignment == grid_id
    if (any(rows)) {
      path_posterior[rows, ] <- matrix(
        local$posterior[rows, ],
        ncol = k,
        byrow = FALSE
      )
    }

    diagnostics$loglik[grid_id] <- local$loglik
    diagnostics$converged[grid_id] <- local$converged
    diagnostics$iterations[grid_id] <- local$iterations
    diagnostics$selected_start[grid_id] <- joint$selected_start
    diagnostics$permutations[grid_id, ] <- local$permutation
    diagnostics$greedy_permutations[grid_id, ] <- local$greedy_permutation %||% local$permutation
    diagnostics$global_permutations[grid_id, ] <- local$global_permutation %||% seq_len(k)
    diagnostics$expert_optimizer_convergence[grid_id, ] <- local$expert_convergence %||% rep(NA_integer_, k)
    diagnostics$expert_gradient_norm[grid_id, ] <- local$expert_gradient_norm %||% rep(NA_real_, k)
    diagnostics$expert_coef_norm[grid_id, ] <- local$expert_coef_norm %||% rep(NA_real_, k)
    diagnostics$alignment_margin[grid_id] <- local$alignment_margin
    diagnostics$posterior_entropy[grid_id] <- .weighted_entropy(
      local$posterior,
      .local_weights(u_values, u_grid[grid_id], bandwidth, estimation_spec)
    )
    diagnostics$start_loglik[[grid_id]] <- joint$start_loglik
  }

  if (length(diagnostics$warnings) && isTRUE(control$warn_ambiguous)) {
    warning(
      sprintf("Label alignment was ambiguous at %d grid point(s); inspect `fit$diagnostics$warnings`.",
              length(diagnostics$warnings)),
      call. = FALSE
    )
  }
  .vcmoe_progress_log(
    progress_state,
    event = "fit_finished",
    status = "completed",
    family = family,
    k = k,
    n = length(y),
    bandwidth = bandwidth,
    n_grid = n_grid,
    n_starts = control$n_starts,
    maxit = control$maxit,
    selected_start = joint$selected_start,
    converged = isTRUE(joint$converged)
  )

  fit_call$progress <- NULL
  parameterization_meta <- .vcmoe_parameterization_metadata(
    family = family,
    k = k,
    bandwidth = bandwidth,
    label = label,
    alignment_method = alignment_method,
    control = control,
    u_scaling = u_scaling,
    parameterization = parameterization,
    estimation_spec = estimation_spec
  )
  parameterization_meta$estimator <- .vcmoe_joint_path_engine_id
  parameterization_meta$engine <- list(
    id = .vcmoe_joint_path_engine_id,
    responsibility_update = "nearest_grid_observation_path",
    diagnostic_objective = "nearest_grid_sample_loglik",
    paper_criterion_match = "nearest_grid_approximation",
    runtime = "Fits one local model per grid point in each EM iteration."
  )

  object <- list(
    call = fit_call,
    formula = formula,
    family = family,
    k = k,
    label = alignment_method,
    engine_id = .vcmoe_joint_path_engine_id,
    parameterization_id = parameterization,
    u_scale = u_scale,
    u_scaling = u_scaling,
    bandwidth = bandwidth,
    u_grid = u_grid,
    u_grid_original = u_grid_original,
    control = control,
    parameterization = parameterization_meta,
    coefficients = list(
      expert = expert,
      expert_slope = expert_slope,
      gating = gating,
      gating_slope = gating_slope,
      sigma = sigma,
      sigma_slope = sigma_slope,
      theta = theta
    ),
    diagnostics = diagnostics,
    posterior = posterior,
    joint_path_posterior = path_posterior,
    terms = list(
      expert = stats::terms(pieces$expert_formula),
      gating = stats::terms(pieces$gating_formula)
    ),
    response = pieces$response,
    response_info = response$info,
    u_name = u_info$name,
    rows_used = which(complete_rows),
    rows_removed = sum(!complete_rows),
    fitted = if (isTRUE(control$keep_data)) {
      list(
        y = response$observed,
        success = response$success,
        trials = response$trials,
        response = response,
        z_design = z_design,
        x_design = x_design,
        expert_offset = expert_offset,
        u = u_values,
        u_original = u_original
      )
    } else {
      NULL
    }
  )
  class(object) <- "vcmoe"
  object
}

#' Fit a varying-coefficient mixture-of-experts model
#'
#' @param formula A formula of the form `y ~ expert_terms | gating_terms`.
#'   For grouped binomial data, use `cbind(success, failure) ~ expert_terms |
#'   gating_terms`. For Negative-Binomial count data, expert-side
#'   `offset(log_size_factor)` terms are supported.
#' @param data A data frame.
#' @param u Continuous index column name or numeric vector.
#' @param k Number of mixture components. `k = 2` is the primary v0 target;
#'   `k = 3:10` are high-k candidate support and require diagnostics before
#'   being treated as stable.
#' @param family Model family. `"gaussian"`, `"binomial"`, and
#'   `"negative-binomial"` are implemented.
#' @param bandwidth Kernel bandwidth. If `NULL`, a Silverman-style default is
#'   used for `u`.
#' @param u_grid Grid where coefficient functions are estimated.
#' @param control Named list overriding EM and label-alignment settings.
#' @param label Label strategy. `"align"` uses exact global alignment for
#'   `k <= 6` and sequential assignment with ambiguity margins for `k >= 7`;
#'   `"global"` requests exact global alignment when feasible and falls back
#'   to the same sequential assignment path for `k >= 7`;
#'   `"greedy"` keeps the older one-step alignment.
#' @param u_scale How to transform `u` before fitting. The default `"unit"`
#'   maps complete-row `u` values to `[0, 1]`; `"none"` leaves `u` on the
#'   supplied scale. `bandwidth` and `u_grid` are interpreted on the transformed
#'   analysis scale.
#' @param parameterization Estimator convention. The public package uses
#'   `"a1_epanechnikov_scaled"`: Epanechnikov density weights `K((u-u0)/h)/h`
#'   and the scaled local-linear basis `(u-u0)/h`.
#' @param engine Fitting engine. The default `"local_grid_em"` preserves the
#'   original independent local-grid EM path. `"joint_path_em"` updates a
#'   shared observation-level responsibility path across grid points.
#' @param progress Joint-path progress reporting. The default `NULL` is silent;
#'   use `TRUE` for messages or a single CSV file path for structured logging.
#'   Iteration frequency is controlled by `control$progress_every`.
#' @return An object of class `vcmoe`.
#' @details Rows with missing or non-finite response, covariates, or `u` are
#'   removed consistently before fitting, with a warning. For single-trial
#'   Bernoulli responses, the default gating ridge is strengthened to
#'   `control$ridge = 1` unless the user explicitly supplies `control$ridge`;
#'   grouped Binomial and other families keep the global default. Joint-path
#'   traces record the sample-level nearest-grid log-likelihood as a diagnostic
#'   criterion. The label-consistent updates do not guarantee that this
#'   diagnostic is monotone; convergence is based on posterior and parameter
#'   deltas instead.
#' @export
vcmoe_fit <- function(formula, data, u, k = 2L, family = "gaussian",
                      bandwidth = NULL, u_grid = NULL, control = list(),
                      label = "align",
                      parameterization = "a1_epanechnikov_scaled",
                      u_scale = c("unit", "none"),
                      engine = c("local_grid_em", "joint_path_em"),
                      progress = NULL) {
  engine <- match.arg(engine)
  if (identical(engine, .vcmoe_joint_path_engine_id)) {
    return(.vcmoe_fit_joint_path(
      formula = formula,
      data = data,
      u = u,
      k = k,
      family = family,
      bandwidth = bandwidth,
      u_grid = u_grid,
      control = control,
      label = label,
      parameterization = parameterization,
      u_scale = u_scale,
      progress = progress,
      fit_call = match.call()
    ))
  }
  if (!is.null(progress) && !identical(progress, FALSE)) {
    stop("`progress` is available only for `engine = \"joint_path_em\"`.", call. = FALSE)
  }
  family <- match.arg(family, c("gaussian", "binomial", "negative-binomial"))
  label <- match.arg(label, c("align", "global", "greedy"))
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  estimation_spec <- .vcmoe_estimation_spec(parameterization)
  k <- as.integer(k)
  if (k < 2L || k > 10L) {
    stop("`k` must be between 2 and 10.", call. = FALSE)
  }
  alignment_method <- if (identical(label, "greedy")) {
    "greedy"
  } else if (identical(label, "global") && k <= 6L) {
    "global"
  } else if (k <= 6L) {
    "global"
  } else {
    "sequential"
  }
  if (identical(family, "binomial")) {
    if (k >= 3L) {
      warning(
        "Binomial VCMoE with k >= 3 is experimental high-k candidate support; inspect diagnostics before interpretation.",
        call. = FALSE
      )
    }
  } else if (identical(family, "negative-binomial") && k >= 3L) {
    warning(
      "Negative-Binomial VCMoE with k >= 3 is high-k candidate support; stable status depends on validation diagnostics.",
      call. = FALSE
    )
  }

  user_control <- control %||% list()
  ridge_user_supplied <- "ridge" %in% names(user_control)
  control <- .vcmoe_default_control(user_control)
  if (!is.null(control$seed)) {
    set.seed(control$seed)
  }

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
  expert_offset_all <- .extract_model_offset(expert_frame_all)
  gating_offset_all <- .extract_model_offset(gating_frame_all)
  if (any(abs(gating_offset_all) > 0)) {
    stop("Gating-side `offset()` terms are not supported in v0.", call. = FALSE)
  }
  if (!identical(family, "negative-binomial") && any(abs(expert_offset_all) > 0)) {
    stop("Expert-side `offset()` terms are currently supported only for `family = \"negative-binomial\"`.", call. = FALSE)
  }
  complete_rows <- .complete_cases_frame(expert_frame_all, nrow(data)) &
    .complete_cases_frame(gating_frame_all, nrow(data)) &
    is.finite(u_values_all) &
    is.finite(expert_offset_all)
  if (!any(complete_rows)) {
    stop("No complete rows remain after checking response, covariates, and `u`.", call. = FALSE)
  }
  if (!all(complete_rows)) {
    warning(
      sprintf("Removed %d row(s) with missing or non-finite response, covariates, or `u`.",
              sum(!complete_rows)),
      call. = FALSE
    )
  }

  expert_frame <- expert_frame_all[complete_rows, , drop = FALSE]
  gating_frame <- gating_frame_all[complete_rows, , drop = FALSE]
  response <- .parse_vcmoe_response(family, stats::model.response(expert_frame), pieces$response)
  if (identical(family, "binomial") &&
      identical(response$type, "bernoulli") &&
      !ridge_user_supplied) {
    control$ridge <- 1
  }
  y <- response$y
  trials <- response$trials
  initial_y <- response$observed
  z_design <- stats::model.matrix(stats::terms(pieces$expert_formula), expert_frame)
  x_design <- stats::model.matrix(stats::terms(pieces$gating_formula), data = gating_frame)
  expert_offset <- expert_offset_all[complete_rows]
  u_original <- u_values_all[complete_rows]
  u_scaling <- .vcmoe_make_u_scaling(u_original, u_scale)
  u_values <- .vcmoe_scale_u(u_original, u_scaling)

  if (is.null(bandwidth)) {
    bandwidth <- .default_bandwidth(u_values)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0) {
    stop("`bandwidth` must be a positive finite number.", call. = FALSE)
  }
  if (is.null(u_grid)) {
    u_grid <- .default_u_grid(u_values)
  }
  u_grid <- sort(unique(as.numeric(u_grid)))
  if (!length(u_grid)) {
    stop("`u_grid` must contain at least one value.", call. = FALSE)
  }
  if (any(!is.finite(u_grid))) {
    stop("`u_grid` must contain finite values on the analysis scale.", call. = FALSE)
  }
  if (identical(u_scale, "unit") && (min(u_grid) < -1e-8 || max(u_grid) > 1 + 1e-8)) {
    warning(
      "`u_scale = \"unit\"` interprets `u_grid` on the scaled [0, 1] analysis scale.",
      call. = FALSE
    )
  }
  u_grid_original <- .vcmoe_unscale_u(u_grid, u_scaling)

  n_grid <- length(u_grid)
  q <- ncol(z_design)
  p <- ncol(x_design)
  expert <- array(NA_real_, dim = c(n_grid, k, q),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(z_design)))
  expert_slope <- expert
  gating <- array(NA_real_, dim = c(n_grid, k, p),
                  dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k)),
                                  term = colnames(x_design)))
  gating_slope <- gating
  sigma <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  sigma_slope <- if (identical(family, "gaussian")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  theta <- if (identical(family, "negative-binomial")) {
    matrix(NA_real_, nrow = n_grid, ncol = k,
           dimnames = list(u = signif(u_grid, 8), component = paste0("component", seq_len(k))))
  } else {
    NULL
  }
  posterior <- array(NA_real_, dim = c(length(y), k, n_grid),
                     dimnames = list(observation = NULL, component = paste0("component", seq_len(k)),
                                     u = signif(u_grid, 8)))

  diagnostics <- list(
    engine = "local_grid_em",
    alignment_method = alignment_method,
    loglik = numeric(n_grid),
    converged = logical(n_grid),
    iterations = integer(n_grid),
    selected_start = integer(n_grid),
    permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    greedy_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    global_permutations = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_optimizer_convergence = matrix(NA_integer_, nrow = n_grid, ncol = k),
    expert_gradient_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    expert_coef_norm = matrix(NA_real_, nrow = n_grid, ncol = k),
    global_path_cost = NA_real_,
    global_transition_cost = rep(NA_real_, n_grid),
    alignment_margin = numeric(n_grid),
    ambiguous = logical(n_grid),
    posterior_entropy = numeric(n_grid),
    warnings = character(0L),
    start_loglik = vector("list", n_grid)
  )

  previous <- NULL
  locals <- vector("list", n_grid)
  for (grid_id in seq_along(u_grid)) {
    local <- .fit_local_best(y, trials, family, initial_y, z_design, x_design, expert_offset,
                             u_values, u_grid[grid_id],
                             bandwidth, k, control, previous, estimation_spec)
    local <- .align_local_fit(local, previous, z_design, x_design, initial_y, control)
    locals[[grid_id]] <- local
    previous <- local
  }

  if (identical(alignment_method, "global")) {
    global_alignment <- .global_align_local_fits(
      locals,
      u_grid,
      z_design,
      x_design,
      initial_y,
      control,
      bandwidth = bandwidth,
      estimation_spec = estimation_spec
    )
    locals <- global_alignment$locals
    diagnostics$global_path_cost <- global_alignment$path_cost
    diagnostics$global_transition_cost <- global_alignment$transition_cost
  }

  for (grid_id in seq_along(u_grid)) {
    local <- locals[[grid_id]]
    diagnostics$ambiguous[grid_id] <- isTRUE(local$ambiguous)
    if (isTRUE(local$ambiguous)) {
      diagnostics$warnings <- c(
        diagnostics$warnings,
        sprintf("Ambiguous label alignment near u = %.4f.", u_grid[grid_id])
      )
    }

    expert[grid_id, , ] <- local$expert_intercept
    expert_slope[grid_id, , ] <- local$expert_slope
    gating[grid_id, , ] <- local$gating_intercept
    gating_slope[grid_id, , ] <- local$gating_slope
    if (!is.null(sigma)) {
      sigma[grid_id, ] <- local$sigma
    }
    if (!is.null(sigma_slope)) {
      sigma_slope[grid_id, ] <- local$sigma_slope
    }
    if (!is.null(theta)) {
      theta[grid_id, ] <- local$theta
    }
    posterior[, , grid_id] <- local$posterior

    diagnostics$loglik[grid_id] <- local$loglik
    diagnostics$converged[grid_id] <- local$converged
    diagnostics$iterations[grid_id] <- local$iterations
    diagnostics$selected_start[grid_id] <- local$selected_start
    diagnostics$permutations[grid_id, ] <- local$permutation
    diagnostics$greedy_permutations[grid_id, ] <- local$greedy_permutation %||% local$permutation
    diagnostics$global_permutations[grid_id, ] <- local$global_permutation %||% seq_len(k)
    diagnostics$expert_optimizer_convergence[grid_id, ] <- local$expert_convergence %||% rep(NA_integer_, k)
    diagnostics$expert_gradient_norm[grid_id, ] <- local$expert_gradient_norm %||% rep(NA_real_, k)
    diagnostics$expert_coef_norm[grid_id, ] <- local$expert_coef_norm %||% rep(NA_real_, k)
    diagnostics$alignment_margin[grid_id] <- local$alignment_margin
    diagnostics$posterior_entropy[grid_id] <- .weighted_entropy(
      local$posterior,
      .local_weights(u_values, u_grid[grid_id], bandwidth, estimation_spec)
    )
    diagnostics$start_loglik[[grid_id]] <- local$start_loglik
  }

  if (length(diagnostics$warnings) && isTRUE(control$warn_ambiguous)) {
    warning(
      sprintf("Label alignment was ambiguous at %d grid point(s); inspect `fit$diagnostics$warnings`.",
              length(diagnostics$warnings)),
      call. = FALSE
    )
  }

  object <- list(
    call = match.call(),
    formula = formula,
    family = family,
    k = k,
    label = alignment_method,
    engine_id = "local_grid_em",
    parameterization_id = parameterization,
    u_scale = u_scale,
    u_scaling = u_scaling,
    bandwidth = bandwidth,
    u_grid = u_grid,
    u_grid_original = u_grid_original,
    control = control,
    parameterization = .vcmoe_parameterization_metadata(
      family = family,
      k = k,
      bandwidth = bandwidth,
      label = label,
      alignment_method = alignment_method,
      control = control,
      u_scaling = u_scaling,
      parameterization = parameterization,
      estimation_spec = estimation_spec
    ),
    coefficients = list(
      expert = expert,
      expert_slope = expert_slope,
      gating = gating,
      gating_slope = gating_slope,
      sigma = sigma,
      sigma_slope = sigma_slope,
      theta = theta
    ),
    diagnostics = diagnostics,
    posterior = posterior,
    terms = list(
      expert = stats::terms(pieces$expert_formula),
      gating = stats::terms(pieces$gating_formula)
    ),
    response = pieces$response,
    response_info = response$info,
    u_name = u_info$name,
    rows_used = which(complete_rows),
    rows_removed = sum(!complete_rows),
    fitted = if (isTRUE(control$keep_data)) {
      list(
        y = response$observed,
        success = response$success,
        trials = response$trials,
        response = response,
        z_design = z_design,
        x_design = x_design,
        expert_offset = expert_offset,
        u = u_values,
        u_original = u_original
      )
    } else {
      NULL
    }
  )
  class(object) <- "vcmoe"
  object
}

#' @export
print.vcmoe <- function(x, ...) {
  cat("VCMoE fit\n")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  components: ", x$k, "\n", sep = "")
  cat("  label alignment: ", x$label %||% x$diagnostics$alignment_method %||% "unknown", "\n", sep = "")
  parameterization <- vcmoe_parameterization(x)
  cat("  kernel: ", parameterization$kernel$name, " (",
      parameterization$kernel$weight_normalization, ")\n", sep = "")
  cat("  local basis: ", parameterization$local_linear_basis$slope_column, "\n", sep = "")
  cat("  u scale: ", parameterization$u_scaling$method %||% x$u_scale %||% "none", "\n", sep = "")
  cat("  grid points: ", length(x$u_grid), "\n", sep = "")
  cat("  bandwidth: ", signif(x$bandwidth, 4), "\n", sep = "")
  cat("  converged grid points: ", sum(x$diagnostics$converged), "/",
      length(x$u_grid), "\n", sep = "")
  if (length(x$diagnostics$warnings)) {
    cat("  label warnings: ", length(x$diagnostics$warnings), "\n", sep = "")
  }
  invisible(x)
}
