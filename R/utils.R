.vcmoe_default_control <- function(control) {
  defaults <- list(
    maxit = 100,
    tol = 1e-5,
    n_starts = 3,
    seed = NULL,
    ridge = 1e-6,
    binomial_ridge = 1,
    binomial_structured_starts = TRUE,
    binomial_start_separation = 1,
    binomial_start_jitter = 0.10,
    negbin_ridge = 1e-4,
    negbin_structured_starts = TRUE,
    negbin_start_separation = 0.6,
    negbin_start_jitter = 0.05,
    negbin_init_theta = 8,
    negbin_theta_ridge = 0,
    negbin_theta_target = NULL,
    negbin_theta_min = 1e-3,
    negbin_theta_max = 1e4,
    negbin_mstep_maxit = 100,
    negbin_eta_min = -20,
    negbin_eta_max = 20,
    negbin_eta_margin = 0.25,
    min_component_weight = 1e-8,
    min_sigma = 1e-4,
    max_sigma = 1e4,
    gaussian_sigma_ridge = 1e-6,
    gaussian_mstep_maxit = 100,
    gating_maxit = 50,
    align_margin_tol = 0.05,
    min_separation = 0.08,
    global_slope_weight = 0.25,
    global_gating_weight = 0.25,
    global_gating_slope_weight = 0.05,
    global_sigma_weight = 0.10,
    global_posterior_weight = 1.00,
    warn_ambiguous = TRUE,
    keep_data = TRUE
  )
  utils::modifyList(defaults, control %||% list())
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.vcmoe_joint_path_engine_id <- "joint_path_em"

.vcmoe_joint_path_control <- function(control) {
  joint_path_defaults <- list(
    progress_every = 1L,
    allow_dense_u_grid = FALSE,
    joint_path_dense_grid_limit = 100L,
    joint_path_dense_grid_ratio = 0.2
  )
  .vcmoe_default_control(
    utils::modifyList(joint_path_defaults, control %||% list())
  )
}

.vcmoe_progress_schema <- c(
  "timestamp_utc", "event", "status", "family", "k", "n", "bandwidth",
  "grid_id", "n_grid", "u", "start_id", "n_starts", "iteration", "maxit",
  "elapsed_seconds", "grid_elapsed_seconds", "start_elapsed_seconds",
  "loglik", "loglik_delta", "converged", "selected_start", "error_message"
)

.vcmoe_progress_setup <- function(progress) {
  if (is.null(progress) || identical(progress, FALSE)) {
    return(list(enabled = FALSE))
  }
  if (isTRUE(progress)) {
    return(list(
      enabled = TRUE,
      mode = "message",
      path = NULL,
      start_time = proc.time()[["elapsed"]]
    ))
  }
  if (is.character(progress) && length(progress) == 1L && nzchar(progress)) {
    parent <- dirname(progress)
    if (!dir.exists(parent)) {
      dir.create(parent, recursive = TRUE, showWarnings = FALSE)
    }
    if (!dir.exists(parent)) {
      stop("Could not create progress log directory: ", parent, call. = FALSE)
    }
    return(list(
      enabled = TRUE,
      mode = "file",
      path = progress,
      start_time = proc.time()[["elapsed"]]
    ))
  }
  stop("`progress` must be NULL, FALSE, TRUE, or a single CSV file path.", call. = FALSE)
}

.vcmoe_progress_due <- function(control, iteration, converged = FALSE) {
  every <- suppressWarnings(as.integer(control$progress_every %||% 1L))
  if (!is.finite(every) || every < 1L) {
    every <- 1L
  }
  iteration == 1L || iteration %% every == 0L ||
    isTRUE(converged) || iteration >= as.integer(control$maxit)
}

.vcmoe_progress_log <- function(progress, event, status = NA_character_,
                                family = NA_character_, k = NA_integer_,
                                n = NA_integer_, bandwidth = NA_real_,
                                grid_id = NA_integer_, n_grid = NA_integer_,
                                u = NA_real_, start_id = NA_integer_,
                                n_starts = NA_integer_, iteration = NA_integer_,
                                maxit = NA_integer_, grid_start_time = NA_real_,
                                start_start_time = NA_real_, loglik = NA_real_,
                                loglik_delta = NA_real_, converged = NA,
                                selected_start = NA_integer_,
                                error_message = NA_character_) {
  if (!is.list(progress) || !isTRUE(progress$enabled)) {
    return(invisible(NULL))
  }
  now <- proc.time()[["elapsed"]]
  elapsed <- now - progress$start_time
  grid_elapsed <- if (is.finite(grid_start_time)) now - grid_start_time else NA_real_
  start_elapsed <- if (is.finite(start_start_time)) now - start_start_time else NA_real_
  row <- data.frame(
    timestamp_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    event = event,
    status = status,
    family = family,
    k = k,
    n = n,
    bandwidth = bandwidth,
    grid_id = grid_id,
    n_grid = n_grid,
    u = u,
    start_id = start_id,
    n_starts = n_starts,
    iteration = iteration,
    maxit = maxit,
    elapsed_seconds = elapsed,
    grid_elapsed_seconds = grid_elapsed,
    start_elapsed_seconds = start_elapsed,
    loglik = loglik,
    loglik_delta = loglik_delta,
    converged = converged,
    selected_start = selected_start,
    error_message = error_message,
    stringsAsFactors = FALSE
  )
  row <- row[.vcmoe_progress_schema]

  if (identical(progress$mode, "message")) {
    message(sprintf(
      "vcmoe progress: %s grid=%s/%s start=%s/%s iter=%s/%s loglik=%s status=%s elapsed=%.1fs",
      event,
      ifelse(is.na(grid_id), "-", grid_id),
      ifelse(is.na(n_grid), "-", n_grid),
      ifelse(is.na(start_id), "-", start_id),
      ifelse(is.na(n_starts), "-", n_starts),
      ifelse(is.na(iteration), "-", iteration),
      ifelse(is.na(maxit), "-", maxit),
      ifelse(is.na(loglik), "-", signif(loglik, 6)),
      ifelse(is.na(status), "-", status),
      elapsed
    ))
  } else if (identical(progress$mode, "file")) {
    append <- file.exists(progress$path) && file.info(progress$path)$size > 0
    utils::write.table(
      row,
      file = progress$path,
      sep = ",",
      row.names = FALSE,
      col.names = !append,
      append = append,
      na = ""
    )
  }
  invisible(NULL)
}

.collapse_deparse <- function(x) {
  paste(deparse(x), collapse = "")
}

.split_vcmoe_formula <- function(formula) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must have the form y ~ expert_terms | gating_terms.", call. = FALSE)
  }
  rhs <- formula[[3L]]
  if (!is.call(rhs) || !identical(rhs[[1L]], as.name("|"))) {
    stop("`formula` must contain `|`, for example y ~ z1 + z2 | x1 + x2.", call. = FALSE)
  }

  lhs <- .collapse_deparse(formula[[2L]])
  expert_rhs <- .collapse_deparse(rhs[[2L]])
  gating_rhs <- .collapse_deparse(rhs[[3L]])
  env <- environment(formula)

  expert_formula <- stats::as.formula(paste(lhs, "~", expert_rhs), env = env)
  gating_formula <- stats::as.formula(paste("~", gating_rhs), env = env)

  list(
    response = lhs,
    expert_formula = expert_formula,
    gating_formula = gating_formula
  )
}

.extract_u <- function(u, data) {
  if (is.character(u) && length(u) == 1L) {
    if (!u %in% names(data)) {
      stop("`u` was provided as a name but is not present in `data`.", call. = FALSE)
    }
    return(list(values = data[[u]], name = u))
  }
  if (length(u) != nrow(data)) {
    stop("`u` must be a column name or a vector with one value per row in `data`.", call. = FALSE)
  }
  list(values = u, name = NULL)
}

.default_bandwidth <- function(u) {
  s <- stats::sd(u)
  if (!is.finite(s) || s <= 0) {
    stop("`u` must vary to estimate varying coefficients.", call. = FALSE)
  }
  1.06 * s * length(u)^(-1 / 5)
}

.default_u_grid <- function(u) {
  unique_u <- sort(unique(u))
  if (length(unique_u) <= 20L) {
    return(unique_u)
  }
  as.numeric(stats::quantile(u, probs = seq(0.05, 0.95, length.out = 20L), names = FALSE))
}

.vcmoe_validate_u_scale <- function(u_scale) {
  match.arg(u_scale %||% "unit", c("unit", "none"))
}

.vcmoe_make_u_scaling <- function(u, u_scale = c("unit", "none")) {
  u_scale <- .vcmoe_validate_u_scale(u_scale)
  u <- as.numeric(u)
  finite <- is.finite(u)
  if (!any(finite)) {
    stop("`u` must contain finite values.", call. = FALSE)
  }
  original_range <- range(u[finite])
  original_domain_length <- diff(original_range)
  if (!is.finite(original_domain_length) || original_domain_length <= 0) {
    stop("`u` must vary to estimate varying coefficients.", call. = FALSE)
  }
  if (identical(u_scale, "unit")) {
    analysis_range <- c(0, 1)
    analysis_domain_length <- 1
  } else {
    analysis_range <- original_range
    analysis_domain_length <- original_domain_length
  }
  list(
    method = u_scale,
    original_range = original_range,
    original_domain_length = original_domain_length,
    analysis_range = analysis_range,
    analysis_domain_length = analysis_domain_length,
    domain_length = analysis_domain_length
  )
}

.vcmoe_scale_u <- function(u, scaling) {
  if (is.null(scaling) || identical(scaling$method, "none")) {
    return(as.numeric(u))
  }
  (as.numeric(u) - scaling$original_range[[1L]]) / scaling$original_domain_length
}

.vcmoe_unscale_u <- function(u, scaling) {
  if (is.null(scaling) || identical(scaling$method, "none")) {
    return(as.numeric(u))
  }
  scaling$original_range[[1L]] + as.numeric(u) * scaling$original_domain_length
}

.kernel_weights <- function(u, u0, bandwidth,
                            kernel = c("gaussian", "epanechnikov"),
                            normalization = c("mean_one", "density")) {
  kernel <- match.arg(kernel)
  normalization <- match.arg(normalization)
  du <- (u - u0) / bandwidth
  w <- if (identical(kernel, "gaussian")) {
    stats::dnorm(du)
  } else {
    0.75 * pmax(0, 1 - du^2) * as.numeric(abs(du) <= 1)
  }
  if (!any(is.finite(w)) || sum(w) <= 0) {
    stop("Kernel weights are degenerate; check `u_grid` and `bandwidth`.", call. = FALSE)
  }
  if (identical(normalization, "density")) {
    return(w / bandwidth)
  }
  w / mean(w)
}

.vcmoe_estimation_spec <- function(parameterization = "a1_epanechnikov_scaled") {
  parameterization <- match.arg(parameterization, "a1_epanechnikov_scaled")
  list(
    id = "a1_epanechnikov_scaled",
    kernel = "epanechnikov",
    weight_normalization = "density",
    weight_expression = "0.75 * (1 - ((u - u0) / bandwidth)^2)_+ / bandwidth",
    local_basis = "scaled",
    slope_column = "(u - u0) / bandwidth",
    slope_storage = "scaled",
    slope_scale = "scaled_u_units",
    slope_scaled_conversion = "stored slope is already for (u - u0) / bandwidth",
    inference_mode = "paper_style_epanechnikov_scaled"
  )
}

.local_weights <- function(u, u0, bandwidth, estimation_spec) {
  .kernel_weights(
    u,
    u0,
    bandwidth,
    kernel = estimation_spec$kernel,
    normalization = estimation_spec$weight_normalization
  )
}

.local_du <- function(u, u0, bandwidth, estimation_spec) {
  if (identical(estimation_spec$local_basis, "scaled")) {
    return((u - u0) / bandwidth)
  }
  u - u0
}

.slope_to_raw <- function(slope, bandwidth, estimation_spec) {
  if (identical(estimation_spec$slope_storage, "scaled")) {
    return(slope / bandwidth)
  }
  slope
}

.slope_to_scaled <- function(slope, bandwidth, estimation_spec) {
  if (identical(estimation_spec$slope_storage, "scaled")) {
    return(slope)
  }
  slope * bandwidth
}

.row_log_sum_exp <- function(x) {
  m <- apply(x, 1L, max)
  m + log(rowSums(exp(x - m)))
}

.row_softmax <- function(eta) {
  m <- apply(eta, 1L, max)
  z <- exp(eta - m)
  z / rowSums(z)
}

.center_logits <- function(beta_full) {
  sweep(beta_full, 2L, colMeans(beta_full), FUN = "-")
}

.gating_prob <- function(design, beta_full) {
  .row_softmax(design %*% t(beta_full))
}

.all_permutations <- function(k) {
  if (k > 6L) {
    stop("Exact label alignment currently supports at most 6 components.", call. = FALSE)
  }
  if (k == 1L) {
    return(matrix(1L, nrow = 1L))
  }
  out <- lapply(seq_len(k), function(first) {
    rest <- setdiff(seq_len(k), first)
    tails <- .all_permutations(k - 1L)
    mapped <- matrix(rest[tails], nrow = nrow(tails))
    cbind(first, mapped)
  })
  unname(do.call(rbind, out))
}

.bit_count <- function(mask) {
  count <- 0L
  while (mask > 0L) {
    count <- count + bitwAnd(mask, 1L)
    mask <- bitwShiftR(mask, 1L)
  }
  count
}

.assignment_permutation <- function(costs) {
  k <- nrow(costs)
  if (ncol(costs) != k) {
    stop("Assignment cost matrix must be square.", call. = FALSE)
  }
  n_states <- 2L^k
  best <- rep(Inf, n_states)
  second <- rep(Inf, n_states)
  parent_state <- rep(NA_integer_, n_states)
  parent_row <- rep(NA_integer_, n_states)
  best[[1L]] <- 0

  update_scores <- function(score, state, previous_state, row) {
    if (!is.finite(score)) {
      return(NULL)
    }
    if (score < best[[state]]) {
      second[[state]] <<- best[[state]]
      best[[state]] <<- score
      parent_state[[state]] <<- previous_state
      parent_row[[state]] <<- row
    } else if (score > best[[state]] && score < second[[state]]) {
      second[[state]] <<- score
    } else if (score == best[[state]] && !is.finite(second[[state]])) {
      second[[state]] <<- score
    }
    NULL
  }

  for (mask in 0:(n_states - 1L)) {
    state <- mask + 1L
    if (!is.finite(best[[state]])) {
      next
    }
    column <- .bit_count(mask) + 1L
    if (column > k) {
      next
    }
    for (row in seq_len(k)) {
      bit <- bitwShiftL(1L, row - 1L)
      if (bitwAnd(mask, bit) != 0L) {
        next
      }
      next_mask <- bitwOr(mask, bit)
      next_state <- next_mask + 1L
      update_scores(best[[state]] + costs[row, column], next_state, state, row)
      update_scores(second[[state]] + costs[row, column], next_state, state, row)
    }
  }

  final_state <- n_states
  permutation <- integer(k)
  for (column in seq(from = k, to = 1L)) {
    permutation[[column]] <- parent_row[[final_state]]
    final_state <- parent_state[[final_state]]
  }
  list(permutation = permutation, score = best[[n_states]], second_score = second[[n_states]])
}

.apply_permutation_local <- function(local, permutation) {
  fields <- c("expert_intercept", "expert_slope", "gating_intercept", "gating_slope")
  for (field in fields) {
    local[[field]] <- local[[field]][permutation, , drop = FALSE]
  }
  if (!is.null(local$sigma)) {
    local$sigma <- local$sigma[permutation]
  }
  if (!is.null(local$sigma_slope)) {
    local$sigma_slope <- local$sigma_slope[permutation]
  }
  if (!is.null(local$theta)) {
    local$theta <- local$theta[permutation]
  }
  if (!is.null(local$posterior)) {
    local$posterior <- local$posterior[, permutation, drop = FALSE]
  }
  if (!is.null(local$pi)) {
    local$pi <- local$pi[, permutation, drop = FALSE]
  }
  vector_fields <- c("expert_convergence", "expert_gradient_norm", "expert_coef_norm")
  for (field in vector_fields) {
    if (!is.null(local[[field]])) {
      local[[field]] <- local[[field]][permutation]
    }
  }
  local
}

.compose_permutation <- function(first, second) {
  first[second]
}

.canonical_order_local <- function(local, z_design) {
  mu_bar <- vapply(seq_len(nrow(local$expert_intercept)), function(component) {
    mean(as.numeric(z_design %*% local$expert_intercept[component, ]))
  }, numeric(1L))
  .apply_permutation_local(local, order(mu_bar))
}

.complete_cases_frame <- function(frame, n) {
  if (ncol(frame) == 0L) {
    return(rep(TRUE, n))
  }
  stats::complete.cases(frame)
}

.extract_model_offset <- function(frame) {
  offset <- stats::model.offset(frame)
  if (is.null(offset)) {
    return(rep(0, nrow(frame)))
  }
  offset <- as.numeric(offset)
  if (length(offset) != nrow(frame)) {
    stop("Model offset must have one value per row.", call. = FALSE)
  }
  if (any(!is.finite(offset))) {
    stop("Model offset values must be finite.", call. = FALSE)
  }
  offset
}

.is_whole_count <- function(x) {
  is.finite(x) & abs(x - round(x)) < 1e-8
}

.parse_vcmoe_response <- function(family, raw_response, response_label) {
  if (identical(family, "gaussian")) {
    if (!is.numeric(raw_response) || is.matrix(raw_response)) {
      stop("The Gaussian implementation requires a numeric vector response.", call. = FALSE)
    }
    y <- as.numeric(raw_response)
    if (any(!is.finite(y))) {
      stop("Gaussian response values must be finite.", call. = FALSE)
    }
    return(list(
      type = "gaussian",
      y = y,
      observed = y,
      success = NULL,
      failure = NULL,
      trials = NULL,
      info = list(type = "gaussian", response = response_label)
    ))
  }

  if (identical(family, "negative-binomial")) {
    if (!is.numeric(raw_response) || is.matrix(raw_response)) {
      stop("Negative-Binomial responses must be numeric count vectors.", call. = FALSE)
    }
    y <- as.numeric(raw_response)
    if (any(!.is_whole_count(y))) {
      stop("Negative-Binomial response values must be finite whole numbers.", call. = FALSE)
    }
    if (any(y < 0)) {
      stop("Negative-Binomial response values must be non-negative.", call. = FALSE)
    }
    return(list(
      type = "negative-binomial",
      y = y,
      observed = log1p(y),
      success = NULL,
      failure = NULL,
      trials = NULL,
      info = list(type = "negative-binomial", response = response_label)
    ))
  }

  if (!identical(family, "binomial")) {
    stop("Unsupported family.", call. = FALSE)
  }

  if (is.matrix(raw_response) || is.data.frame(raw_response)) {
    response_matrix <- as.matrix(raw_response)
    if (ncol(response_matrix) != 2L) {
      stop("Grouped binomial responses must have two columns: success and failure.", call. = FALSE)
    }
    success <- as.numeric(response_matrix[, 1L])
    failure <- as.numeric(response_matrix[, 2L])
    if (any(!.is_whole_count(success)) || any(!.is_whole_count(failure))) {
      stop("Grouped binomial success and failure counts must be finite whole numbers.", call. = FALSE)
    }
    if (any(success < 0) || any(failure < 0)) {
      stop("Grouped binomial success and failure counts must be non-negative.", call. = FALSE)
    }
    trials <- success + failure
    if (any(trials <= 0)) {
      stop("Grouped binomial rows must have positive trial counts.", call. = FALSE)
    }
    response_names <- colnames(response_matrix)
    success_name <- if (!is.null(response_names) && length(response_names) >= 1L) response_names[[1L]] else NULL
    failure_name <- if (!is.null(response_names) && length(response_names) >= 2L) response_names[[2L]] else NULL
    return(list(
      type = "grouped",
      y = success,
      observed = success / trials,
      success = success,
      failure = failure,
      trials = trials,
      info = list(
        type = "grouped",
        response = response_label,
        success = success_name,
        failure = failure_name
      )
    ))
  }

  if (is.logical(raw_response)) {
    raw_response <- as.numeric(raw_response)
  }
  if (!is.numeric(raw_response)) {
    stop("Bernoulli responses must be numeric, integer, or logical 0/1 values.", call. = FALSE)
  }
  success <- as.numeric(raw_response)
  if (any(!is.finite(success))) {
    stop("Bernoulli response values must be finite 0/1 values.", call. = FALSE)
  }
  if (any(!(success %in% c(0, 1)))) {
    stop("Bernoulli response values must be 0/1.", call. = FALSE)
  }
  trials <- rep(1, length(success))
  list(
    type = "bernoulli",
    y = success,
    observed = success,
    success = success,
    failure = trials - success,
    trials = trials,
    info = list(type = "bernoulli", response = response_label)
  )
}

.component_separation <- function(local, z_design, y) {
  k <- nrow(local$expert_intercept)
  if (k < 2L) {
    return(Inf)
  }
  scale_y <- max(stats::sd(y), 1e-6)
  values <- numeric(0L)
  for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      mu_a <- as.numeric(z_design %*% local$expert_intercept[a, ])
      mu_b <- as.numeric(z_design %*% local$expert_intercept[b, ])
      values <- c(values, sqrt(mean((mu_a - mu_b)^2)) / scale_y)
    }
  }
  min(values)
}

.alignment_cost_matrix <- function(current, reference, z_design, x_design, y) {
  k <- nrow(current$expert_intercept)
  scale_y <- max(stats::sd(y), 1e-6)
  costs <- matrix(0, nrow = k, ncol = k)
  for (cur in seq_len(k)) {
    mu_cur <- as.numeric(z_design %*% current$expert_intercept[cur, ])
    gate_cur <- as.numeric(x_design %*% current$gating_intercept[cur, ])
    for (ref in seq_len(k)) {
      mu_ref <- as.numeric(z_design %*% reference$expert_intercept[ref, ])
      gate_ref <- as.numeric(x_design %*% reference$gating_intercept[ref, ])
      expert_cost <- mean((mu_cur - mu_ref)^2) / (scale_y^2)
      gating_cost <- mean((gate_cur - gate_ref)^2)
      sigma_cost <- 0
      if (!is.null(current$sigma) && !is.null(reference$sigma)) {
        sigma_cost <- (log(current$sigma[cur]) - log(reference$sigma[ref]))^2
        if (!is.null(current$sigma_slope) && !is.null(reference$sigma_slope)) {
          sigma_cost <- sigma_cost + 0.25 * (current$sigma_slope[cur] - reference$sigma_slope[ref])^2
        }
      }
      theta_cost <- 0
      if (!is.null(current$theta) && !is.null(reference$theta)) {
        theta_cost <- (log(current$theta[cur]) - log(reference$theta[ref]))^2
      }
      posterior_cost <- 0
      if (!is.null(current$posterior) && !is.null(reference$posterior)) {
        posterior_cost <- mean((current$posterior[, cur] - reference$posterior[, ref])^2)
      }
      costs[cur, ref] <- expert_cost + 0.25 * gating_cost +
        0.1 * sigma_cost + 0.1 * theta_cost + posterior_cost
    }
  }
  costs
}

.global_transition_cost <- function(previous, current, delta_u, z_design, x_design, y, control,
                                    bandwidth = 1, estimation_spec = .vcmoe_estimation_spec()) {
  scale_y <- max(stats::sd(y), 1e-6)
  k <- nrow(previous$expert_intercept)
  expert_cost <- 0
  expert_slope_cost <- 0
  gating_cost <- 0
  gating_slope_cost <- 0
  previous_expert_slope <- .slope_to_raw(previous$expert_slope, bandwidth, estimation_spec)
  current_expert_slope <- .slope_to_raw(current$expert_slope, bandwidth, estimation_spec)
  previous_gating_slope <- .slope_to_raw(previous$gating_slope, bandwidth, estimation_spec)
  current_gating_slope <- .slope_to_raw(current$gating_slope, bandwidth, estimation_spec)

  for (component in seq_len(k)) {
    predicted_expert <- previous$expert_intercept[component, ] +
      delta_u * previous_expert_slope[component, ]
    expert_cost <- expert_cost + mean(
      as.numeric(z_design %*% (current$expert_intercept[component, ] - predicted_expert))^2
    ) / (scale_y^2)
    expert_slope_cost <- expert_slope_cost + mean(
      as.numeric(z_design %*% (current_expert_slope[component, ] - previous_expert_slope[component, ]))^2
    ) / (scale_y^2)

    predicted_gating <- previous$gating_intercept[component, ] +
      delta_u * previous_gating_slope[component, ]
    gating_cost <- gating_cost + mean(
      as.numeric(x_design %*% (current$gating_intercept[component, ] - predicted_gating))^2
    )
    gating_slope_cost <- gating_slope_cost + mean(
      as.numeric(x_design %*% (current_gating_slope[component, ] - previous_gating_slope[component, ]))^2
    )
  }

  sigma_cost <- 0
  if (!is.null(current$sigma) && !is.null(previous$sigma)) {
    previous_sigma_slope <- if (!is.null(previous$sigma_slope)) {
      .slope_to_raw(previous$sigma_slope, bandwidth, estimation_spec)
    } else {
      rep(0, length(previous$sigma))
    }
    current_sigma_slope <- if (!is.null(current$sigma_slope)) {
      .slope_to_raw(current$sigma_slope, bandwidth, estimation_spec)
    } else {
      rep(0, length(current$sigma))
    }
    predicted_log_sigma <- log(previous$sigma) + delta_u * previous_sigma_slope
    sigma_cost <- mean((log(current$sigma) - predicted_log_sigma)^2) +
      0.25 * mean((current_sigma_slope - previous_sigma_slope)^2)
  }
  theta_cost <- 0
  if (!is.null(current$theta) && !is.null(previous$theta)) {
    theta_cost <- mean((log(current$theta) - log(previous$theta))^2)
  }
  posterior_cost <- 0
  if (!is.null(current$posterior) && !is.null(previous$posterior)) {
    posterior_cost <- mean((current$posterior - previous$posterior)^2)
  }

  expert_cost / k +
    control$global_slope_weight * expert_slope_cost / k +
    control$global_gating_weight * gating_cost / k +
    control$global_gating_slope_weight * gating_slope_cost / k +
    control$global_sigma_weight * sigma_cost +
    control$global_sigma_weight * theta_cost +
    control$global_posterior_weight * posterior_cost
}

.global_align_local_fits <- function(locals, u_grid, z_design, x_design, y, control,
                                     bandwidth = 1,
                                     estimation_spec = .vcmoe_estimation_spec()) {
  n_grid <- length(locals)
  if (n_grid <= 1L) {
    locals[[1L]]$global_permutation <- seq_len(nrow(locals[[1L]]$expert_intercept))
    locals[[1L]]$global_transition_cost <- NA_real_
    return(list(locals = locals, path_cost = 0, transition_cost = NA_real_))
  }

  k <- nrow(locals[[1L]]$expert_intercept)
  perms <- .all_permutations(k)
  n_states <- nrow(perms)
  identity_state <- which(apply(perms, 1L, function(perm) all(perm == seq_len(k))))[1L]

  permuted <- lapply(seq_len(n_grid), function(grid_id) {
    lapply(seq_len(n_states), function(state_id) {
      .apply_permutation_local(locals[[grid_id]], perms[state_id, ])
    })
  })

  dp <- matrix(Inf, nrow = n_grid, ncol = n_states)
  back <- matrix(NA_integer_, nrow = n_grid, ncol = n_states)
  dp[1L, identity_state] <- 0

  transition_cache <- vector("list", n_grid)
  for (grid_id in 2:n_grid) {
    delta_u <- u_grid[grid_id] - u_grid[grid_id - 1L]
    transition_cache[[grid_id]] <- matrix(NA_real_, nrow = n_states, ncol = n_states)
    previous_states <- which(is.finite(dp[grid_id - 1L, ]))
    for (state_id in seq_len(n_states)) {
      best_value <- Inf
      best_previous <- NA_integer_
      for (previous_state in previous_states) {
        if (is.na(transition_cache[[grid_id]][previous_state, state_id])) {
          transition_cache[[grid_id]][previous_state, state_id] <- .global_transition_cost(
            permuted[[grid_id - 1L]][[previous_state]],
            permuted[[grid_id]][[state_id]],
            delta_u,
            z_design,
            x_design,
            y,
            control,
            bandwidth,
            estimation_spec
          )
        }
        candidate <- dp[grid_id - 1L, previous_state] +
          transition_cache[[grid_id]][previous_state, state_id]
        if (candidate < best_value) {
          best_value <- candidate
          best_previous <- previous_state
        }
      }
      dp[grid_id, state_id] <- best_value
      back[grid_id, state_id] <- best_previous
    }
  }

  path_states <- integer(n_grid)
  path_states[n_grid] <- which.min(dp[n_grid, ])
  for (grid_id in seq(from = n_grid, to = 2L)) {
    path_states[grid_id - 1L] <- back[grid_id, path_states[grid_id]]
  }

  transition_cost <- rep(NA_real_, n_grid)
  aligned <- vector("list", n_grid)
  for (grid_id in seq_len(n_grid)) {
    extra_perm <- perms[path_states[grid_id], ]
    aligned[[grid_id]] <- .apply_permutation_local(locals[[grid_id]], extra_perm)
    aligned[[grid_id]]$global_permutation <- extra_perm
    if (!is.null(locals[[grid_id]]$permutation)) {
      aligned[[grid_id]]$permutation <- .compose_permutation(locals[[grid_id]]$permutation, extra_perm)
    }
    if (grid_id > 1L) {
      transition_cost[grid_id] <- transition_cache[[grid_id]][
        path_states[grid_id - 1L],
        path_states[grid_id]
      ]
    }
    aligned[[grid_id]]$global_transition_cost <- transition_cost[grid_id]
  }

  list(
    locals = aligned,
    path_cost = dp[n_grid, path_states[n_grid]],
    transition_cost = transition_cost
  )
}

.align_local_fit <- function(current, reference, z_design, x_design, y, control) {
  if (is.null(reference)) {
    permutation <- order(vapply(seq_len(nrow(current$expert_intercept)), function(component) {
      mean(as.numeric(z_design %*% current$expert_intercept[component, ]))
    }, numeric(1L)))
    current <- .apply_permutation_local(current, permutation)
    current$permutation <- permutation
    current$greedy_permutation <- permutation
    current$alignment_score <- NA_real_
    current$alignment_margin <- Inf
    current$ambiguous <- .component_separation(current, z_design, y) < control$min_separation
    return(current)
  }

  costs <- .alignment_cost_matrix(current, reference, z_design, x_design, y)
  if (nrow(costs) <= 6L) {
    perms <- .all_permutations(nrow(costs))
    scores <- apply(perms, 1L, function(perm) {
      sum(costs[cbind(perm, seq_len(ncol(costs)))])
    })
    ord <- order(scores)
    best_perm <- perms[ord[1L], ]
    best_score <- scores[ord[1L]]
    second_score <- if (length(ord) > 1L) scores[ord[2L]] else Inf
  } else {
    assignment <- .assignment_permutation(costs)
    best_perm <- assignment$permutation
    best_score <- assignment$score
    second_score <- assignment$second_score
  }
  margin <- second_score - best_score
  relative_margin <- margin / (abs(best_score) + 1e-8)
  separation <- .component_separation(current, z_design, y)

  aligned <- .apply_permutation_local(current, best_perm)
  aligned$permutation <- best_perm
  aligned$greedy_permutation <- best_perm
  aligned$alignment_score <- best_score
  aligned$alignment_margin <- margin
  aligned$ambiguous <- (is.finite(relative_margin) &&
    relative_margin < control$align_margin_tol) || separation < control$min_separation
  aligned
}

.weighted_entropy <- function(posterior, weights = NULL) {
  p <- pmax(posterior, 1e-12)
  entropy <- -rowSums(p * log(p))
  if (is.null(weights)) {
    return(mean(entropy))
  }
  weights <- pmax(as.numeric(weights), 0)
  if (length(weights) != length(entropy) || sum(weights) <= 0) {
    return(NA_real_)
  }
  sum(weights * entropy) / sum(weights)
}

.fit_wls <- function(y, design, weights, ridge) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= 1e-10) {
    return(rep(0, ncol(design)))
  }
  xtx <- crossprod(design, design * weights)
  rhs <- crossprod(design, y * weights)
  diag(xtx) <- diag(xtx) + ridge
  out <- tryCatch(
    solve(xtx, rhs),
    error = function(e) qr.solve(xtx, rhs)
  )
  as.numeric(out)
}
