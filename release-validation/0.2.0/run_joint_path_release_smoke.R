#!/usr/bin/env Rscript

# Lightweight paired release smoke for the local-grid and joint-path engines.
# Run from the package root with:
#   Rscript release-validation/0.2.0/run_joint_path_release_smoke.R

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
description_path <- file.path(root, "DESCRIPTION")
if (!file.exists(description_path)) {
  stop("Run this script from the VCMoE package root.", call. = FALSE)
}
description <- read.dcf(description_path, fields = c("Package", "Version"))
if (!identical(unname(description[1L, "Package"]), "VCMoE")) {
  stop("The current directory is not the VCMoE package root.", call. = FALSE)
}
release_version <- unname(description[1L, "Version"])
if (!identical(release_version, "0.2.0")) {
  warning(sprintf("Expected VCMoE 0.2.0, found %s.", release_version), call. = FALSE)
}

results_dir <- file.path(root, "release-validation", "0.2.0", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
replicate_path <- file.path(results_dir, "joint_path_release_smoke_replicates.csv")
summary_path <- file.path(results_dir, "joint_path_release_smoke_summary.csv")
plot_path <- file.path(results_dir, "joint_path_release_smoke.png")

load_release_source <- function(attempts = 5L, delay_seconds = 2) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package `pkgload` is required to run this source-tree smoke test.", call. = FALSE)
  }
  errors <- character(0L)
  for (attempt in seq_len(attempts)) {
    loaded <- tryCatch(
      {
        pkgload::load_all(
          root,
          reset = TRUE,
          recompile = FALSE,
          export_all = FALSE,
          helpers = FALSE,
          quiet = TRUE
        )
        TRUE
      },
      error = function(error) {
        errors <<- c(errors, conditionMessage(error))
        FALSE
      }
    )
    if (loaded) {
      return(invisible(TRUE))
    }
    if (attempt < attempts) {
      message(sprintf(
        "Package source did not load on attempt %d/%d; retrying in %.1f seconds.",
        attempt, attempts, delay_seconds
      ))
      Sys.sleep(delay_seconds)
    }
  }
  stop(
    sprintf("Could not load the package source after %d attempts: %s",
            attempts, paste(unique(errors), collapse = " | ")),
    call. = FALSE
  )
}

load_release_source()

scenarios <- list(
  list(
    id = "gaussian_k2", label = "Gaussian K=2", family = "gaussian",
    k = 2L, n = 120L, repeats = 2L, separation = 1.7,
    bandwidth = 0.45, u_grid = seq(0.15, 0.85, length.out = 5L), maxit = 40L
  ),
  list(
    id = "grouped_binomial_k2", label = "Grouped Binomial K=2", family = "binomial",
    k = 2L, n = 120L, repeats = 2L, separation = 1.7, trials = 8L,
    bandwidth = 0.45, u_grid = seq(0.15, 0.85, length.out = 5L), maxit = 80L
  ),
  list(
    id = "negbin_offset_k2", label = "NegBin + offset K=2", family = "negative-binomial",
    k = 2L, n = 120L, repeats = 2L, separation = 1.5,
    bandwidth = 0.45, u_grid = seq(0.20, 0.80, length.out = 4L), maxit = 60L
  ),
  list(
    id = "gaussian_k3", label = "Gaussian K=3", family = "gaussian",
    k = 3L, n = 120L, repeats = 1L, separation = 1.8,
    bandwidth = 0.45, u_grid = seq(0.20, 0.80, length.out = 4L), maxit = 80L
  ),
  list(
    id = "gaussian_k10", label = "Gaussian K=10", family = "gaussian",
    k = 10L, n = 180L, repeats = 1L, separation = 1.9,
    bandwidth = 0.50, u_grid = c(0.25, 0.50, 0.75), maxit = 20L
  )
)
engines <- c("local_grid_em", "joint_path_em")

safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) else NA_real_
}

safe_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

path_curvature_rms <- function(value, u_grid) {
  if (is.null(value) || length(u_grid) < 3L || dim(value)[1L] != length(u_grid)) {
    return(NA_real_)
  }
  path <- matrix(value, nrow = length(u_grid))
  delta_u <- diff(u_grid)
  slope <- diff(path) / delta_u
  midpoint_delta <- (delta_u[-1L] + delta_u[-length(delta_u)]) / 2
  curvature <- diff(slope) / midpoint_delta
  safe_values <- curvature[is.finite(curvature)]
  if (length(safe_values)) sqrt(mean(safe_values^2)) else NA_real_
}

adjacent_path_rms <- function(value, order_index) {
  if (is.null(value) || nrow(value) < 2L) {
    return(NA_real_)
  }
  delta <- diff(value[order_index, , drop = FALSE])
  safe_values <- delta[is.finite(delta)]
  if (length(safe_values)) sqrt(mean(safe_values^2)) else NA_real_
}

data_signature <- function(data) {
  numeric_columns <- vapply(data, is.numeric, logical(1L))
  values <- unlist(data[numeric_columns], use.names = FALSE)
  weights <- (seq_along(values) %% 997L) + 1L
  sprintf("%.12e", sum(values * weights))
}

simulate_scenario <- function(scenario, seed) {
  if (identical(scenario$family, "gaussian")) {
    return(VCMoE::simulate_vcmoe_gaussian(
      n = scenario$n,
      k = scenario$k,
      seed = seed,
      separation = scenario$separation
    ))
  }
  if (identical(scenario$family, "binomial")) {
    return(VCMoE::simulate_vcmoe_binomial(
      n = scenario$n,
      k = scenario$k,
      seed = seed,
      separation = scenario$separation,
      trials = scenario$trials
    ))
  }
  VCMoE::simulate_vcmoe_negbin(
    n = scenario$n,
    k = scenario$k,
    seed = seed,
    separation = scenario$separation
  )
}

fit_formula <- function(family) {
  if (identical(family, "binomial")) {
    return(cbind(success, failure) ~ z1 | x1)
  }
  if (identical(family, "negative-binomial")) {
    return(y ~ z1 + offset(log_size_factor) | x1)
  }
  y ~ z1 | x1
}

truth_marginal_mean <- function(simulation, family) {
  if (identical(family, "binomial")) {
    return(rowSums(simulation$truth$probability * simulation$truth$success_probability))
  }
  rowSums(simulation$truth$probability * simulation$truth$mean)
}

response_free_data <- function(data) {
  response_columns <- intersect(c("y", "success", "failure"), names(data))
  data[response_columns] <- NULL
  data
}

empty_result <- function(scenario, replicate, pair_id, data_seed, fit_seed,
                         signature, engine) {
  data.frame(
    release_version = release_version,
    scenario_id = scenario$id,
    scenario_label = scenario$label,
    pair_id = pair_id,
    replicate = replicate,
    data_seed = data_seed,
    fit_seed = fit_seed,
    data_signature = signature,
    family = scenario$family,
    k = scenario$k,
    n = scenario$n,
    n_grid = length(scenario$u_grid),
    bandwidth = scenario$bandwidth,
    maxit = scenario$maxit,
    engine = engine,
    success = FALSE,
    finite_outputs = FALSE,
    runtime_seconds = NA_real_,
    converged = FALSE,
    grid_convergence_fraction = NA_real_,
    mean_iterations = NA_real_,
    max_iterations = NA_real_,
    min_effective_component_size = NA_real_,
    posterior_entropy = NA_real_,
    prediction_rmse = NA_real_,
    expert_path_roughness = NA_real_,
    gating_path_roughness = NA_real_,
    posterior_path_roughness = NA_real_,
    prior_path_roughness = NA_real_,
    posterior_row_sum_max_error = NA_real_,
    warning_count = 0L,
    warnings = "",
    error = "",
    stringsAsFactors = FALSE
  )
}

run_fit <- function(simulation, scenario, replicate, pair_id, data_seed,
                    fit_seed, signature, engine) {
  result <- empty_result(
    scenario, replicate, pair_id, data_seed, fit_seed, signature, engine
  )
  warnings <- character(0L)
  started <- proc.time()[["elapsed"]]

  evaluated <- tryCatch(
    withCallingHandlers(
      {
        fit <- VCMoE::vcmoe_fit(
          fit_formula(scenario$family),
          data = simulation$data,
          u = "u",
          k = scenario$k,
          family = scenario$family,
          bandwidth = scenario$bandwidth,
          u_grid = scenario$u_grid,
          engine = engine,
          control = list(
            maxit = scenario$maxit,
            tol = 1e-4,
            n_starts = 1L,
            seed = fit_seed,
            warn_ambiguous = FALSE,
            keep_data = TRUE,
            gating_maxit = 40L,
            gaussian_mstep_maxit = 60L,
            negbin_mstep_maxit = 60L,
            negbin_theta_ridge = 0.05,
            negbin_theta_target = 8
          )
        )
        diagnostics <- VCMoE::vcmoe_diagnostics(fit)
        prediction_data <- response_free_data(simulation$data)
        marginal_prediction <- stats::predict(
          fit, newdata = prediction_data, type = "mean"
        )
        posterior <- stats::predict(fit, type = "posterior")
        prior <- stats::predict(
          fit, newdata = prediction_data, type = "prior"
        )
        coefficient_values <- unlist(
          fit$coefficients[!vapply(fit$coefficients, is.null, logical(1L))],
          recursive = TRUE,
          use.names = FALSE
        )
        finite_outputs <- length(coefficient_values) > 0L &&
          all(is.finite(coefficient_values)) &&
          all(is.finite(marginal_prediction)) &&
          all(is.finite(posterior)) &&
          all(is.finite(prior)) &&
          all(is.finite(diagnostics$min_component_effective_n))
        grid_convergence <- diagnostics$converged %in% TRUE
        engine_converged <- if (identical(engine, "joint_path_em")) {
          isTRUE(fit$diagnostics$joint_path_converged)
        } else {
          all(grid_convergence)
        }
        order_index <- order(simulation$data$u)
        list(
          finite_outputs = finite_outputs,
          converged = engine_converged,
          grid_convergence_fraction = mean(grid_convergence),
          mean_iterations = safe_mean(diagnostics$iterations),
          max_iterations = safe_max(diagnostics$iterations),
          min_effective_component_size = safe_min(
            diagnostics$min_component_effective_n
          ),
          posterior_entropy = safe_mean(diagnostics$posterior_entropy),
          prediction_rmse = sqrt(mean(
            (marginal_prediction - truth_marginal_mean(
              simulation, scenario$family
            ))^2
          )),
          expert_path_roughness = path_curvature_rms(
            fit$coefficients$expert, fit$u_grid
          ),
          gating_path_roughness = path_curvature_rms(
            fit$coefficients$gating, fit$u_grid
          ),
          posterior_path_roughness = adjacent_path_rms(
            posterior, order_index
          ),
          prior_path_roughness = adjacent_path_rms(prior, order_index),
          posterior_row_sum_max_error = max(abs(rowSums(posterior) - 1))
        )
      },
      warning = function(warning_condition) {
        warnings <<- c(warnings, conditionMessage(warning_condition))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error_condition) error_condition
  )

  result$runtime_seconds <- proc.time()[["elapsed"]] - started
  warnings <- unique(gsub("[\r\n]+", " ", warnings))
  result$warning_count <- length(warnings)
  result$warnings <- paste(warnings, collapse = " | ")
  if (inherits(evaluated, "error")) {
    result$error <- gsub("[\r\n]+", " ", conditionMessage(evaluated))
    return(result)
  }

  result$success <- TRUE
  for (name in names(evaluated)) {
    result[[name]] <- evaluated[[name]]
  }
  result
}

rows <- list()
row_id <- 0L
for (scenario_index in seq_along(scenarios)) {
  scenario <- scenarios[[scenario_index]]
  for (replicate in seq_len(scenario$repeats)) {
    data_seed <- 20000L + 100L * scenario_index + replicate
    fit_seed <- 50000L + 100L * scenario_index + replicate
    pair_id <- sprintf("%s_rep%02d", scenario$id, replicate)
    simulation <- simulate_scenario(scenario, data_seed)
    signature <- data_signature(simulation$data)
    message(sprintf(
      "[%s] n=%d, K=%d, seed=%d", pair_id, scenario$n, scenario$k, data_seed
    ))
    for (engine in engines) {
      message(sprintf("  fitting %s", engine))
      row_id <- row_id + 1L
      rows[[row_id]] <- run_fit(
        simulation = simulation,
        scenario = scenario,
        replicate = replicate,
        pair_id = pair_id,
        data_seed = data_seed,
        fit_seed = fit_seed,
        signature = signature,
        engine = engine
      )
    }
  }
}
replicates <- do.call(rbind, rows)
rownames(replicates) <- NULL

pair_counts <- table(replicates$pair_id)
if (any(pair_counts != length(engines))) {
  stop("Internal error: each simulation pair must contain exactly two engines.", call. = FALSE)
}
signature_counts <- vapply(
  split(replicates$data_signature, replicates$pair_id),
  function(value) length(unique(value)),
  integer(1L)
)
if (any(signature_counts != 1L)) {
  stop("Internal error: paired engines did not receive identical data.", call. = FALSE)
}

summary_groups <- split(
  replicates,
  interaction(replicates$scenario_id, replicates$engine, drop = TRUE)
)
summary_rows <- lapply(summary_groups, function(group) {
  data.frame(
    release_version = group$release_version[[1L]],
    scenario_id = group$scenario_id[[1L]],
    scenario_label = group$scenario_label[[1L]],
    family = group$family[[1L]],
    k = group$k[[1L]],
    n = group$n[[1L]],
    n_grid = group$n_grid[[1L]],
    bandwidth = group$bandwidth[[1L]],
    maxit = group$maxit[[1L]],
    engine = group$engine[[1L]],
    n_replicates = nrow(group),
    success_rate = mean(group$success %in% TRUE),
    finite_output_rate = mean(group$finite_outputs %in% TRUE),
    convergence_rate = mean(group$converged %in% TRUE),
    grid_convergence_fraction_mean = safe_mean(group$grid_convergence_fraction),
    runtime_seconds_mean = safe_mean(group$runtime_seconds),
    runtime_seconds_sd = safe_sd(group$runtime_seconds),
    min_effective_component_size_mean = safe_mean(
      group$min_effective_component_size
    ),
    min_effective_component_size_min = safe_min(
      group$min_effective_component_size
    ),
    posterior_entropy_mean = safe_mean(group$posterior_entropy),
    prediction_rmse_mean = safe_mean(group$prediction_rmse),
    prediction_rmse_sd = safe_sd(group$prediction_rmse),
    expert_path_roughness_mean = safe_mean(group$expert_path_roughness),
    gating_path_roughness_mean = safe_mean(group$gating_path_roughness),
    posterior_path_roughness_mean = safe_mean(group$posterior_path_roughness),
    prior_path_roughness_mean = safe_mean(group$prior_path_roughness),
    warning_fit_fraction = mean(group$warning_count > 0L),
    stringsAsFactors = FALSE
  )
})
summary_results <- do.call(rbind, summary_rows)
scenario_order <- vapply(scenarios, `[[`, character(1L), "id")
summary_results <- summary_results[
  order(match(summary_results$scenario_id, scenario_order),
        match(summary_results$engine, engines)),
  , drop = FALSE
]
rownames(summary_results) <- NULL

utils::write.csv(replicates, replicate_path, row.names = FALSE, na = "")
utils::write.csv(summary_results, summary_path, row.names = FALSE, na = "")

engine_colors <- c(local_grid_em = "#2364AA", joint_path_em = "#C44536")
short_labels <- c(
  gaussian_k2 = "Gaussian K=2",
  grouped_binomial_k2 = "Binomial K=2",
  negbin_offset_k2 = "NegBin K=2",
  gaussian_k3 = "Gaussian K=3",
  gaussian_k10 = "Gaussian K=10"
)

plot_metric <- function(metric, y_label, log10_scale = FALSE, limits = NULL,
                        legend = FALSE) {
  values <- replicates[[metric]]
  transformed <- if (log10_scale) log10(pmax(values, 1e-10)) else values
  finite_values <- transformed[is.finite(transformed)]
  if (is.null(limits)) {
    if (length(finite_values)) {
      limits <- range(finite_values)
      padding <- max(diff(limits) * 0.12, 0.05)
      limits <- limits + c(-padding, padding)
    } else {
      limits <- c(0, 1)
    }
  }
  graphics::plot(
    NA_real_, NA_real_,
    xlim = c(0.5, length(scenario_order) + 0.5),
    ylim = limits,
    xaxt = "n",
    xlab = "",
    ylab = if (log10_scale) paste0("log10 ", y_label) else y_label,
    bty = "l"
  )
  graphics::axis(
    1,
    at = seq_along(scenario_order),
    labels = unname(short_labels[scenario_order]),
    las = 2,
    cex.axis = 0.75
  )
  for (pair_id in unique(replicates$pair_id)) {
    pair_rows <- replicates[replicates$pair_id == pair_id, , drop = FALSE]
    pair_rows <- pair_rows[match(engines, pair_rows$engine), , drop = FALSE]
    scenario_position <- match(pair_rows$scenario_id[[1L]], scenario_order)
    x <- scenario_position + c(-0.11, 0.11)
    y <- if (log10_scale) {
      log10(pmax(pair_rows[[metric]], 1e-10))
    } else {
      pair_rows[[metric]]
    }
    if (all(is.finite(y))) {
      graphics::segments(x[[1L]], y[[1L]], x[[2L]], y[[2L]], col = "#A7A7A7")
    }
    graphics::points(
      x,
      y,
      pch = 19,
      col = unname(engine_colors[pair_rows$engine]),
      cex = 0.9
    )
  }
  graphics::grid(nx = NA, ny = NULL, col = "#E7E7E7")
  if (legend) {
    graphics::legend(
      "topright",
      legend = c("local_grid_em", "joint_path_em"),
      col = unname(engine_colors[engines]),
      pch = 19,
      bty = "n",
      cex = 0.75
    )
  }
}

grDevices::png(plot_path, width = 1800, height = 1200, res = 160)
old_par <- graphics::par(
  mfrow = c(2, 3),
  mar = c(6.8, 4.2, 2.2, 0.8),
  oma = c(0.5, 0.5, 2.0, 0.5),
  las = 1
)
plot_metric("prediction_rmse", "prediction RMSE", log10_scale = TRUE, legend = TRUE)
plot_metric("runtime_seconds", "runtime (seconds)", log10_scale = TRUE)
plot_metric("grid_convergence_fraction", "converged grid fraction", limits = c(-0.05, 1.05))
plot_metric("min_effective_component_size", "minimum component ESS")
plot_metric("posterior_entropy", "posterior entropy")
plot_metric("expert_path_roughness", "expert path roughness", log10_scale = TRUE)
graphics::mtext(
  sprintf("VCMoE %s paired engine release smoke", release_version),
  outer = TRUE,
  side = 3,
  line = 0.4,
  font = 2
)
graphics::par(old_par)
grDevices::dev.off()

message(sprintf("Wrote %s", replicate_path))
message(sprintf("Wrote %s", summary_path))
message(sprintf("Wrote %s", plot_path))

required_ok <- replicates$success %in% TRUE & replicates$finite_outputs %in% TRUE
if (!all(required_ok)) {
  failed <- replicates[!required_ok, c("pair_id", "engine", "error"), drop = FALSE]
  print(failed, row.names = FALSE)
  stop("Release smoke had failed or non-finite fits; inspect the replicate CSV.", call. = FALSE)
}
posterior_ok <- is.finite(replicates$posterior_row_sum_error) &
  replicates$posterior_row_sum_error <= 1e-8
component_ok <- is.finite(replicates$min_effective_component_size) &
  replicates$min_effective_component_size >= 1
entropy_ok <- is.finite(replicates$posterior_entropy) &
  replicates$posterior_entropy >= 0
if (!all(posterior_ok & component_ok & entropy_ok)) {
  failed <- replicates[
    !(posterior_ok & component_ok & entropy_ok),
    c(
      "pair_id", "engine", "posterior_row_sum_error",
      "min_effective_component_size", "posterior_entropy"
    ),
    drop = FALSE
  ]
  print(failed, row.names = FALSE)
  stop("Release smoke failed posterior/component stability checks.", call. = FALSE)
}

paired_rmse <- reshape(
  replicates[, c("pair_id", "engine", "prediction_rmse")],
  idvar = "pair_id",
  timevar = "engine",
  direction = "wide"
)
rmse_ratio <- paired_rmse$prediction_rmse.joint_path_em /
  paired_rmse$prediction_rmse.local_grid_em
if (!all(is.finite(rmse_ratio)) || stats::median(rmse_ratio) > 1.10) {
  stop(
    "Joint-path median paired prediction RMSE exceeds the 10% release tolerance.",
    call. = FALSE
  )
}
message(sprintf(
  paste0(
    "Release smoke passed: %d/%d finite fits; posterior/component checks passed; ",
    "median joint/local prediction RMSE ratio %.3f."
  ),
  sum(required_ok), length(required_ok), stats::median(rmse_ratio)
))
