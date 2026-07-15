#!/usr/bin/env Rscript

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The release fixture generator requires pkgload.", call. = FALSE)
}

release_root <- normalizePath(getwd(), mustWork = TRUE)
source_root <- Sys.getenv("VCMOE_V010_SOURCE_ROOT", "/private/tmp/VCMoE-v010-fixture")
source_root <- normalizePath(source_root, mustWork = TRUE)
source_commit <- trimws(system2(
  "git", c("-C", shQuote(source_root), "rev-parse", "HEAD"), stdout = TRUE
))
if (!identical(source_commit, "559d717343f24f5794148de6c9e9cf8e12b4b53a")) {
  stop("Fixture source must be the exact v0.1.0 commit 559d717.", call. = FALSE)
}
description <- read.dcf(file.path(source_root, "DESCRIPTION"))
if (!identical(unname(description[1L, "Version"]), "0.1.0")) {
  stop("Fixture source DESCRIPTION is not VCMoE 0.1.0.", call. = FALSE)
}

pkgload::load_all(source_root, quiet = TRUE)
options(digits = 17)

fixture_dir <- file.path(release_root, "tests", "testthat", "fixtures", "v010")
dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)

array_long <- function(x, case, set, term = NULL) {
  if (is.null(x)) {
    return(data.frame())
  }
  if (length(dim(x)) == 3L) {
    grid <- expand.grid(
      u = as.numeric(dimnames(x)[[1L]]),
      component = dimnames(x)[[2L]],
      term = dimnames(x)[[3L]],
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    grid <- expand.grid(
      u = as.numeric(dimnames(x)[[1L]]),
      component = dimnames(x)[[2L]],
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    grid$term <- term
  }
  grid$case <- case
  grid$coefficient_set <- set
  grid$estimate <- as.vector(x)
  grid[, c("case", "coefficient_set", "u", "component", "term", "estimate")]
}

prediction_long <- function(fit, data, dataset) {
  pieces <- VCMoE:::.predict_components(fit, newdata = data)
  posterior <- predict(fit, newdata = data, type = "posterior")
  mean_value <- predict(fit, newdata = data, type = "mean")
  component_mean <- predict(fit, newdata = data, type = "component")
  component_names <- colnames(pieces$prior)
  rows <- list(data.frame(
    dataset = dataset,
    row_id = seq_len(nrow(data)),
    prediction_type = "mean",
    component = NA_character_,
    value = as.numeric(mean_value),
    stringsAsFactors = FALSE
  ))
  for (type in c("component", "prior", "posterior")) {
    value <- switch(
      type,
      component = component_mean,
      prior = pieces$prior,
      posterior = posterior
    )
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset,
      row_id = rep(seq_len(nrow(data)), times = ncol(value)),
      prediction_type = type,
      component = rep(component_names, each = nrow(data)),
      value = as.vector(value),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

write_csv <- function(x, path) {
  utils::write.table(
    x, path, sep = ",", row.names = FALSE, col.names = TRUE,
    quote = TRUE, na = "NA"
  )
}

sha256 <- function(path) {
  output <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
}

case_specs <- list(
  gaussian = list(
    family = "gaussian",
    formula = y ~ z1 | x1,
    n = 180L,
    heldout_n = 80L,
    simulation_seed = 61001L,
    heldout_seed = 61101L,
    separation = 2,
    trials = NA_integer_,
    mean_count = NA_real_,
    fit_seed = 62001L,
    bandwidth = 0.45,
    u_grid = c(0.25, 0.50, 0.75),
    n_starts = 3L,
    maxit = 200L
  ),
  binomial = list(
    family = "binomial",
    formula = cbind(success, failure) ~ z1 | x1,
    n = 180L,
    heldout_n = 80L,
    simulation_seed = 61002L,
    heldout_seed = 61102L,
    separation = 1.8,
    trials = 12L,
    mean_count = NA_real_,
    fit_seed = 62002L,
    bandwidth = 0.45,
    u_grid = c(0.25, 0.50, 0.75),
    n_starts = 3L,
    maxit = 200L
  ),
  negative_binomial = list(
    family = "negative-binomial",
    formula = y ~ z1 + offset(log_size_factor) | x1,
    n = 360L,
    heldout_n = 120L,
    simulation_seed = 61013L,
    heldout_seed = 61113L,
    separation = 1.8,
    trials = NA_integer_,
    mean_count = 8,
    fit_seed = 62013L,
    bandwidth = 0.55,
    u_grid = c(0.30, 0.50, 0.70),
    n_starts = 2L,
    maxit = 150L
  )
)

manifest_rows <- list()
for (case in names(case_specs)) {
  spec <- case_specs[[case]]
  simulate_case <- function(n, seed) {
    switch(
      case,
      gaussian = simulate_vcmoe_gaussian(
        n = n, k = 2L, seed = seed, separation = spec$separation
      ),
      binomial = simulate_vcmoe_binomial(
        n = n, k = 2L, seed = seed, separation = spec$separation,
        trials = spec$trials
      ),
      negative_binomial = simulate_vcmoe_negbin(
        n = n, k = 2L, seed = seed, separation = spec$separation,
        mean_count = spec$mean_count
      )
    )$data
  }
  training <- simulate_case(spec$n, spec$simulation_seed)
  heldout <- simulate_case(spec$heldout_n, spec$heldout_seed)
  fit <- suppressWarnings(vcmoe_fit(
    formula = spec$formula,
    data = training,
    u = "u",
    k = 2L,
    family = spec$family,
    bandwidth = spec$bandwidth,
    u_grid = spec$u_grid,
    label = "global",
    control = list(
      maxit = spec$maxit,
      n_starts = spec$n_starts,
      seed = spec$fit_seed,
      warn_ambiguous = FALSE
    )
  ))
  if (!all(fit$diagnostics$converged)) {
    stop(case, " fixture fit did not converge at every grid point.", call. = FALSE)
  }
  if (any(fit$diagnostics$ambiguous)) {
    stop(case, " fixture fit has ambiguous label alignment.", call. = FALSE)
  }
  if (!is.null(fit$coefficients$theta) &&
      any(fit$coefficients$theta >= 0.99 * fit$control$negbin_theta_max)) {
    stop(case, " fixture fit has boundary-dominated theta.", call. = FALSE)
  }

  training_path <- file.path(fixture_dir, paste0(case, "-training.csv"))
  heldout_path <- file.path(fixture_dir, paste0(case, "-heldout.csv"))
  coefficient_path <- file.path(fixture_dir, paste0(case, "-coefficients.csv"))
  prediction_path <- file.path(fixture_dir, paste0(case, "-predictions.csv"))
  diagnostic_path <- file.path(fixture_dir, paste0(case, "-diagnostics.csv"))
  write_csv(training, training_path)
  write_csv(heldout, heldout_path)

  coefficients <- do.call(rbind, list(
    array_long(fit$coefficients$expert, case, "expert"),
    array_long(fit$coefficients$expert_slope, case, "expert_slope"),
    array_long(fit$coefficients$gating, case, "gating"),
    array_long(fit$coefficients$gating_slope, case, "gating_slope"),
    array_long(fit$coefficients$sigma, case, "sigma", "log_sigma_scale"),
    array_long(fit$coefficients$sigma_slope, case, "sigma_slope", "log_sigma_slope"),
    array_long(fit$coefficients$theta, case, "theta", "theta")
  ))
  write_csv(coefficients, coefficient_path)
  predictions <- rbind(
    prediction_long(fit, training, "training"),
    prediction_long(fit, heldout, "heldout")
  )
  predictions$case <- case
  predictions <- predictions[, c(
    "case", "dataset", "row_id", "prediction_type", "component", "value"
  )]
  write_csv(predictions, prediction_path)
  diagnostics <- data.frame(
    case = case,
    grid_id = seq_along(fit$u_grid),
    u = fit$u_grid,
    loglik = fit$diagnostics$loglik,
    converged = fit$diagnostics$converged,
    iterations = fit$diagnostics$iterations,
    selected_start = fit$diagnostics$selected_start,
    permutation = apply(fit$diagnostics$permutations, 1L, paste, collapse = ","),
    stringsAsFactors = FALSE
  )
  write_csv(diagnostics, diagnostic_path)

  manifest_rows[[case]] <- data.frame(
    case = case,
    package_version = "0.1.0",
    source_tag = "v0.1.0",
    source_commit = source_commit,
    family = spec$family,
    formula = paste(deparse(spec$formula), collapse = ""),
    n = spec$n,
    heldout_n = spec$heldout_n,
    simulation_seed = spec$simulation_seed,
    heldout_seed = spec$heldout_seed,
    fit_seed = spec$fit_seed,
    separation = spec$separation,
    trials = spec$trials,
    mean_count = spec$mean_count,
    bandwidth = spec$bandwidth,
    u_grid = paste(spec$u_grid, collapse = ";"),
    n_starts = spec$n_starts,
    maxit = spec$maxit,
    R_version = R.version.string,
    training_sha256 = sha256(training_path),
    heldout_sha256 = sha256(heldout_path),
    coefficients_sha256 = sha256(coefficient_path),
    predictions_sha256 = sha256(prediction_path),
    diagnostics_sha256 = sha256(diagnostic_path),
    stringsAsFactors = FALSE
  )
}

write_csv(
  do.call(rbind, manifest_rows),
  file.path(fixture_dir, "manifest.csv")
)
message("Wrote element-level v0.1.0 regression fixtures to ", fixture_dir)
