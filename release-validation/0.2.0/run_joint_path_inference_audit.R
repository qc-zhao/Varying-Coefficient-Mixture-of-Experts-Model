#!/usr/bin/env Rscript

# Lightweight structural audit for engine-aware inference in VCMoE 0.2.0.
# Run from the package root.

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
description <- read.dcf(file.path(root, "DESCRIPTION"), fields = c("Package", "Version"))
if (!identical(unname(description[1L, "Package"]), "VCMoE")) {
  stop("Run this script from the VCMoE package root.", call. = FALSE)
}
release_version <- unname(description[1L, "Version"])

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package `pkgload` is required.", call. = FALSE)
}
pkgload::load_all(root, reset = TRUE, recompile = FALSE, export_all = FALSE,
                  helpers = FALSE, quiet = TRUE)

results_dir <- file.path(root, "release-validation", "0.2.0", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
replicate_path <- file.path(results_dir, "joint_path_inference_audit.csv")
summary_path <- file.path(results_dir, "joint_path_inference_audit_summary.csv")
plot_path <- file.path(results_dir, "joint_path_inference_audit.png")

scenarios <- list(
  list(id = "gaussian", family = "gaussian", seed = 9201L, n = 90L,
       maxit = 60L),
  list(id = "negbin_offset", family = "negative-binomial", seed = 9211L,
       n = 90L, maxit = 80L)
)
engines <- c("local_grid_em", "joint_path_em")

simulate_case <- function(scenario) {
  if (identical(scenario$family, "gaussian")) {
    VCMoE::simulate_vcmoe_gaussian(
      n = scenario$n, k = 2L, seed = scenario$seed, separation = 1.8
    )
  } else {
    VCMoE::simulate_vcmoe_negbin(
      n = scenario$n, k = 2L, seed = scenario$seed, separation = 1.8,
      mean_count = 8
    )
  }
}

fit_formula <- function(family) {
  if (identical(family, "negative-binomial")) {
    y ~ z1 + offset(log_size_factor) | x1
  } else {
    y ~ z1 | x1
  }
}

capture_warnings <- function(expr) {
  warnings <- character(0L)
  value <- withCallingHandlers(
    expr,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

empty_row <- function(scenario, engine) {
  data.frame(
    release_version = release_version,
    scenario = scenario$id,
    family = scenario$family,
    engine = engine,
    fit_success = FALSE,
    fit_converged = FALSE,
    finite_predictions = FALSE,
    posterior_row_sum_max_error = NA_real_,
    glrt_status = "not_run",
    glrt_null_converged = FALSE,
    glrt_null_engine = NA_character_,
    statistic_criterion = NA_character_,
    paper_criterion_match = NA_character_,
    null_averaging_measure = NA_character_,
    constrained_mle = NA,
    calibration = NA_character_,
    calibration_status = NA_character_,
    theory_status = NA_character_,
    confband_success = FALSE,
    confband_finite_fraction = NA_real_,
    covariance_target = NA_character_,
    estimator_covariance_match = NA_character_,
    coverage_theory = NA_character_,
    bootstrap_success = FALSE,
    bootstrap_engine_preserved = FALSE,
    warning_count = 0L,
    warnings = "",
    error = "",
    stringsAsFactors = FALSE
  )
}

run_case <- function(scenario, engine) {
  row <- empty_row(scenario, engine)
  evaluated <- tryCatch({
    simulation <- simulate_case(scenario)
    fitted <- capture_warnings(VCMoE::vcmoe_fit(
      fit_formula(scenario$family),
      data = simulation$data,
      u = "u",
      k = 2L,
      family = scenario$family,
      bandwidth = 0.50,
      u_grid = c(0.25, 0.50, 0.75),
      engine = engine,
      control = list(
        maxit = scenario$maxit,
        n_starts = 1L,
        seed = scenario$seed + 1L,
        warn_ambiguous = FALSE
      )
    ))
    fit <- fitted$value
    row$fit_success <- TRUE
    row$fit_converged <- if (identical(engine, "joint_path_em")) {
      isTRUE(fit$diagnostics$joint_path_converged)
    } else {
      all(fit$diagnostics$converged %in% TRUE)
    }
    prediction <- predict(fit, type = "mean")
    posterior <- predict(fit, type = "posterior")
    row$finite_predictions <- all(is.finite(prediction)) && all(is.finite(posterior))
    row$posterior_row_sum_max_error <- max(abs(rowSums(posterior) - 1))

    banded <- capture_warnings(VCMoE::vcmoe_confband(
      fit,
      level = 0.90,
      coefficient_set = if (identical(scenario$family, "gaussian")) {
        c("expert", "gating", "sigma")
      } else {
        c("expert", "gating", "theta")
      },
      strict = FALSE
    ))
    band <- banded$value
    row$confband_success <- inherits(band, "vcmoe_confband")
    row$confband_finite_fraction <- mean(
      is.finite(band$intervals$estimate) & is.finite(band$intervals$se)
    )
    row$covariance_target <- band$settings$covariance_target
    row$estimator_covariance_match <- band$settings$estimator_covariance_match
    row$coverage_theory <- band$settings$coverage_theory

    tested <- capture_warnings(VCMoE::vcmoe_glrt(
      fit,
      simulation$data,
      test = "coefficient",
      coefficient_set = "gating",
      component = 1L,
      term = "x1",
      control = list(maxit = 120L, reltol = 1e-5, strict = FALSE)
    ))
    glrt <- tested$value
    row$glrt_status <- glrt$status
    row$glrt_null_converged <- is.finite(glrt$null_convergence) &&
      glrt$null_convergence == 0L
    row$glrt_null_engine <- glrt$settings$null_engine_id
    row$statistic_criterion <- glrt$settings$statistic_criterion
    row$paper_criterion_match <- glrt$settings$paper_criterion_match
    row$null_averaging_measure <- glrt$settings$null_averaging_measure
    row$constrained_mle <- glrt$settings$constrained_mle
    row$calibration <- glrt$settings$calibration
    row$calibration_status <- glrt$settings$calibration_status
    row$theory_status <- glrt$settings$theory_status

    booted <- capture_warnings(VCMoE::vcmoe_bootstrap(
      fit,
      simulation$data,
      B = 2L,
      seed = scenario$seed + 2L,
      min_successful = 2L,
      keep_fits = TRUE,
      control = list(
        maxit = 12L,
        n_starts = 1L,
        warn_ambiguous = FALSE
      )
    ))
    bootstrap <- booted$value
    row$bootstrap_success <- inherits(bootstrap, "vcmoe_bootstrap") &&
      length(bootstrap$fits) == 2L
    row$bootstrap_engine_preserved <- row$bootstrap_success &&
      all(vapply(
        bootstrap$fits,
        function(refit) identical(refit$engine_id, engine),
        logical(1L)
      ))
    all_warnings <- unique(c(fitted$warnings, banded$warnings,
                             tested$warnings, booted$warnings))
    row$warning_count <- length(all_warnings)
    row$warnings <- paste(all_warnings, collapse = " | ")
    row
  }, error = function(error) {
    row$error <- conditionMessage(error)
    row
  })
  evaluated
}

rows <- list()
row_id <- 0L
for (scenario in scenarios) {
  for (engine in engines) {
    row_id <- row_id + 1L
    message(sprintf("Inference audit %s / %s", scenario$id, engine))
    rows[[row_id]] <- run_case(scenario, engine)
  }
}
results <- do.call(rbind, rows)

results$structural_pass <- with(results,
  fit_success & finite_predictions & posterior_row_sum_max_error < 1e-8 &
    confband_success & confband_finite_fraction > 0 & bootstrap_success &
    bootstrap_engine_preserved & calibration == "none" &
    calibration_status %in% c("not_requested", "blocked") &
    (glrt_null_converged | glrt_status == "blocked")
)
joint <- results$engine == "joint_path_em"
results$metadata_pass <- TRUE
results$metadata_pass[joint] <- with(results[joint, ],
  statistic_criterion == "nearest_grid_sample_loglik" &
    paper_criterion_match == "nearest_grid_approximation" &
    null_averaging_measure == "nearest_grid_assignment_frequency" &
    !constrained_mle &
    estimator_covariance_match == "asymptotic_plugin_not_exact_finite_grid_match"
)
results$metadata_pass[!joint] <- with(results[!joint, ],
  statistic_criterion == "summed_kernel_weighted_local_pseudologlik" &
    paper_criterion_match == "no" & constrained_mle
)
results$pass <- results$structural_pass & results$metadata_pass

summary <- data.frame(
  release_version = release_version,
  cases = nrow(results),
  cases_passed = sum(results$pass),
  fit_successes = sum(results$fit_success),
  finite_prediction_cases = sum(results$finite_predictions),
  confband_successes = sum(results$confband_success),
  bootstrap_engine_preserved_cases = sum(results$bootstrap_engine_preserved),
  converged_nulls = sum(results$glrt_null_converged),
  correctly_blocked_nonconverged_nulls = sum(
    !results$glrt_null_converged & results$glrt_status == "blocked"
  ),
  release_gate = if (all(results$pass)) "PASS" else "FAIL",
  stringsAsFactors = FALSE
)

utils::write.csv(results, replicate_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)

grDevices::png(plot_path, width = 1400, height = 800, res = 150)
old_par <- graphics::par(mar = c(7, 4, 3, 1))
on.exit({graphics::par(old_par); grDevices::dev.off()}, add = TRUE)
labels <- paste(results$scenario, results$engine, sep = "\n")
graphics::barplot(
  as.integer(results$pass), names.arg = labels, las = 2, ylim = c(0, 1.1),
  col = ifelse(results$pass, "#2F855A", "#C53030"),
  ylab = "Release audit pass", main = "VCMoE 0.2.0 inference consistency audit"
)
graphics::abline(h = 1, lty = 3, col = "grey40")

print(summary)
if (!identical(summary$release_gate, "PASS")) {
  stop("Joint-path inference audit failed; inspect the result CSV.", call. = FALSE)
}
