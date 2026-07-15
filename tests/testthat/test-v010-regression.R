.v010_fixture_dir <- function() {
  testthat::test_path("fixtures", "v010")
}

.v010_case_spec <- function(case) {
  switch(
    case,
    gaussian = list(
      formula = y ~ z1 | x1,
      family = "gaussian",
      bandwidth = 0.45,
      u_grid = c(0.25, 0.50, 0.75),
      maxit = 200L,
      n_starts = 3L,
      seed = 62001L,
      tolerance = 2e-6
    ),
    binomial = list(
      formula = cbind(success, failure) ~ z1 | x1,
      family = "binomial",
      bandwidth = 0.45,
      u_grid = c(0.25, 0.50, 0.75),
      maxit = 200L,
      n_starts = 3L,
      seed = 62002L,
      tolerance = 2e-6
    ),
    negative_binomial = list(
      formula = y ~ z1 + offset(log_size_factor) | x1,
      family = "negative-binomial",
      bandwidth = 0.55,
      u_grid = c(0.30, 0.50, 0.70),
      maxit = 150L,
      n_starts = 2L,
      seed = 62013L,
      tolerance = 5e-6
    )
  )
}

.v010_fit_case <- function(case, engine = NULL) {
  spec <- .v010_case_spec(case)
  data <- utils::read.csv(
    file.path(.v010_fixture_dir(), paste0(case, "-training.csv")),
    check.names = FALSE
  )
  args <- list(
    formula = spec$formula,
    data = data,
    u = "u",
    k = 2L,
    family = spec$family,
    bandwidth = spec$bandwidth,
    u_grid = spec$u_grid,
    label = "global",
    control = list(
      maxit = spec$maxit,
      n_starts = spec$n_starts,
      seed = spec$seed,
      warn_ambiguous = FALSE
    )
  )
  if (!is.null(engine)) {
    args$engine <- engine
  }
  suppressWarnings(do.call(vcmoe_fit, args))
}

.relative_max_error <- function(actual, expected) {
  max(abs(actual - expected) / pmax(1, abs(expected)))
}

.actual_prediction_long <- function(fit, data, dataset) {
  pieces <- VCMoE:::.predict_components(fit, newdata = data)
  posterior <- predict(fit, newdata = data, type = "posterior")
  component_mean <- predict(fit, newdata = data, type = "component")
  rows <- list(data.frame(
    dataset = dataset,
    row_id = seq_len(nrow(data)),
    prediction_type = "mean",
    component = NA_character_,
    value = as.numeric(predict(fit, newdata = data, type = "mean")),
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
      component = rep(colnames(pieces$prior), each = nrow(data)),
      value = as.vector(value),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

test_that("default Gaussian, Binomial, and NB fits preserve the 0.1.0 numeric contract", {
  for (case in c("gaussian", "binomial", "negative_binomial")) {
    spec <- .v010_case_spec(case)
    fit <- .v010_fit_case(case)
    expect_identical(fit$engine_id %||% "local_grid_em", "local_grid_em")
    expect_identical(fit$family, spec$family)
    expect_identical(fit$k, 2L)
    expect_equal(fit$u_grid, spec$u_grid, tolerance = 0)

    expected_coef <- utils::read.csv(
      file.path(.v010_fixture_dir(), paste0(case, "-coefficients.csv")),
      check.names = FALSE
    )
    for (set in unique(expected_coef$coefficient_set)) {
      expected <- expected_coef$estimate[expected_coef$coefficient_set == set]
      actual <- as.vector(fit$coefficients[[set]])
      tolerance <- if (identical(set, "theta")) 2e-5 else spec$tolerance
      expect_lte(
        .relative_max_error(actual, expected),
        tolerance
      )
    }

    expected_diag <- utils::read.csv(
      file.path(.v010_fixture_dir(), paste0(case, "-diagnostics.csv")),
      check.names = FALSE
    )
    expect_lte(
      .relative_max_error(fit$diagnostics$loglik, expected_diag$loglik),
      spec$tolerance
    )
    expect_identical(fit$diagnostics$converged, expected_diag$converged)
    expect_identical(fit$diagnostics$selected_start, expected_diag$selected_start)
    expect_identical(
      apply(fit$diagnostics$permutations, 1L, paste, collapse = ","),
      expected_diag$permutation
    )

    training <- utils::read.csv(
      file.path(.v010_fixture_dir(), paste0(case, "-training.csv")),
      check.names = FALSE
    )
    heldout <- utils::read.csv(
      file.path(.v010_fixture_dir(), paste0(case, "-heldout.csv")),
      check.names = FALSE
    )
    actual_prediction <- rbind(
      .actual_prediction_long(fit, training, "training"),
      .actual_prediction_long(fit, heldout, "heldout")
    )
    expected_prediction <- utils::read.csv(
      file.path(.v010_fixture_dir(), paste0(case, "-predictions.csv")),
      check.names = FALSE,
      na.strings = "NA"
    )
    expect_identical(
      actual_prediction[, c("dataset", "row_id", "prediction_type", "component")],
      expected_prediction[, c("dataset", "row_id", "prediction_type", "component")]
    )
    expect_lte(
      .relative_max_error(actual_prediction$value, expected_prediction$value),
      spec$tolerance
    )
    posterior <- predict(fit, type = "posterior")
    prior <- VCMoE:::.predict_components(fit, newdata = training)$prior
    expect_equal(rowSums(posterior), rep(1, nrow(training)), tolerance = 1e-10)
    expect_equal(rowSums(prior), rep(1, nrow(training)), tolerance = 1e-10)
  }
})

test_that("omitting engine remains exactly equivalent to explicit local-grid EM", {
  default_fit <- .v010_fit_case("gaussian")
  explicit_fit <- .v010_fit_case("gaussian", engine = "local_grid_em")
  expect_equal(default_fit$coefficients, explicit_fit$coefficients, tolerance = 0)
  expect_equal(
    predict(default_fit, type = "mean"),
    predict(explicit_fit, type = "mean"),
    tolerance = 0
  )
  expect_equal(
    predict(default_fit, type = "posterior"),
    predict(explicit_fit, type = "posterior"),
    tolerance = 0
  )
})
