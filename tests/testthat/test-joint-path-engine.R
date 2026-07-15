test_that("joint-path Gaussian fits preserve the public vcmoe ABI for k = 2:10", {
  for (k in 2L:10L) {
    sim <- simulate_vcmoe_gaussian(
      n = max(70L, 12L * k),
      k = k,
      seed = 8100L + k,
      separation = 1.7
    )
    fit <- suppressWarnings(vcmoe_fit(
      y ~ 1 | 1,
      data = sim$data,
      u = "u",
      k = k,
      bandwidth = 0.55,
      u_grid = 0.5,
      engine = "joint_path_em",
      control = list(
        maxit = 2L,
        n_starts = 1L,
        seed = 8200L + k,
        warn_ambiguous = FALSE
      )
    ))
    expect_s3_class(fit, "vcmoe")
    expect_identical(fit$engine_id, "joint_path_em")
    expect_identical(fit$diagnostics$engine, "joint_path_em")
    expect_equal(dim(coef(fit, "expert")), c(1L, k, 1L))
    expect_equal(dim(predict(fit, type = "posterior")), c(nrow(sim$data), k))
    expect_equal(
      rowSums(predict(fit, type = "posterior")),
      rep(1, nrow(sim$data)),
      tolerance = 1e-8
    )
    expect_true(all(is.finite(predict(fit, type = "mean"))))
  }
})

test_that("joint-path Binomial and NB offset fits return finite family outputs", {
  grouped <- simulate_vcmoe_binomial(
    n = 90, k = 2, seed = 8301, separation = 1.5, trials = 8
  )
  grouped_fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ z1 | x1,
    data = grouped$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.45,
    u_grid = c(0.3, 0.7),
    engine = "joint_path_em",
    control = list(maxit = 6L, n_starts = 1L, seed = 8302, warn_ambiguous = FALSE)
  ))
  expect_identical(grouped_fit$engine_id, "joint_path_em")
  expect_true(all(predict(grouped_fit, type = "mean") >= 0))
  expect_true(all(predict(grouped_fit, type = "mean") <= 1))

  nb <- simulate_vcmoe_negbin(n = 90, k = 2, seed = 8401, separation = 1.5)
  nb_fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = nb$data,
    u = "u",
    family = "negative-binomial",
    bandwidth = 0.45,
    u_grid = c(0.3, 0.7),
    engine = "joint_path_em",
    control = list(maxit = 6L, n_starts = 1L, seed = 8402, warn_ambiguous = FALSE)
  ))
  expect_identical(nb_fit$engine_id, "joint_path_em")
  expect_true(all(is.finite(coef(nb_fit, "theta"))))
  expect_true(all(coef(nb_fit, "theta") > 0))
  expect_true(all(is.finite(predict(nb_fit, type = "mean"))))
})

test_that("joint-path dense-grid and progress guardrails are active", {
  sim <- simulate_vcmoe_gaussian(n = 30, k = 2, seed = 8501)
  expect_error(
    vcmoe_fit(
      y ~ z1 | x1,
      data = sim$data,
      u = "u",
      bandwidth = 0.4,
      u_grid = sort(unique(sim$data$u)),
      engine = "joint_path_em",
      control = list(maxit = 1L, n_starts = 1L, seed = 8502)
    ),
    "dense-grid guard"
  )

  progress <- tempfile(fileext = ".csv")
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.4,
    u_grid = c(0.3, 0.7),
    engine = "joint_path_em",
    progress = progress,
    control = list(maxit = 2L, n_starts = 1L, seed = 8503, warn_ambiguous = FALSE)
  ))
  expect_true(file.exists(progress))
  progress_rows <- utils::read.csv(progress)
  expect_true(all(c("fit_started", "joint_path_iteration", "fit_finished") %in%
    progress_rows$event))
  expect_equal(sum(fit$diagnostics$joint_path_assignment$n_assigned), nrow(sim$data))
})

test_that("joint-path reports a maxit hit explicitly", {
  sim <- simulate_vcmoe_gaussian(n = 45, k = 2, seed = 8551)
  warnings <- character(0L)
  fit <- withCallingHandlers(
    vcmoe_fit(
      y ~ z1 | x1,
      data = sim$data,
      u = "u",
      bandwidth = 0.45,
      u_grid = c(0.3, 0.7),
      engine = "joint_path_em",
      control = list(
        maxit = 1L,
        n_starts = 1L,
        seed = 8552,
        warn_ambiguous = FALSE
      )
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_false(fit$diagnostics$joint_path_converged)
  expect_true(any(grepl("did not converge", warnings, fixed = TRUE)))
})

test_that("joint-path analytic and bootstrap inference preserve the engine", {
  sim <- simulate_vcmoe_gaussian(n = 75, k = 2, seed = 8601, separation = 1.7)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = c(0.3, 0.7),
    engine = "joint_path_em",
    control = list(maxit = 8L, n_starts = 1L, seed = 8602, warn_ambiguous = FALSE)
  ))
  band <- NULL
  expect_warning(
    band <- vcmoe_confband(
      fit,
      level = 0.90,
      coefficient_set = c("expert", "gating", "sigma"),
      strict = FALSE
    ),
    "cross-grid responsibility uncertainty"
  )
  expect_s3_class(band, "vcmoe_confband")
  expect_identical(band$settings$engine_id, "joint_path_em")
  expect_identical(
    band$settings$covariance_target,
    "observed_local_marginal_likelihood"
  )
  expect_identical(
    band$settings$estimator_covariance_match,
    "asymptotic_plugin_not_exact_finite_grid_match"
  )
  expect_false(band$settings$shared_path_uncertainty_accounted)
  expect_false(band$settings$label_uncertainty_accounted)
  expect_identical(
    band$settings$coverage_theory,
    "local_likelihood_asymptotic_plugin_not_finite_grid_joint_path_theory"
  )
  expect_true(all(c(
    "covariance_target",
    "estimator_covariance_match",
    "shared_path_uncertainty_accounted",
    "label_uncertainty_accounted",
    "coverage_theory"
  ) %in% names(band$intervals)))
  expect_identical(
    band$settings$estimating_equation,
    "jasa_observed_local_likelihood_plugin"
  )
  expect_identical(
    band$settings$covariance_scope,
    "local_asymptotic_plugin_excludes_finite_grid_cross_grid_coupling"
  )
  expect_true(all(c("score_sum_max", "score_imbalance_max") %in%
    names(band$diagnostics)))
  expect_true(all(band$intervals$status %in% c("ok", "warning", "blocked")))
  expect_gt(nrow(band$intervals), 0L)

  boot <- suppressWarnings(vcmoe_bootstrap(
    fit,
    sim$data,
    B = 2L,
    seed = 8603,
    min_successful = 2L,
    keep_fits = TRUE,
    control = list(maxit = 3L, n_starts = 1L, warn_ambiguous = FALSE)
  ))
  expect_identical(boot$settings$engine_id, "joint_path_em")
  expect_true(all(vapply(
    boot$fits,
    function(x) identical(x$engine_id, "joint_path_em"),
    logical(1L)
  )))
})

test_that("joint-path reduced and GLRT null fits use projected EM", {
  sim <- simulate_vcmoe_gaussian(n = 80, k = 2, seed = 8701, separation = 1.7)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = c(0.25, 0.5, 0.75),
    engine = "joint_path_em",
    control = list(maxit = 8L, n_starts = 1L, seed = 8702, warn_ambiguous = FALSE)
  ))
  reduced <- suppressWarnings(vcmoe_fit_reduced(
    fit,
    constrain = "gating_constant",
    control = list(maxit = 4L, strict = FALSE)
  ))
  expect_identical(
    reduced$diagnostics$glrt_null_engine,
    "joint_path_em_constrained_null"
  )
  expect_true(all(is.finite(predict(reduced, type = "mean"))))
  expect_true(all(apply(reduced$coefficients$gating, c(2L, 3L), stats::sd) < 1e-10))

  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    calibration = "none",
    control = list(maxit = 4L, strict = FALSE)
  ))
  expect_identical(
    glrt$null_fit$diagnostics$glrt_null_engine,
    "joint_path_em_constrained_null"
  )
  expect_identical(
    glrt$null_fit$diagnostics$glrt_statistic_criterion,
    "nearest_grid_sample_loglik"
  )
  expect_identical(
    glrt$settings$null_averaging_measure,
    "nearest_grid_assignment_frequency"
  )
  expect_identical(glrt$settings$paper_criterion_match, "nearest_grid_approximation")
  expect_false(glrt$settings$constrained_mle)
  expect_identical(
    glrt$settings$theory_status,
    "nearest_grid_approximation_not_exact_manuscript_theory"
  )
  expect_equal(
    glrt$full_loglik,
    VCMoE:::.vcmoe_glrt_global_pseudologlik(fit),
    tolerance = 1e-8
  )
  expect_equal(
    utils::tail(fit$diagnostics$joint_path_trace$objective, 1L),
    glrt$full_loglik,
    tolerance = 1e-8
  )
  expect_true(is.data.frame(glrt$null_fit$diagnostics$joint_path_trace))
  expect_gt(nrow(glrt$null_fit$diagnostics$joint_path_trace), 0L)
})

test_that("joint-path projected null uses observation-assignment weights", {
  grid_parameters <- list(
    c(shared = 1, fixed = 7),
    c(shared = 4, fixed = 8),
    c(shared = 10, fixed = 9)
  )
  constraint <- list(
    shared_parameters = "shared",
    fixed_zero_parameters = "fixed"
  )
  projected <- VCMoE:::.vcmoe_glrt_project_grid_parameters(
    grid_parameters,
    constraint,
    projection_weights = c(8, 2, 0)
  )
  expect_equal(
    vapply(projected$parameters, `[[`, numeric(1L), "shared"),
    rep(1.6, 3L)
  )
  expect_equal(
    vapply(projected$parameters, `[[`, numeric(1L), "fixed"),
    rep(0, 3L)
  )
  expect_identical(
    projected$projection$averaging_measure,
    "nearest_grid_assignment_frequency"
  )
  expect_equal(projected$projection$weight_sum, 10)
  expect_equal(projected$projection$nonzero_grid_points, 2L)
  expect_error(
    VCMoE:::.vcmoe_glrt_project_grid_parameters(
      grid_parameters,
      constraint,
      projection_weights = c(0, 0, 0)
    ),
    "positive sum"
  )
})

test_that("joint-path null failures propagate and nonconvergence blocks GLRT", {
  sim <- simulate_vcmoe_gaussian(n = 80, k = 2, seed = 8721, separation = 1.7)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = c(0.25, 0.5, 0.75),
    engine = "joint_path_em",
    control = list(maxit = 8L, n_starts = 1L, seed = 8722, warn_ambiguous = FALSE)
  ))
  constraint <- VCMoE:::.vcmoe_glrt_constraint(
    fit,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1"
  )
  broken <- fit
  broken$fitted <- NULL
  expect_error(
    VCMoE:::.vcmoe_glrt_null_fit(broken, constraint),
    "fit\\$fitted"
  )

  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    calibration = "none",
    control = list(maxit = 1L, strict = FALSE)
  ))
  expect_identical(glrt$status, "blocked")
  expect_match(glrt$block_reason, "null_optimizer_not_fully_converged")
  expect_error(
    vcmoe_fit_reduced(
      fit,
      constrain = "expert_constant",
      control = list(maxit = 1L, strict = TRUE)
    ),
    "did not converge"
  )
})

test_that("local-grid GLRT retains the constrained BFGS null engine", {
  sim <- simulate_vcmoe_gaussian(n = 55, k = 2, seed = 8751, separation = 1.5)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 5L, n_starts = 1L, seed = 8752, warn_ambiguous = FALSE)
  ))
  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    calibration = "none",
    control = list(maxit = 3L, strict = FALSE)
  ))
  expect_identical(
    glrt$null_fit$diagnostics$glrt_null_engine,
    "local_grid_bfgs_constrained_null"
  )
  expect_identical(
    glrt$settings$statistic_criterion,
    "summed_kernel_weighted_local_pseudologlik"
  )
  expect_identical(glrt$settings$paper_criterion_match, "no")
  expect_true(glrt$settings$constrained_mle)
})

test_that("GLRT defaults to an uncalibrated statistic", {
  sim <- simulate_vcmoe_gaussian(n = 55, k = 2, seed = 8761, separation = 1.5)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.45,
    u_grid = c(0.35, 0.65),
    control = list(maxit = 20L, n_starts = 1L, seed = 8762, warn_ambiguous = FALSE)
  ))
  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    control = list(maxit = 100L, strict = FALSE)
  ))
  expect_identical(glrt$settings$calibration, "none")
  expect_identical(
    glrt$settings$calibration_status,
    if (identical(glrt$status, "blocked")) "blocked" else "not_requested"
  )
  expect_true(is.na(glrt$p_value))
})

test_that("joint-path GLRT bootstrap preserves full and null engines", {
  sim <- simulate_vcmoe_gaussian(n = 60, k = 2, seed = 8771, separation = 2.0)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth = 0.60,
    u_grid = 0.50,
    engine = "joint_path_em",
    control = list(maxit = 60L, n_starts = 1L, seed = 8772, warn_ambiguous = FALSE)
  ))
  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    calibration = "bootstrap",
    B = 2L,
    seed = 8773,
    control = list(maxit = 200L, reltol = 1e-5, strict = FALSE),
    refit_control = list(
      maxit = 60L,
      n_starts = 1L,
      warn_ambiguous = FALSE
    )
  ))
  expect_false(identical(glrt$status, "blocked"))
  expect_identical(glrt$settings$calibration_status, "empirical_parametric_bootstrap")
  expect_equal(nrow(glrt$replicate_summary), 2L)
  expect_true(all(glrt$replicate_summary$status == "ok"))
  expect_true(all(glrt$replicate_summary$engine_id == "joint_path_em"))
  expect_true(all(
    glrt$replicate_summary$null_engine_id == "joint_path_em_constrained_null"
  ))
  expect_true(all(glrt$replicate_summary$null_convergence == 0L))
})

test_that("bandwidth selection uses one engine consistently", {
  sim <- simulate_vcmoe_gaussian(n = 60, k = 2, seed = 8801, separation = 1.5)
  selection <- suppressWarnings(vcmoe_select_bandwidth(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    bandwidth_grid = 0.45,
    folds = 2L,
    u_grid = c(0.35, 0.65),
    engine = "joint_path_em",
    control = list(maxit = 3L, n_starts = 1L, seed = 8802, warn_ambiguous = FALSE),
    seed = 8803
  ))
  expect_identical(selection$settings$engine, "joint_path_em")
  expect_true(all(selection$cv_details$engine_id == "joint_path_em"))
  expect_identical(selection$fit$engine_id, "joint_path_em")
})

test_that("Binomial and NB bootstrap refits preserve joint-path EM", {
  grouped <- simulate_vcmoe_binomial(
    n = 80, k = 2, seed = 8901, separation = 1.8, trials = 8
  )
  grouped_fit <- suppressWarnings(vcmoe_fit(
    cbind(success, failure) ~ z1 | x1,
    data = grouped$data,
    u = "u",
    family = "binomial",
    bandwidth = 0.50,
    u_grid = 0.50,
    engine = "joint_path_em",
    control = list(maxit = 6L, n_starts = 1L, seed = 8902, warn_ambiguous = FALSE)
  ))
  grouped_boot <- suppressWarnings(vcmoe_bootstrap(
    grouped_fit,
    grouped$data,
    B = 2L,
    seed = 8903,
    min_successful = 2L,
    keep_fits = TRUE,
    control = list(maxit = 3L, n_starts = 1L, warn_ambiguous = FALSE)
  ))
  expect_length(grouped_boot$fits, 2L)
  expect_true(all(vapply(
    grouped_boot$fits,
    function(x) identical(x$engine_id, "joint_path_em"),
    logical(1L)
  )))

  nb <- simulate_vcmoe_negbin(
    n = 90, k = 2, seed = 8911, separation = 1.8, mean_count = 8
  )
  nb_fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 + offset(log_size_factor) | x1,
    data = nb$data,
    u = "u",
    family = "negative-binomial",
    bandwidth = 0.55,
    u_grid = 0.50,
    engine = "joint_path_em",
    control = list(maxit = 6L, n_starts = 1L, seed = 8912, warn_ambiguous = FALSE)
  ))
  nb_band <- suppressWarnings(vcmoe_confband(
    nb_fit,
    level = 0.90,
    coefficient_set = c("expert", "gating", "theta"),
    strict = FALSE
  ))
  expect_identical(nb_band$settings$engine_id, "joint_path_em")
  expect_true(any(
    nb_band$intervals$coefficient_set == "nuisance" &
      nb_band$intervals$term == "log_theta"
  ))

  nb_boot <- suppressWarnings(vcmoe_bootstrap(
    nb_fit,
    nb$data,
    B = 2L,
    seed = 8913,
    min_successful = 2L,
    keep_fits = TRUE,
    control = list(maxit = 3L, n_starts = 1L, warn_ambiguous = FALSE)
  ))
  expect_length(nb_boot$fits, 2L)
  expect_true(all(vapply(
    nb_boot$fits,
    function(x) identical(x$engine_id, "joint_path_em"),
    logical(1L)
  )))
})

test_that("K=3 joint-path inference uses the common fitted-object ABI", {
  sim <- simulate_vcmoe_gaussian(n = 100, k = 3, seed = 8921, separation = 1.8)
  fit <- suppressWarnings(vcmoe_fit(
    y ~ z1 | x1,
    data = sim$data,
    u = "u",
    k = 3,
    bandwidth = 0.50,
    u_grid = c(0.30, 0.70),
    engine = "joint_path_em",
    control = list(maxit = 6L, n_starts = 1L, seed = 8922, warn_ambiguous = FALSE)
  ))
  band <- suppressWarnings(vcmoe_confband(
    fit,
    coefficient_set = c("expert", "gating", "sigma"),
    strict = FALSE
  ))
  expect_s3_class(band, "vcmoe_confband")
  expect_identical(band$settings$engine_id, "joint_path_em")

  reduced <- suppressWarnings(vcmoe_fit_reduced(
    fit,
    constrain = "gating_constant",
    control = list(maxit = 3L, strict = FALSE)
  ))
  expect_identical(
    reduced$diagnostics$glrt_null_engine,
    "joint_path_em_constrained_null"
  )
  expect_true(all(is.finite(predict(reduced, type = "mean"))))

  glrt <- suppressWarnings(vcmoe_glrt(
    fit,
    sim$data,
    test = "coefficient",
    coefficient_set = "expert",
    component = 1,
    term = "z1",
    calibration = "none",
    control = list(maxit = 3L, strict = FALSE)
  ))
  expect_identical(
    glrt$null_fit$diagnostics$glrt_null_engine,
    "joint_path_em_constrained_null"
  )

  boot <- suppressWarnings(vcmoe_bootstrap(
    fit,
    sim$data,
    B = 2L,
    seed = 8923,
    min_successful = 2L,
    keep_fits = TRUE,
    control = list(maxit = 3L, n_starts = 1L, warn_ambiguous = FALSE)
  ))
  expect_length(boot$fits, 2L)
  expect_true(all(vapply(
    boot$fits,
    function(x) identical(x$engine_id, "joint_path_em"),
    logical(1L)
  )))
})
