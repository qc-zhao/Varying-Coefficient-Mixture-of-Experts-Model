test_that("0.2.0 exposes stable engines without research-only arguments", {
  fit_formals <- names(formals(vcmoe_fit))
  expect_true("engine" %in% fit_formals)
  expect_false("confirm_experimental" %in% fit_formals)
  expect_false("composite_alpha" %in% fit_formals)
  expect_false("response_weights" %in% fit_formals)
  expect_false("ordered_u_starts" %in% fit_formals)

  exports <- getNamespaceExports("VCMoE")
  expect_true("vcmoe_fit_reduced" %in% exports)
  expect_false("vcmoe_fit_experimental" %in% exports)

  sim <- simulate_vcmoe_gaussian(n = 30, k = 2, seed = 7101)
  expect_error(
    vcmoe_fit(
      formula = list(y = y ~ z1),
      data = sim$data,
      u = "u"
    ),
    "formula"
  )
})

test_that("package implementation has no removed research API residue", {
  package_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  scan_paths <- c(
    file.path(package_root, "R"),
    file.path(package_root, "man"),
    file.path(package_root, "vignettes"),
    file.path(package_root, "NAMESPACE"),
    file.path(package_root, "README.md")
  )
  files <- unlist(lapply(scan_paths, function(path) {
    if (dir.exists(path)) {
      list.files(path, recursive = TRUE, full.names = TRUE)
    } else if (file.exists(path)) {
      path
    } else {
      character(0L)
    }
  }), use.names = FALSE)
  text <- unlist(lapply(files, readLines, warn = FALSE), use.names = FALSE)
  removed <- c(
    "composite_alpha",
    "response_weights",
    "confirm_experimental",
    "experimental_joint_path_em",
    "allow_experimental_joint_path_inference",
    "ordered_u_starts",
    "legacy_gaussian_raw"
  )
  for (token in removed) {
    expect_false(any(grepl(token, text, fixed = TRUE)), info = token)
  }
})
