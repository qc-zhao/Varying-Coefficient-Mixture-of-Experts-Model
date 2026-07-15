<div id="main" class="col-md-9" role="main">

# Changelog

<div class="section level2">

## VCMoE 0.2.0

<div class="section level3">

### Fitting engines

-   Adds stable single-response joint-path EM through
    `vcmoe_fit(..., engine = "joint_path_em")`.
-   Keeps `engine = "local_grid_em"` as the default and preserves the
    0.1.0 local-grid numerical path when `engine` is omitted.
-   Adds dense-grid cost guardrails, nearest-grid observation
    assignment, joint-path convergence traces, and optional progress
    logging.

</div>

<div class="section level3">

### Inference and selection

-   Makes analytic confidence bands available through the common
    fitted-object interface for both engines. Joint-path bands identify
    their JASA local-likelihood sandwich plug-in scope and expose
    score-imbalance diagnostics.
-   Makes parametric bootstrap refits and bandwidth selection
    engine-aware.
-   Adds sample-weighted projected joint-path constrained-null fitting
    for GLRT and `vcmoe_fit_reduced()` while retaining the established
    local-grid BFGS null optimizer. Joint-path GLRT uses sample-level
    nearest-grid likelihood contributions and blocks failed or
    nonconverged null fits.
-   Changes the GLRT calibration default to `"none"`. The statistic is
    available without a reference distribution; analytic Epanechnikov
    calibration is an explicit approximation, and parametric bootstrap
    calibration is engine-preserving.
-   Adds criterion, projection, constrained-MLE, calibration, and
    theory-status metadata to GLRT results. Analytic confidence bands
    now state whether their covariance target matches the fitted
    estimator and which shared-path and label uncertainties are
    excluded.

</div>

<div class="section level3">

### Scope

-   Version 0.2.0 remains a single-response package for Gaussian,
    Binomial, and Negative-Binomial models with 2 through 10 components.

</div>

</div>

</div>
