<div id="main" class="col-md-9" role="main">

# Analytic-style confidence bands for a VCMoE fit

<div class="ref-description section level2">

Analytic-style confidence bands for a VCMoE fit

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_confband(
  fit,
  data = NULL,
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  strict = TRUE,
  control = list()
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit with `k = 2:10`.

-   data:

    Optional original data frame. The current implementation uses the
    data stored in `fit$fitted`; refit with `keep_data = TRUE` if
    needed.

-   level:

    Confidence level.

-   type:

    Interval columns to expose as `lower` and `upper`.

-   coefficient\_set:

    Coefficient blocks to return.

-   strict:

    Whether weak local fits should return blocked intervals.

-   control:

    Optional development inference controls. HC0 is the only active
    covariance adjustment.

</div>

<div class="section level2">

## Value

A `vcmoe_confband` object with interval and diagnostic data frames.

</div>

<div class="section level2">

## Details

For `engine = "joint_path_em"`, the covariance follows the JASA observed
local-likelihood asymptotic sandwich plug-in. It does not include
shared-path, label-selection, or finite-grid cross-grid responsibility
uncertainty. Joint-path convergence and the returned score-imbalance
diagnostics should therefore be inspected. The returned metadata
identifies the covariance target, estimator/covariance match, omitted
uncertainty, and coverage-theory scope; no bias or boundary correction
is applied.

</div>

</div>
