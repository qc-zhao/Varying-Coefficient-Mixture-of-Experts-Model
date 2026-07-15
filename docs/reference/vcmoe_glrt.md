<div id="main" class="col-md-9" role="main">

# Generalized likelihood-ratio test for VCMoE coefficient variation

<div class="ref-description section level2">

Generalized likelihood-ratio test for VCMoE coefficient variation

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_glrt(
  fit,
  data,
  test = c("coefficient", "constant_block", "constant_all"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  component = NULL,
  term = NULL,
  calibration = c("none", "bootstrap", "analytic_epanechnikov", "both",
    "parametric_bootstrap"),
  B = 200L,
  seed = NULL,
  control = list(),
  refit_control = list(),
  verbose = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit.

-   data:

    Original data frame used to fit `fit`.

-   test:

    Test type. `"coefficient"` tests one coefficient function;
    `"constant_block"` tests all expert or gating functions jointly;
    `"constant_all"` tests all fitted coefficient functions jointly.

-   coefficient\_set:

    Coefficient block for coefficient-specific or block-constant tests.

-   component:

    Component label or index for coefficient-specific tests.

-   term:

    Term name for coefficient-specific tests.

-   calibration:

    Calibration method. The default `"none"` returns the statistic
    without attaching a reference distribution.
    `"analytic_epanechnikov"` uses the Epanechnikov modified chi-square
    calibration; `"bootstrap"` uses parametric bootstrap calibration;
    `"both"` reports both. The analytic calibration is retained as an
    explicitly requested approximation because the implemented statistic
    is not identical to the manuscript criterion.

-   B:

    Number of bootstrap calibration replicates.

-   seed:

    Optional random seed.

-   control:

    Controls for constrained null optimization and diagnostics.

-   refit\_control:

    Controls overriding bootstrap full-model refits.

-   verbose:

    Whether to message bootstrap progress.

</div>

<div class="section level2">

## Value

A `vcmoe_glrt` object.

</div>

<div class="section level2">

## Details

Local-grid fits retain the 0.1.0 constrained BFGS null optimizer.
Joint-path fits use a paper-inspired sample-weighted grid-projected
null: after every M-step, each constrained coefficient path is replaced
by its mean weighted by the number of observations assigned to each
nearest grid point, and constrained local slopes are set to zero. Its
statistic compares sample-level likelihood contributions evaluated at
each observation's nearest grid point. The projected update is not a
generic constrained optimizer and its diagnostic likelihood trace need
not be monotone. Bootstrap calibration preserves both the full-fit
engine and its matching null engine.

</div>

</div>
