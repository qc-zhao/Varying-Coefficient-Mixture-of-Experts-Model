<div id="main" class="col-md-9" role="main">

# Select a VCMoE bandwidth by K-fold cross-validation

<div class="ref-description section level2">

Select a VCMoE bandwidth by K-fold cross-validation

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_select_bandwidth(
  formula,
  data,
  u,
  k = 2L,
  family = "gaussian",
  bandwidth_grid = NULL,
  folds = 5L,
  u_grid = NULL,
  control = list(),
  label = "align",
  parameterization = "a1_epanechnikov_scaled",
  u_scale = c("unit", "none"),
  seed = NULL,
  refit = TRUE,
  engine = c("local_grid_em", "joint_path_em")
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    A formula of the form `y ~ expert_terms | gating_terms`.

-   data:

    A data frame.

-   u:

    Continuous index column name or numeric vector.

-   k:

    Number of mixture components. Values from 2 through 10 are
    supported.

-   family:

    Model family. `"gaussian"`, `"binomial"`, and `"negative-binomial"`
    are supported.

-   bandwidth\_grid:

    Candidate bandwidth values. If `NULL`, uses multiples of the default
    bandwidth.

-   folds:

    Number of random cross-validation folds.

-   u\_grid:

    Grid where coefficient functions are estimated.

-   control:

    Named list passed to the selected fitting engine.

-   label:

    Label strategy passed to the selected fitting engine.

-   parameterization:

    Estimator convention passed to the selected fitting engine.

-   u\_scale:

    `u` scaling strategy passed to the selected fitting engine.

-   seed:

    Optional random seed for fold assignment and, when `control$seed` is
    absent, deterministic CV refits.

-   refit:

    Whether to refit the final model on all data using the selected
    bandwidth.

-   engine:

    Fitting engine used for every cross-validation fold and the optional
    final refit. The default `"local_grid_em"` preserves the 0.1.0
    behavior; `"joint_path_em"` uses the joint-path engine of
    `vcmoe_fit()`.

</div>

<div class="section level2">

## Value

An object of class `vcmoe_bandwidth_selection`.

</div>

</div>
