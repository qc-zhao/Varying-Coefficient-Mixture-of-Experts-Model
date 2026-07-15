<div id="main" class="col-md-9" role="main">

# Fit a varying-coefficient mixture-of-experts model

<div class="ref-description section level2">

Fit a varying-coefficient mixture-of-experts model

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_fit(
  formula,
  data,
  u,
  k = 2L,
  family = "gaussian",
  bandwidth = NULL,
  u_grid = NULL,
  control = list(),
  label = "align",
  parameterization = "a1_epanechnikov_scaled",
  u_scale = c("unit", "none"),
  engine = c("local_grid_em", "joint_path_em"),
  progress = NULL
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    A formula of the form `y ~ expert_terms | gating_terms`. For grouped
    binomial data, use
    `cbind(success, failure) ~ expert_terms | gating_terms`. For
    Negative-Binomial count data, expert-side `offset(log_size_factor)`
    terms are supported.

-   data:

    A data frame.

-   u:

    Continuous index column name or numeric vector.

-   k:

    Number of mixture components. `k = 2` is the primary v0 target;
    `k = 3:10` are high-k candidate support and require diagnostics
    before being treated as stable.

-   family:

    Model family. `"gaussian"`, `"binomial"`, and `"negative-binomial"`
    are implemented.

-   bandwidth:

    Kernel bandwidth. If `NULL`, a Silverman-style default is used for
    `u`.

-   u\_grid:

    Grid where coefficient functions are estimated.

-   control:

    Named list overriding EM and label-alignment settings.

-   label:

    Label strategy. `"align"` uses exact global alignment for `k <= 6`
    and sequential assignment with ambiguity margins for `k >= 7`;
    `"global"` requests exact global alignment when feasible and falls
    back to the same sequential assignment path for `k >= 7`; `"greedy"`
    keeps the older one-step alignment.

-   parameterization:

    Estimator convention. The public package uses
    `"a1_epanechnikov_scaled"`: Epanechnikov density weights
    `K((u-u0)/h)/h` and the scaled local-linear basis `(u-u0)/h`.

-   u\_scale:

    How to transform `u` before fitting. The default `"unit"` maps
    complete-row `u` values to `[0, 1]`; `"none"` leaves `u` on the
    supplied scale. `bandwidth` and `u_grid` are interpreted on the
    transformed analysis scale.

-   engine:

    Fitting engine. The default `"local_grid_em"` preserves the original
    independent local-grid EM path. `"joint_path_em"` updates a shared
    observation-level responsibility path across grid points.

-   progress:

    Joint-path progress reporting. The default `NULL` is silent; use
    `TRUE` for messages or a single CSV file path for structured
    logging. Iteration frequency is controlled by
    `control$progress_every`.

</div>

<div class="section level2">

## Value

An object of class `vcmoe`.

</div>

<div class="section level2">

## Details

Rows with missing or non-finite response, covariates, or `u` are removed
consistently before fitting, with a warning. For single-trial Bernoulli
responses, the default gating ridge is strengthened to
`control$ridge = 1` unless the user explicitly supplies `control$ridge`;
grouped Binomial and other families keep the global default. Joint-path
traces record the sample-level nearest-grid log-likelihood as a
diagnostic criterion. The label-consistent updates do not guarantee that
this diagnostic is monotone; convergence is based on posterior and
parameter deltas instead.

</div>

</div>
