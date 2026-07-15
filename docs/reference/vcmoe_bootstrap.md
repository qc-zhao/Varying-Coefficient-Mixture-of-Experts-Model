<div id="main" class="col-md-9" role="main">

# Parametric bootstrap inference for a VCMoE fit

<div class="ref-description section level2">

Parametric bootstrap inference for a VCMoE fit

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_bootstrap(
  fit,
  data,
  u = NULL,
  B = 200L,
  coefficient_set = c("expert", "gating"),
  seed = NULL,
  control = list(),
  min_successful = max(20L, ceiling(0.5 * B)),
  keep_fits = FALSE,
  verbose = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit with `k = 2:10`.

-   data:

    Original data frame used to fit `fit`.

-   u:

    Optional original `u` values or column name.

-   B:

    Number of bootstrap replicates.

-   coefficient\_set:

    Coefficient sets to store.

-   seed:

    Optional random seed.

-   control:

    Control overrides for bootstrap refits.

-   min\_successful:

    Minimum successful replicates for reliable inference.

-   keep\_fits:

    Whether to store successful bootstrap fit objects.

-   verbose:

    Whether to message progress.

</div>

<div class="section level2">

## Value

An object of class `vcmoe_bootstrap`.

</div>

<div class="section level2">

## Details

Bootstrap refits preserve the reference fitting engine. A joint-path
reference is therefore refitted with
`vcmoe_fit(..., engine = "joint_path_em")` rather than silently falling
back to local-grid EM.

</div>

</div>
