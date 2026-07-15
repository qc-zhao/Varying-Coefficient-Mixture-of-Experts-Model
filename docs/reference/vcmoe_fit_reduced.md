<div id="main" class="col-md-9" role="main">

# Fit a block-constant reduced VCMoE model

<div class="ref-description section level2">

Refits a VCMoE object under a block-constant coefficient constraint. A
local-grid reference uses the established constrained BFGS optimizer; a
joint-path reference uses a paper-inspired sample-weighted grid
projection and applies the selected constraint after every M-step. This
projection is not a generic constrained optimizer, and its diagnostic
likelihood need not be monotone.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
vcmoe_fit_reduced(
  fit,
  constrain = c("gating_constant", "expert_constant", "all_constant"),
  control = list()
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   fit:

    A `vcmoe` fit.

-   constrain:

    Constraint to impose. `"gating_constant"` freezes the gating
    contrast functions in `u`; `"expert_constant"` freezes expert mean
    functions and family dispersion paths; `"all_constant"` freezes all
    fitted coefficient functions.

-   control:

    Controls for constrained null fitting. With the default
    `strict = TRUE`, a nonconverged reduced fit is rejected.

</div>

<div class="section level2">

## Value

A reduced object of class `vcmoe` with recomputed posterior and
likelihood caches.

</div>

</div>
