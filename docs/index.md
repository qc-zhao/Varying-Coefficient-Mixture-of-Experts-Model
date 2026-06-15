<div id="main" class="col-md-9" role="main">

# VCMoE

<div class="section level1">

`VCMoE` fits varying-coefficient mixture-of-experts models in R for
Gaussian, Binomial, and Negative-Binomial responses.

The package is designed for settings where both the expert response
model and the gating probabilities may vary along a continuous
coordinate such as time, pseudotime, dose, or space.

<div class="section level2">

## Model Overview

VCMoE represents each observation as coming from one of `k` latent
components. For component (c), the expert model describes the response
distribution, while the gating model describes the probability of
belonging to that component.

For a continuous coordinate `u`, the gating model has the form:

``` text
Pr(C_i = c | X_i, U_i = u) = softmax_c( X_i gamma_c(u) )
```

The expert mean is modeled by family-specific component functions. For a
Gaussian response, the expert model has the form:

``` text
E(Y_i | C_i = c, Z_i, U_i = u) = Z_i alpha_c(u)
```

The current implementation uses Epanechnikov kernel local fitting with a
scaled local-linear basis. The public interface supports `k = 2:10` with
diagnostic outputs for convergence, label alignment, component overlap,
and inference quality.

</div>

<div class="section level2">

## Installation

<div id="cb3" class="sourceCode">

``` r
install.packages("remotes")
remotes::install_github("qc-zhao/VCMoE")
```

</div>

Then load the package:

<div id="cb4" class="sourceCode">

``` r
library(VCMoE)
```

</div>

</div>

<div class="section level2">

## Quick Gaussian Example

<div id="cb5" class="sourceCode">

``` r
set.seed(1)

sim <- simulate_vcmoe_gaussian(
  n = 300,
  k = 2,
  scenario = "well_separated"
)

fit <- vcmoe_fit(
  y ~ z1 | x1,
  data = sim$data,
  u = sim$data$u,
  k = 2,
  family = "gaussian",
  bandwidth = 0.25
)

coef(fit, "expert")
predict(fit, type = "posterior")
plot_coefficients(fit)
```

</div>

</div>

<div class="section level2">

## Main Workflows

| Task                                                      | Function                                                                            |
|-----------------------------------------------------------|-------------------------------------------------------------------------------------|
| Fit a VCMoE model                                         | `vcmoe_fit()`                                                                       |
| Extract coefficient functions                             | `coef()`                                                                            |
| Predict mean, component means, or posterior probabilities | `predict()`                                                                         |
| Select bandwidth by held-out likelihood                   | `vcmoe_select_bandwidth()`                                                          |
| Bootstrap inference                                       | `vcmoe_bootstrap()`                                                                 |
| Analytic-style confidence bands                           | `vcmoe_confband()`                                                                  |
| Generalized likelihood-ratio tests                        | `vcmoe_glrt()`                                                                      |
| Simulate example data                                     | `simulate_vcmoe_gaussian()`, `simulate_vcmoe_binomial()`, `simulate_vcmoe_negbin()` |

</div>

<div class="section level2">

## Tutorials

Start with the Gaussian tutorial:

[Gaussian inference
tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-gaussian-no-offset.md)

The tutorial shows how to simulate Gaussian data, fit a VCMoE model,
inspect posterior probabilities, plot analytic simultaneous confidence
bands, run bootstrap inference, and choose a bandwidth.

Count-model tutorials are separate:

-   [Binomial model
    tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-binomial.md)
-   [Negative-Binomial model
    tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-negbin.md)

</div>

<div class="section level2">

## Function Reference

The full function reference is available here:

[Reference index](https://qc-zhao.github.io/VCMoE/reference/index.md)

</div>

<div class="section level2">

## Citation

Please cite:

Zhao Q, Greenwood CMT, Zhang Q. *Varying-Coefficient Mixture of Experts
Model*. arXiv:2601.01699. <https://arxiv.org/abs/2601.01699>

</div>

</div>

</div>
