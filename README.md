# VCMoE

<!-- badges: start -->
[![R](https://img.shields.io/badge/R-package-276DC3.svg)](https://www.r-project.org/)
[![CRAN status](https://www.r-pkg.org/badges/version/VCMoE)](https://cran.r-project.org/package=VCMoE)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://qc-zhao.github.io/VCMoE/)
[![arXiv](https://img.shields.io/badge/arXiv-2601.01699-B31B1B.svg)](https://arxiv.org/abs/2601.01699)
[![Issues welcome](https://img.shields.io/badge/issues-welcome-brightgreen.svg)](https://github.com/qc-zhao/VCMoE/issues)
<!-- badges: end -->

## Varying-Coefficient Mixture-of-Experts Models

`VCMoE` is an R package for fitting varying-coefficient
mixture-of-experts models. It supports Gaussian, Binomial, and
Negative-Binomial responses, with local-linear estimation, component label
alignment, bandwidth selection, diagnostics, confidence bands, bootstrap
inference, and generalized likelihood-ratio tests.

Version 0.2.0 provides two fitting engines. The default
`engine = "local_grid_em"` preserves the original behavior, while
`engine = "joint_path_em"` updates one observation-level responsibility path
across all grid points. Joint-path fits use the same prediction, diagnostics,
confidence-band, bootstrap, GLRT, and bandwidth-selection interfaces.

The package is intended for problems where component-specific response
relationships and component probabilities change along a continuous coordinate,
such as time, pseudotime, dose, or spatial location.

## Installation

Install the released version from CRAN:

```r
install.packages("VCMoE")
```

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("qc-zhao/VCMoE")
```

Load the package:

```r
library(VCMoE)
```

Need help with installation or usage? Please open a GitHub issue:

<https://github.com/qc-zhao/VCMoE/issues>

## Quick Start

```r
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

Select joint-path EM explicitly when one responsibility path should be shared
across the local coefficient models:

```r
joint_fit <- vcmoe_fit(
  y ~ z1 | x1,
  data = sim$data,
  u = sim$data$u,
  k = 2,
  family = "gaussian",
  bandwidth = 0.25,
  engine = "joint_path_em"
)
```

Joint-path runtime grows with the number of observations, grid points, and EM
iterations. Dense grids are rejected by default; override the guard only after
estimating the computational cost. Its nearest-grid sample log-likelihood trace
is diagnostic and need not increase at every iteration; convergence is based on
posterior and parameter deltas, which should always be inspected.

Joint-path analytic-style bands use an observed local-likelihood sandwich
plug-in. They report score-imbalance diagnostics and do not model shared-path,
label-selection, or finite-grid cross-grid responsibility uncertainty.
Joint-path GLRT uses a paper-inspired sample-weighted grid-projected null and
evaluates each observation once at its nearest grid point. This is a documented
grid approximation, not an exact constrained MLE or exact manuscript
criterion. `vcmoe_glrt()` therefore returns an uncalibrated statistic by
default; analytic or bootstrap calibration must be requested explicitly. A
failed or nonconverged null fit is never reported as valid inference.

## Documentation

The full documentation website includes Gaussian, Binomial, and
Negative-Binomial tutorials plus the function reference:

<https://qc-zhao.github.io/VCMoE/>

Useful links:

- [Gaussian simulation tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-gaussian-no-offset.html)
- [Binomial VCMoE Tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-binomial.html)
- [Negative-Binomial VCMoE Tutorial](https://qc-zhao.github.io/VCMoE/articles/vcmoe-negbin.html)
- [Function reference](https://qc-zhao.github.io/VCMoE/reference/index.html)
- [GitHub issues](https://github.com/qc-zhao/VCMoE/issues)

## Citation

Please cite:

Zhao Q, Greenwood CMT, Zhang Q. *Varying-Coefficient Mixture of Experts Model*.
arXiv:2601.01699. <https://arxiv.org/abs/2601.01699>
