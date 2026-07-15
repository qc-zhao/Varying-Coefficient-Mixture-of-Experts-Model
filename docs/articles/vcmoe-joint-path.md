<div id="main" class="col-md-9" role="main">

# Joint-Path EM

`VCMoE` provides two fitting engines. The default `local_grid_em` engine
fits each grid point independently and aligns component labels
afterward. The `joint_path_em` engine instead maintains one
observation-level responsibility matrix across its outer EM loop. At
each iteration it updates every local model from that shared path and
then refreshes each observation from its nearest grid point.

The engines share the same fitted-object and prediction interfaces.
Joint-path EM is opt-in and the default behavior of `vcmoe_fit()` is
unchanged.

<div class="section level2">

## Fit a joint-path model

<div id="cb1" class="sourceCode">

``` r
library(VCMoE)
```

</div>

<div id="cb2" class="sourceCode">

``` r
sim <- simulate_vcmoe_gaussian(
  n = 90,
  k = 2,
  seed = 21,
  separation = 1.7,
  scenario = "well_separated"
)

fit <- vcmoe_fit(
  y ~ z1 | x1,
  data = sim$data,
  u = "u",
  k = 2,
  family = "gaussian",
  bandwidth = 0.40,
  u_grid = c(0.25, 0.50, 0.75),
  engine = "joint_path_em",
  control = list(
    maxit = 12,
    n_starts = 1,
    seed = 22,
    warn_ambiguous = FALSE
  )
)

fit$engine_id
#> [1] "joint_path_em"
head(predict(fit, type = "posterior"))
#>              [,1]         [,2]
#> [1,] 1.849120e-02 9.815088e-01
#> [2,] 9.931178e-01 6.882247e-03
#> [3,] 2.740944e-01 7.259056e-01
#> [4,] 9.979856e-01 2.014450e-03
#> [5,] 3.426019e-09 1.000000e+00
#> [6,] 1.000000e+00 1.809757e-16
```

</div>

Joint-path diagnostics include the selected start, outer iteration
trace, and the number of observations assigned to each grid point.

<div id="cb3" class="sourceCode">

``` r
fit$diagnostics$joint_path_converged
#> [1] FALSE
fit$diagnostics$joint_path_assignment
#>   grid_id    u u_original n_assigned
#> 1       1 0.25  0.2562670         26
#> 2       2 0.50  0.5011344         20
#> 3       3 0.75  0.7460017         44
tail(fit$diagnostics$joint_path_trace)
#>    iteration objective objective_delta posterior_delta parameter_delta
#> 7          7 -98.90620       0.6030760      0.12659834      0.12075461
#> 8          8 -98.49352       0.4126819      0.15574634      0.08947613
#> 9          9 -98.22225       0.2712709      0.14582814      0.06661256
#> 10        10 -97.99489       0.2273618      0.11564046      0.06812650
#> 11        11 -97.78657       0.2083137      0.07866552      0.07119413
#> 12        12 -97.61032       0.1762492      0.04935146      0.06572002
#>    converged
#> 7      FALSE
#> 8      FALSE
#> 9      FALSE
#> 10     FALSE
#> 11     FALSE
#> 12     FALSE
```

</div>

The trace’s `objective` column is the sample-level nearest-grid
log-likelihood. It is a diagnostic, not a monotonicity guarantee: the
label-consistent M-step and nearest-grid responsibility update can
decrease it. Convergence is therefore based on posterior and parameter
deltas.

</div>

<div class="section level2">

## Guardrails

Runtime grows approximately with the number of observations, grid
points, EM iterations, and starts. By default, joint-path EM rejects
more than 100 grid points or a grid-to-observation ratio above 0.2. Set
`control$allow_dense_u_grid = TRUE` only after explicitly accepting that
cost. A grid close to `unique(u)` also produces a warning because it
approaches one local model per observation.

If no start converges by `control$maxit`, the fit is returned with an
explicit warning so that its finite interim estimates can still be
inspected. Treat such estimates as provisional and review
`joint_path_trace` before inference.

</div>

<div class="section level2">

## Inference and bandwidth selection

Analytic confidence bands use an observed local-likelihood sandwich
plug-in evaluated at the fitted coefficient arrays. For joint-path fits
this is a local-curvature approximation, not the covariance of the
complete finite-grid shared-path estimator. The returned metadata states
that shared-path, label-selection, and finite-grid cross-grid
responsibility uncertainty are not included and reports score imbalance.
Parametric bootstrap refits, GLRT full and null fits, reduced fits, and
bandwidth-selection refits preserve the selected engine.

For joint-path GLRT, the constant-coefficient null uses a paper-inspired
update: constrained local paths are replaced after every M-step by their
mean weighted by nearest-grid observation-assignment counts. This is a
projected estimator, not a generic constrained MLE. The statistic
evaluates every observation once at its nearest grid point, which is a
grid approximation to the sample-point criterion. Failed or nonconverged
null fits are blocked. GLRT calibration defaults to `"none"`; request
analytic calibration only as a documented approximation, or use
engine-preserving parametric bootstrap calibration.

<div id="cb4" class="sourceCode">

``` r
band <- vcmoe_confband(fit, strict = FALSE)

boot <- vcmoe_bootstrap(
  fit,
  data = sim$data,
  B = 200,
  seed = 23
)

test <- vcmoe_glrt(
  fit,
  data = sim$data,
  test = "coefficient",
  coefficient_set = "expert",
  component = 1,
  term = "z1",
  calibration = "none"
)

selection <- vcmoe_select_bandwidth(
  y ~ z1 | x1,
  data = sim$data,
  u = "u",
  family = "gaussian",
  bandwidth_grid = c(0.30, 0.40, 0.50),
  folds = 3,
  u_grid = c(0.25, 0.50, 0.75),
  engine = "joint_path_em"
)
```

</div>

Always inspect convergence, component sizes, label ambiguity, and
null-fit diagnostics before interpreting inference from either engine.

</div>

</div>
