# VCMoE 0.2.0 release validation record

## 0.1.0 regression contract

The reference contract is generated from tag `v0.1.0` at commit `559d717`, not
from the release worktree. Frozen training and held-out CSVs cover Gaussian,
grouped Binomial, and Negative-Binomial offset models. Element-level expected
files record expert/gating/nuisance paths, training and held-out marginal,
component, prior, and posterior predictions, local likelihoods, convergence,
selected starts, and label permutations.

All reference fits must converge at every grid point, have no ambiguous label
alignment, and avoid Negative-Binomial theta boundaries. The 0.2.0 default
call, with no engine argument, must reproduce these values within strict
cross-BLAS tolerances. An explicit `engine = "local_grid_em"` fit must equal
the default call.

## Joint-path smoke matrix

The release smoke study covers Gaussian component counts 2 through 10 plus
single-response Binomial and Negative-Binomial offset fits. It checks finite
coefficients and predictions, normalized posterior rows, positive nuisance
parameters, dense-grid guards, progress records, analytic confidence bands,
engine-preserving bootstrap refits, projected-null GLRT/reduced fits, and
engine-aware bandwidth selection.

The inference audit additionally verifies assignment-frequency weighted null
projection, engine-specific likelihood criteria, rejection of failed or
nonconverged null fits, uncalibrated GLRT defaults, and explicit confidence-band
covariance/coverage metadata. Analytic GLRT calibration is treated as an
opt-in approximation rather than a release acceptance target.

## Evaluation outputs

Lightweight replicate and summary CSVs and diagnostic plots are stored under
`release-validation/0.2.0/`, which is excluded from the source tarball. Final
outcomes, package checks, vignette build, pkgdown build, and tarball audit are
recorded after execution.

## Final smoke outcome

The paired release smoke completed 16/16 finite fits across Gaussian K=2, K=3,
and K=10, grouped Binomial K=2, and Negative-Binomial-with-offset K=2. Posterior
normalization and component-size checks passed in every case. Joint-path
prediction RMSE was lower than local-grid RMSE in 7/8 paired datasets; the
median joint/local RMSE ratio was 0.905. Joint-path runtime was higher in the
small K=2 cases, and one of two NB joint-path fits plus both intentionally short
K=10 fits reached their iteration limits. These are recorded warnings rather
than hidden convergence successes.

The inference audit completed 4/4 structural cases (Gaussian and NB, each with
local-grid and joint-path engines). All produced finite predictions and
confidence-band objects; all bootstrap refits retained their source engine; all
four constrained null fits converged. Joint-path nulls reported nearest-grid
sample likelihood, assignment-frequency weighted projection, non-MLE status,
and nearest-grid manuscript approximation metadata. GLRT calibration defaulted
to `none`. The audit release gate was `PASS`.

## Package validation outcome

- The complete `testthat` suite passed, including exact 0.1.0 Gaussian,
  grouped-Binomial, and Negative-Binomial coefficient/prediction fixtures.
- The joint-path GLRT bootstrap test completed 2/2 replicates with
  `joint_path_em` full refits, `joint_path_em_constrained_null` null refits, and
  converged null diagnostics.
- Vignettes and the pkgdown site built successfully;
  `pkgdown::check_pkgdown()` reported no problems.
- The built `VCMoE_0.2.0.tar.gz` passed all 26 source-content checks, including
  exclusion of private paths, runners, result directories, and removed APIs.
- The online `R CMD check --as-cran` completed with 0 errors, 0 warnings, and
  one environment note because the check machine's HTML Tidy is too old for
  validation. Tests, vignette rebuilding, and the PDF manual all passed.
