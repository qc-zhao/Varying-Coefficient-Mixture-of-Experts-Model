# VCMoE 0.2.0 model record

## Purpose

VCMoE estimates single-response Gaussian, Binomial, and Negative-Binomial
mixture-of-experts models whose expert, gating, and nuisance parameters vary
smoothly with an index `u`.

## Engines

`local_grid_em` remains the default and retains the 0.1.0 implementation: each
grid point is fitted independently, followed by path-level label alignment.

`joint_path_em` maintains one observation-level responsibility matrix across
the outer EM loop. Each outer M-step updates every local model using that
shared responsibility path. Each observation then takes its updated
responsibilities from its nearest grid point. The completed local path is
label-aligned with the same public coefficient and prediction ABI as the
default engine.

## Revision rationale

The joint-path engine implements the legacy/JASA path-level EM interpretation
without changing the default estimator. It records outer iteration traces,
posterior and parameter deltas, grid assignment counts, and selected starts.
Dense grids are guarded because runtime grows roughly with
`iterations * grid points * observations`.

Inference uses the common fitted-object ABI. Bootstrap and bandwidth refits
preserve the source engine. Joint-path GLRT null fits apply the selected
constraint after every M-step by replacing constrained paths with their
nearest-grid assignment-frequency weighted mean; local-grid null fits keep the
0.1.0 BFGS path. This is a paper-inspired projected estimator, not an exact
constrained MLE. GLRT calibration defaults to `none`; analytic calibration is
an explicitly requested approximation, while parametric bootstrap refits keep
the source fitting and null engines.

Analytic confidence bands keep the established local marginal-likelihood
sandwich calculation. For joint-path fits this is a local-curvature plug-in: it
does not include shared responsibility-path, label-selection, or finite-grid
cross-grid coupling uncertainty. Returned metadata states this scope directly.

## Deliberate exclusions

Version 0.2.0 does not include response-vector/list-formula fitting, composite
objectives, response weights, legacy parameterization modes, ordered-start
research options, marker assignment, or private real-data runners/results.
