\name{vcmoe_confband}
\alias{vcmoe_confband}
\title{Analytic-style confidence bands for a VCMoE fit}
\usage{
vcmoe_confband(
  fit,
  data = NULL,
  level = 0.95,
  type = c("pointwise", "simultaneous"),
  coefficient_set = c("expert", "gating", "sigma", "theta"),
  strict = TRUE,
  control = list()
)
}
\arguments{
\item{fit}{A \code{vcmoe} fit with \code{k = 2:10}.}
\item{data}{Optional original data frame. The current implementation uses
the data stored in \code{fit$fitted}; refit with \code{keep_data = TRUE} if
needed.}
\item{level}{Confidence level.}
\item{type}{Interval columns to expose as \code{lower} and \code{upper}.}
\item{coefficient_set}{Coefficient blocks to return.}
\item{strict}{Whether weak local fits should return blocked intervals.}
\item{control}{Optional development inference controls. HC0 is the only
active covariance adjustment.}
}
\value{
A \code{vcmoe_confband} object with interval and diagnostic data frames.
}
\description{
Analytic-style confidence bands for a VCMoE fit
}
\details{
For \code{engine = "joint_path_em"}, the covariance follows the JASA
observed local-likelihood asymptotic sandwich plug-in. It does not include
shared-path, label-selection, or finite-grid cross-grid responsibility
uncertainty. Joint-path convergence and the returned score-imbalance
diagnostics should therefore be inspected. The returned metadata identifies
the covariance target, estimator/covariance match, omitted uncertainty, and
coverage-theory scope; no bias or boundary correction is applied.
}
