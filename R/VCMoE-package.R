#' VCMoE: Varying-Coefficient Mixture-of-Experts Models
#'
#' Fit and evaluate single-response Gaussian, Binomial, and Negative-Binomial
#' varying-coefficient mixture-of-experts models.
#'
#' @keywords internal
#' @importFrom stats approx as.formula coef confint dbinom delete.response
#' @importFrom stats dnorm dnbinom kmeans model.frame model.matrix model.response
#' @importFrom stats optim optimize plogis predict quantile rbinom rnbinom
#' @importFrom stats rnorm runif sd setNames terms
#' @importFrom utils globalVariables
"_PACKAGE"
