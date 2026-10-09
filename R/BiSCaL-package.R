#' BiSCaL: bivariate separable causal effects in longitudinal studies
#'
#' Causal inference for two ordinal (including binary) responses measured
#' repeatedly over time, for example two cholesterol measures in a clinical
#' trial, when the responses may be missing intermittently and follow-up may be
#' censored. `BiSCaL` implements the methods of Hu and Yi (2026). The treatment
#' \eqn{A} is conceptualized as two components, \eqn{A_1} acting on the first
#' response and \eqn{A_2} on the second (\eqn{A = a} corresponds to
#' \eqn{A_1 = A_2 = a}). Effects are defined at the end time point \eqn{K}
#' through generalized odds ratios (GOR) of the two responses: the total
#' treatment effect \eqn{TTE_K} compares \eqn{(A_1, A_2) = (1, 1)} with
#' \eqn{(0, 0)}, and the separable treatment effect \eqn{STE_{j,K}} changes
#' only the component acting on response \eqn{j}; they satisfy
#' \eqn{STE_{1,K} + STE_{2,K} = 2\,TTE_K}.
#'
#' @section Workflow:
#' 1. Arrange the data in long format, one row per subject and visit, with the
#'    visit index `0, 1, ..., K`, the treatment and baseline covariates; code
#'    missing responses as `NA` (see [simtrial] for an example).
#' 2. Fit the transition and baseline models with [biscal()]; inspect them
#'    with `summary()` ([biscal-methods]).
#' 3. Estimate and test the effects for given covariate values with [aste()],
#'    or averaged over a covariate distribution with [easte()]; [gor()] and
#'    [cf_prob()] give the underlying generalized odds ratios and
#'    counterfactual distributions. Plot effects over follow-up with
#'    [plot.biscal_aste()].
#'
#' [simulate_biscal()], [biscal_coefs()] and [biscal_truth()] reproduce the
#' simulation designs of the paper. See `vignette("BiSCaL")` for a worked
#' analysis.
#'
#' @template ref
#' @keywords internal
#' @importFrom stats coef confint logLik nobs vcov
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
