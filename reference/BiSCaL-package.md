# BiSCaL: bivariate separable causal effects in longitudinal studies

Causal inference for two ordinal (including binary) responses measured
repeatedly over time, for example two cholesterol measures in a clinical
trial, when the responses may be missing intermittently and follow-up
may be censored. `BiSCaL` implements the methods of Hu and Yi (2026).
The treatment \\A\\ is conceptualized as two components, \\A_1\\ acting
on the first response and \\A_2\\ on the second (\\A = a\\ corresponds
to \\A_1 = A_2 = a\\). Effects are defined at the end time point \\K\\
through generalized odds ratios (GOR) of the two responses: the total
treatment effect \\TTE_K\\ compares \\(A_1, A_2) = (1, 1)\\ with \\(0,
0)\\, and the separable treatment effect \\STE\_{j,K}\\ changes only the
component acting on response \\j\\; they satisfy \\STE\_{1,K} +
STE\_{2,K} = 2\\TTE_K\\.

## Workflow

1.  Arrange the data in long format, one row per subject and visit, with
    the visit index `0, 1, ..., K`, the treatment and baseline
    covariates; code missing responses as `NA` (see
    [simtrial](https://pingbohu43.github.io/BiSCaL/reference/simtrial.md)
    for an example).

2.  Fit the transition and baseline models with
    [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md);
    inspect them with [`summary()`](https://rdrr.io/r/base/summary.html)
    ([biscal-methods](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)).

3.  Estimate and test the effects for given covariate values with
    [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md),
    or averaged over a covariate distribution with
    [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md);
    [`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md) and
    [`cf_prob()`](https://pingbohu43.github.io/BiSCaL/reference/cf_prob.md)
    give the underlying generalized odds ratios and counterfactual
    distributions. Plot effects over follow-up with
    [`plot.biscal_aste()`](https://pingbohu43.github.io/BiSCaL/reference/plot.biscal_aste.md).

[`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md),
[`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md)
and
[`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md)
reproduce the simulation designs of the paper. See
[`vignette("BiSCaL")`](https://pingbohu43.github.io/BiSCaL/articles/BiSCaL.md)
for a worked analysis.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

Useful links:

- <https://github.com/pingbohu43/BiSCaL>

- <https://pingbohu43.github.io/BiSCaL/>

- Report bugs at <https://github.com/pingbohu43/BiSCaL/issues>

## Author

**Maintainer**: Pingbo Hu <pingbo.hu@gmail.com> \[copyright holder\]
