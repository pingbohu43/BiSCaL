# Changelog

## BiSCaL 0.1.0

- First release, implementing the methods of Hu and Yi (2026,
  *Statistics in Medicine*, <doi:10.1002/sim.70320>):
  - [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md)
    fits the transition and baseline models by maximum likelihood, for
    binary or ordinal responses, with intermittent missingness and
    right-censoring;
  - [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md)
    and
    [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md)
    estimate the total and separable treatment effects, and test the
    difference of the absolute separable effects, given covariate values
    or averaged over a covariate distribution;
  - [`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md) and
    [`cf_prob()`](https://pingbohu43.github.io/BiSCaL/reference/cf_prob.md)
    give the generalized odds ratios and the counterfactual joint
    distributions;
  - [`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md),
    [`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md)
    and
    [`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md)
    reproduce the simulation designs of the paper;
  - the simulated example data set `simtrial`.
- Checked against the original research code of the paper (regression
  tests in `tests/testthat/test-regression.R`).
