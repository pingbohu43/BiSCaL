# BiSCaL 0.1.0

* First release, implementing the methods of Hu and Yi (2026,
  *Statistics in Medicine*, <doi:10.1002/sim.70320>):
  * `biscal()` fits the transition and baseline models by maximum likelihood,
    for binary or ordinal responses, with intermittent missingness and
    right-censoring;
  * `aste()` and `easte()` estimate the total and separable treatment effects,
    and test the difference of the absolute separable effects, given covariate
    values or averaged over a covariate distribution;
  * `gor()` and `cf_prob()` give the generalized odds ratios and the
    counterfactual joint distributions;
  * `simulate_biscal()`, `biscal_coefs()` and `biscal_truth()` reproduce the
    simulation designs of the paper;
  * the simulated example data set `simtrial`.
* Checked against the original research code of the paper (regression tests
  in `tests/testthat/test-regression.R`).
