# Package index

## Model fitting

Fit the transition and baseline models by maximum likelihood.

- [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md)
  : Fit the regression models for separable treatment effects
- [`coef(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`vcov(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`nobs(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`logLik(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`confint(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`print(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`summary(`*`<biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  [`print(`*`<summary.biscal>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
  : Methods for fitted biscal models

## Separable treatment effects

Estimates, standard errors and tests, given covariate values or averaged
over them.

- [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) :
  Separable treatment effects and tests given covariate values
- [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md) :
  Population-level separable treatment effects and tests
- [`print(`*`<biscal_aste>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-effects.md)
  [`print(`*`<biscal_easte>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-effects.md)
  [`as.data.frame(`*`<biscal_aste>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-effects.md)
  [`as.data.frame(`*`<biscal_easte>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/biscal-effects.md)
  : Print and convert separable-effect estimates
- [`plot(`*`<biscal_aste>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/plot.biscal_aste.md)
  [`plot(`*`<biscal_easte>`*`)`](https://pingbohu43.github.io/BiSCaL/reference/plot.biscal_aste.md)
  : Plot separable treatment effects

## Counterfactual quantities

- [`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md) :
  Generalized odds ratios of the two responses
- [`cf_prob()`](https://pingbohu43.github.io/BiSCaL/reference/cf_prob.md)
  : Counterfactual joint distribution of the two responses

## Data and simulation

- [`simtrial`](https://pingbohu43.github.io/BiSCaL/reference/simtrial.md)
  : Simulated longitudinal trial with two binary responses
- [`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md)
  : Simulate longitudinal bivariate responses
- [`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md)
  : Coefficients of the simulation designs of the paper
- [`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md)
  : True effects of a simulation design
