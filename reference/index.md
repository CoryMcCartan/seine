# Package index

## Ecological inference

Estimate global and local quantities and perform sensitivity analyses

- [`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md)
  [`as.matrix(`*`<ei_est>`*`)`](https://corymccartan.com/seine/reference/ei_est.md)
  [`vcov(`*`<ei_est>`*`)`](https://corymccartan.com/seine/reference/ei_est.md)
  [`nobs(`*`<ei_est>`*`)`](https://corymccartan.com/seine/reference/ei_est.md)
  : Estimate ecological quantities
- [`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md)
  [`as.array(`*`<ei_est_local>`*`)`](https://corymccartan.com/seine/reference/ei_est_local.md)
  : Produce local ecological estimates
- [`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md) :
  Conduct a sensitivity analysis for estimated ecological quantities
- [`ei_sens_rv()`](https://corymccartan.com/seine/reference/ei_sens_rv.md)
  : Robustness values for ecological inference
- [`ei_bench()`](https://corymccartan.com/seine/reference/ei_bench.md) :
  Benchmark sensitivity parameters from observed covariates
- [`plot(`*`<ei_sens>`*`)`](https://corymccartan.com/seine/reference/plot.ei_sens.md)
  : Bias contour plot for ecological inference estimates

## Ecological modeling

Fit regression and weighting models to aggregate data to use in
inference

- [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) :
  Fit an ecological inference regression model

- [`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md) :
  Estimate Riesz representer for ecological inference

- [`ei_resid_cov()`](https://corymccartan.com/seine/reference/ei_resid_cov.md)
  : Estimate the residual covariance of the local estimands

- [`ei_wgt_unif()`](https://corymccartan.com/seine/reference/ei_wgt.md)
  [`ei_wgt_prop()`](https://corymccartan.com/seine/reference/ei_wgt.md)
  [`ei_wgt_sqrt()`](https://corymccartan.com/seine/reference/ei_wgt.md)
  : Estimation weighting models

- [`predict(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  [`fitted(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  [`residuals(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  [`vcov(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  [`summary(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  [`weights(`*`<ei_ridge>`*`)`](https://corymccartan.com/seine/reference/ridge-methods.md)
  :

  Methods for `ei_ridge` models

- [`weights(`*`<ei_riesz>`*`)`](https://corymccartan.com/seine/reference/weights.ei_riesz.md)
  : Extract Riesz representer weights

- [`ei_ridge_impl()`](https://corymccartan.com/seine/reference/ei-impl.md)
  [`ei_riesz_impl()`](https://corymccartan.com/seine/reference/ei-impl.md)
  :

  Low-level implementations of
  [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md)
  and
  [`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md)

- [`ei_wrap_model()`](https://corymccartan.com/seine/reference/ei_wrap_model.md)
  :

  Wrap another predictive model for use in `ei_est`

## Preprocessing

Functions to set up EI problems and tidy data

- [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  [`weights(`*`<ei_spec>`*`)`](https://corymccartan.com/seine/reference/ei_spec.md)
  : Specify an ecological inference problem
- [`plot(`*`<ei_spec>`*`)`](https://corymccartan.com/seine/reference/plot.ei_spec.md)
  : Plot an EI specification
- [`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
  : Convert counts to proportions

## Data

- [`elec_1968`](https://corymccartan.com/seine/reference/elec_1968.md) :
  1968 U.S. presidential election data
- [`ei_synthetic()`](https://corymccartan.com/seine/reference/ei_synthetic.md)
  : Generate synthetic ecological data
