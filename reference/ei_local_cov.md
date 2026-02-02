# Estimate the residual covariance of the local estimands

Under a slightly stronger coarsening at random assumption (applying to
second moments), *and* an assumption of homoskedasticity in the
covariates, this function estimates the covariance matrix of the local
estimands \\\beta\_{gj}=\mathbb{E}\[Y \| X_j=1, Z=z_g, G=g\]\\ around
their local mean. See the reference for more detail.

## Usage

``` r
ei_local_cov(regr, data)
```

## Arguments

- regr:

  A fitted regression model, from
  [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md),
  or another kind of regression model wrapped with
  [`ei_wrap_model()`](https://corymccartan.com/seine/reference/ei_wrap_model.md).

- data:

  The data frame, matrix, or
  [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md) object
  that was used to fit the regression.

## Value

A covariance matrix. The variables are ordered by predictor within
outcome, e.g. (Y1\|X1, Y1\|X2, ..., Y2\|X1, Y2\|X2, ...).

## Details

Homoskedasticity in the covariates implies that the variance of the
residuals depends linearly on the entries of \\\bar{X}\bar{X}^\top\\.
This function fits an auto-tuned ridge regression of the empirical
second moments of the residuals on these predictors, and uses the
polarization identity discussed in the references to estimate the
covariance for each local estimand. When the estiamated covariance is
not positive semidefinite, it is projected onto the cone of positive
semidefinite matrices.

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).
