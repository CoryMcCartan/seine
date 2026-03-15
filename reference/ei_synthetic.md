# Generate synthetic ecological data

Samples data from a following truncated Normal ecological model. The
data can be generated completely at random, or can be generated
conditional on provided predictors `x` and/or covariates `z`.

## Usage

``` r
ei_synthetic(
  n,
  p = 0,
  n_x = 2,
  x = n_x:1,
  z = 0.25 * exp(-(seq_len(p) - 1)/2),
  r2_xz = 0.5,
  r2_bz = 0.5,
  b_loc = NULL,
  b_cov = NULL
)
```

## Arguments

- n:

  The number of rows (geographies) to generate. Defaults to the number
  of rows in `x` or `z`, if they are a matrix or data frame.

- p:

  The number of covariates. Defaults to the number of columns in `z`, if
  it is a matrix or data frame, or the length of `z`, if it is a vector
  of singular values.

- n_x:

  The number of predictor variables. Defaults to the number of columns
  in `x`, if it is a matrix or data frame, or the length of `x`, if it
  is a vector of mean parameters for the softmax-transformed Normal
  distribution.

- x:

  Either a matrix or data frame containing the predictor percentages in
  each row, or a vector containing Dirichlet parameters to use in
  sampling predictor percentages.

- z:

  A matrix or data frame containing geography-level covariates, or a
  vector of values to form a Toeplitz covariance matrix for the random
  covariates.

- r2_xz:

  The approximate \\R^2\\ of the covariates `z` and predictors `x`. See
  the model specification for details. If either `r2_xz` or `r2_bz` are
  zero, then there is no confounding, and an unadjusted Goodman
  regression will estimate the global parameters correctly.

- r2_bz:

  The approximate \\R^2\\ of the covariates `z` and unit-level
  parameters `b`. See the model specification for details. If either
  `r2_xz` or `r2_bz` are zero, then there is no confounding, and an
  unadjusted Goodman regression will estimate the global parameters
  correctly.

- b_loc:

  The center of the distribution of geography-level parameters. Defaults
  to a linearly spaced sequence across groups from 0.5 to 0.9. Because
  of the truncation, this will not exactly be the mean of the
  geography-level parameters.

- b_cov:

  The residual covariance matrix for geography-level parameters.
  Defaults to `0.02 * (1 + diag(n_x))`.

## Value

An `ei_spec` object with additional attributes:

- `b_loc` and `b_cov`

- `Lambda` with the coefficients of `z`

- `eta`, the linear predictor for `b`

- `est_true`, the mean of the geography-level parameters, formatted
  similarly to the output from
  [`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md)

- `r2_xz_act` and `r2_bz_act`, containing the actual (sample) \\R^2\\
  values for `x` and `z`, and `b` and `z`, respectively.

## Details

This function samples data from the following truncated Normal
ecological model: \$\$ \begin{pmatrix}x_i\\ z_i\end{pmatrix}
\stackrel{\text{iid}}{\sim} \mathcal{N}\_{\[0,1\]^{n_x} \times
\mathbb{R}^p}\left( \begin{pmatrix}\mu_x\\ 0\end{pmatrix},
\begin{pmatrix}\Sigma_x & \Gamma \\ \Gamma & T\end{pmatrix}\right) \$\$
\$\$ \eta = z_i^\top \Lambda + \mathtt{b\_{loc}} \$\$ \$\$ b_i
\stackrel{\text{iid}}{\sim} \mathcal{N}\_{\[0, 1\]^{n_x}}(\eta,
\mathtt{B\_{cov}}) \$\$ \$\$ y_i = b_i^\top x_i, \$\$ where \\\mu_x\\
and \\\Sigma_x\\ are the mean and covariance of the Normal approximation
to a Dirichlet distribution with parameters supplied by the `x` argument
below, and \\\Gamma\\, \\T\\, and \\\Gamma\\ are matrices sampled to
have certain properties, as described below. The subscripts on
\\\mathcal{N}\\ indicate truncation; i.e., both the predictors `x` and
the unit-level parameters `b` are truncated to the *n_x*-dimensional
hypercube.

The matrix \\T\\ is a symmetric Toeplitz matrix with diagonals provided
by the `z` argument. Generally, a decreasing set of nonnegative values
will be sufficient for a positive definite \\T\\.

The matrices \\\Gamma\\ and \\\Lambda\\ are initially filled with
independent samples from a standard Normal distribution. \\\Gamma\\ is
then projected so that its rows sum to zero, preserving the sum-to-1
requirement on `x`, and so that its columns are scaled to produce the
correct \\R^2\\ value matching `r2_xz`. The matrix \\\Lambda\\ is
likewise scaled to produce the correct \\R^2\\ value matching `r2_bz`.
Due to the truncation in the sampling of `x` and `b`, the in-sample
\\R^2\\ values will generally be slightly smaller than the provided
arguments.

Aspects of the model can be replaced with data provided to the function.
If `x` or `z` is provided as a matrix or data frame, then the other
value is sampled from its marginal distribution. If both are provided,
then the first row of the model is skipped.

## Examples

``` r
ei_synthetic(n = 10)
#> EI Specification
#> • Predictors: `x1` and `x2`
#> • Outcome: `y`
#> • Covariates: none
#> # A tibble: 10 × 3
#>       y    x1    x2
#>   <dbl> <dbl> <dbl>
#> 1 0.305 0.650 0.350
#> 2 0.735 0.737 0.263
#> 3 0.713 0.562 0.438
#> 4 0.544 0.470 0.530
#> 5 0.737 0.824 0.176
#> # ℹ 5 more rows

ei_synthetic(n = 10, p = 2, n_x = 3)
#> EI Specification
#> • Predictors: `x1`, `x2`, and `x3`
#> • Outcome: `y`
#> • Covariates: `z1` and `z2`
#> # A tibble: 10 × 6
#>       y    x1    x2     x3      z1     z2
#>   <dbl> <dbl> <dbl>  <dbl>   <dbl>  <dbl>
#> 1 0.358 0.489 0.339 0.172   0.0777 -0.190
#> 2 0.614 0.624 0.220 0.156  -0.169  -0.256
#> 3 0.353 0.616 0.292 0.0921 -0.598  -0.457
#> 4 0.193 0.583 0.365 0.0525  0.0779  0.479
#> 5 0.614 0.526 0.150 0.324   1.07    0.382
#> # ℹ 5 more rows

# Manual hyperparameters: x2 dominant and z1, z2 very correlated
ei_synthetic(n = 10, x = c(1, 95, 4), z = c(10, 9.999))
#> EI Specification
#> • Predictors: `x1`, `x2`, and `x3`
#> • Outcome: `y`
#> • Covariates: `z1` and `z2`
#> # A tibble: 10 × 6
#>       y      x1    x2     x3     z1     z2
#>   <dbl>   <dbl> <dbl>  <dbl>  <dbl>  <dbl>
#> 1 0.549 0.00275 0.978 0.0195 -2.48  -2.45 
#> 2 0.396 0.00315 0.948 0.0486 -1.40  -1.44 
#> 3 0.955 0.00279 0.926 0.0709  3.17   3.12 
#> 4 0.717 0.0109  0.957 0.0317  0.630  0.654
#> 5 0.505 0.00447 0.960 0.0359 -4.79  -4.79 
#> # ℹ 5 more rows

# Condition on provided x but not z
data(elec_1968)
ei_synthetic(
    x = cbind(elec_1968$pop_white, 1 - elec_1968$pop_white),
    p = 5,
    b_loc = c(0.3, 0.9),
    b_cov = matrix(c(0.02, 0.016, 0.016, 0.2), nrow=2)
)
#> EI Specification
#> • Predictors: `x1` and `x2`
#> • Outcome: `y`
#> • Covariates: `z1`, `z2`, `z3`, `z4`, and `z5`
#> # A tibble: 1,143 × 8
#>       y    x1     x2      z1       z2       z3      z4      z5
#>   <dbl> <dbl>  <dbl>   <dbl>    <dbl>    <dbl>   <dbl>   <dbl>
#> 1 0.519 0.716 0.284   0.171  -0.0160   0.215    0.0210  0.840 
#> 2 0.319 0.819 0.181  -0.869  -0.639   -0.384   -0.355  -0.0614
#> 3 0.518 0.538 0.462  -0.446  -0.504    0.00911 -0.118  -0.815 
#> 4 0.302 0.721 0.279  -0.428  -0.527   -0.108   -0.0846 -0.682 
#> 5 0.151 0.976 0.0241 -0.0945  0.00858  0.0118   0.185  -0.202 
#> # ℹ 1,138 more rows
```
