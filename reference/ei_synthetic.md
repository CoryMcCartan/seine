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
#> 1 0.282 0.833 0.167
#> 2 0.427 0.771 0.229
#> 3 0.709 0.573 0.427
#> 4 0.563 0.568 0.432
#> 5 0.559 0.649 0.351
#> # ℹ 5 more rows

ei_synthetic(n = 10, p = 2, n_x = 3)
#> EI Specification
#> • Predictors: `x1`, `x2`, and `x3`
#> • Outcome: `y`
#> • Covariates: `z1` and `z2`
#> # A tibble: 10 × 6
#>       y    x1    x2    x3     z1     z2
#>   <dbl> <dbl> <dbl> <dbl>  <dbl>  <dbl>
#> 1 0.719 0.247 0.528 0.225 -0.811 -0.443
#> 2 0.483 0.523 0.244 0.234 -0.217  0.141
#> 3 0.508 0.411 0.443 0.147  0.231  0.172
#> 4 0.634 0.490 0.203 0.307  0.643  0.546
#> 5 0.793 0.343 0.378 0.279 -0.471 -0.594
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
#> 1 0.448 0.0129  0.956 0.0311 -0.238 -0.251
#> 2 0.661 0.00965 0.960 0.0308  0.504  0.482
#> 3 0.356 0.0188  0.963 0.0184  3.26   3.30 
#> 4 0.629 0.00670 0.947 0.0464 -0.306 -0.308
#> 5 0.591 0.00573 0.952 0.0420 -3.67  -3.69 
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
#>       y    x1     x2      z1      z2       z3      z4     z5
#>   <dbl> <dbl>  <dbl>   <dbl>   <dbl>    <dbl>   <dbl>  <dbl>
#> 1 0.459 0.716 0.284  -0.201   0.290  -0.00837 -0.160   0.140
#> 2 0.411 0.819 0.181   0.456  -0.572  -0.393   -0.355  -0.322
#> 3 0.120 0.538 0.462  -0.0575 -0.0419  0.234    0.215   0.364
#> 4 0.493 0.721 0.279   0.182  -0.248  -0.524   -0.0946 -0.536
#> 5 0.313 0.976 0.0241  0.765  -0.436   0.705    0.419   0.643
#> # ℹ 1,138 more rows
```
