# Low-level implementations of `ei_ridge()` and `ei_riesz()`

No checks are performed on the inputs. Use of
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) and
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md) is
strongly recommended unless many regressions must be fit, e.g., within a
tight loop. Only works for a single outcome, i.e., `y` must be a vector,
not a matrix.

## Usage

``` r
ei_ridge_impl(
  x,
  y,
  z,
  weights = rep(1, nrow(x)),
  bounds = c(-Inf, Inf),
  penalty = NULL,
  vcov = TRUE
)

ei_riesz_impl(x, z, total, weights = rep(1, nrow(x)), penalty)
```

## Arguments

- x:

  A matrix of predictors

- y:

  A vector of outcomes

- z:

  A matrix of covariates

- weights:

  A vector of estimation weights

- bounds:

  A vector `c(min, max)` of bounds for the outcome.

- penalty:

  The ridge penalty (a non-negative scalar), which must be specified for
  `ei_riesz_impl()` but can be automatically estimated with
  `ei_ridge_impl()` by providing `penalty=NULL`.

- vcov:

  If `TRUE`, calculate and return the covariance matrix of the estimated
  coefficients. Ignored when `bounds` are provided.

- total:

  A vector of total observations per unit.

## Value

A list with model components.
