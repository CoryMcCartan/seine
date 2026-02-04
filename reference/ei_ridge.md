# Fit an ecological inference regression model

Fits a penalized regression model for ecological inference, allowing for
overall and unit-level estimates of conditional means using
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md).

## Usage

``` r
ei_ridge(
  x,
  ...,
  weights,
  bounds = FALSE,
  sum_one = FALSE,
  penalty = NULL,
  scale = TRUE,
  vcov = TRUE
)

# S3 method for class 'formula'
ei_ridge(
  formula,
  data,
  weights,
  bounds = FALSE,
  sum_one = FALSE,
  penalty = NULL,
  scale = TRUE,
  vcov = TRUE,
  ...
)

# S3 method for class 'ei_spec'
ei_ridge(
  x,
  weights,
  bounds = FALSE,
  sum_one = FALSE,
  penalty = NULL,
  scale = TRUE,
  vcov = TRUE,
  ...
)

# S3 method for class 'data.frame'
ei_ridge(
  x,
  y,
  z,
  weights,
  bounds = FALSE,
  sum_one = FALSE,
  penalty = NULL,
  scale = TRUE,
  vcov = TRUE,
  ...
)

# S3 method for class 'matrix'
ei_ridge(
  x,
  y,
  z,
  weights,
  bounds = FALSE,
  sum_one = FALSE,
  penalty = NULL,
  scale = TRUE,
  vcov = TRUE,
  ...
)

# Default S3 method
ei_ridge(x, ...)
```

## Arguments

- x:

  Depending on the context:

  - A **data frame** of predictors.

  - A **matrix** of predictors.

  - An [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md)
    object containing the outcome, predictor, and covariates.

  Predictors must be proportions that sum to 1 across rows. You can use
  [`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
  to assist in preparing predictor variables. Covariates in an
  [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md) object
  are shifted to have mean zero. If `scale=TRUE` (the default), they are
  also scaled to have unit variance.

- ...:

  Not currently used, but required for extensibility.

- weights:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A vector of unit weights for estimation. These may be the same or
  different from the total number of observations in each aggregate unit
  (see the `total` argument to
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)).
  See the discussion below under 'Weights' for choosing this parameter.
  The default, uniform weights, makes a slightly stronger-than-necessary
  assumption about the relationship between the unit totals and the
  unknown data.

- bounds:

  A vector `c(min, max)` of bounds for the outcome. If `bounds = NULL`,
  they will be inferred from the outcome variable: if it is contained
  within \\\[0, 1\]\\, for instance, then the bounds will be `c(0, 1)`.
  The default `bounds = FALSE` uses an unbounded outcome.

- sum_one:

  If `TRUE`, the outcome variables are constrained to sum to one. Can
  only apply when `bounds` are enforced and there is more than one
  outcome variable. If `NULL`, infers `sum_one = TRUE` when the bounds
  are `c(0, 1)` the outcome variables sum to 1.

- penalty:

  The ridge penalty (a non-negative scalar). Set to `NULL` to
  automatically determine the penalty which minimizes mean-square error,
  via an efficient leave-one-out cross validation procedure. The ridge
  regression solution is \$\$\hat\beta = (X^\top X + \lambda
  I)^{-1}X^\top y,\$\$ where \\\lambda\\ is the value of `penalty`. One
  can equivalently think of the penalty as imposing a \\\mathcal{N}(0,
  \sigma^2/\lambda^2)\\ prior on the \\\beta\\. Keep in mind when
  choosing `penalty` manually that covariates in `z` are scaled to have
  mean zero and unit variance before fitting.

- scale:

  If `TRUE`, scale covariates `z` to have unit variance.

- vcov:

  If `TRUE`, calculate and return the a scaled covariance matrix of the
  estimated coefficients. When `bounds` are provided, the (scaled)
  covariance matrix for the unbounded estimate is returned as a
  conservative approximation. The covariance matrix is "scaled" because
  it does not include the residual variance. For the covariance for a
  particular outcome variable, multiply the returned `$vcov_u` by
  `sigma2` for that outcome.

- formula:

  A formula such as `y ~ x0 + x1 | z` specifying the outcome `y`
  regressed on the predictors of interest `x` and any covariates `z`.
  The predictors should form a partition, that is, `x0 + x1 = 1` for
  each observation. Users can be include more than two predictors as
  well, e.g. `pct_white + pct_black + pct_hisp + pct_other`. If there
  are just two predictors, it is acceptable to only include one in the
  formula; the other will be formed as 1 minus the provided predictor.
  Include additional covariates separated by a vertical bar `|`. These
  covariates are strongly recommended for reliable ecological inference.
  Covariates are shifted to have mean zero. If `scale=TRUE` (the
  default), they are also scaled to have unit variance.

- data:

  When a **formula** is used, `data` is a **data frame** containing both
  the predictors and the outcome.

- y:

  When `x` is a **data frame** or **matrix**, `y` is the outcome
  specified as:

  - A **data frame** with numeric columns.

  - A **matrix**

  - A numeric **vector**.

  When the outcome is a proportion, you can use
  [`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
  to assist in preparing it.

- z:

  When `x` is a **data frame** or **matrix**, `w` are any covariates,
  specified as:

  - A **data frame** with numeric columns.

  - A **matrix**

  These are shifted to have mean zero. If `scale=TRUE` (the default),
  they are also scaled to have unit variance.

## Value

An `ei_ridge` object, which supports various
[ridge-methods](https://corymccartan.com/seine/reference/ridge-methods.md).

## Details

The regression is calculated using the singular value decomposition,
which allows for efficient recalculation under different `penalty`
values as part of leave-one-out cross-validation. When `bounds` are
provided, the regression is calculated via quadratic programming, as
there is no closed-form solution. The unbounded regression is run to
select the `penalty` automatically in this case, if it is not provided.
Estimation is still efficient, though somewhat slower than in the
unbounded case. The covariance matrix of the estimates is not available
when bounds are applied.

## Weights

The weakest identification result for ecological inference makes no
assumption about the number of observations per aggregate unit (the
totals). It requires, however, weighting the estimation steps according
to the totals. This may reduce efficiency when the totals are variable
and a slightly stronger condition holds.

Specifically, if the totals are conditionally mean-independent of the
missing data (the aggregation-unit level means of the outcome within
each predictor level), given covariates, then it is appropriate to use
uniform weights in estimation, or any fixed set of weights.

In general, estimation efficiency is improved when units with larger
variance in the outcome receive less weight. Various bulit-in options
are provided by the helper functions in
[`ei_wgt()`](https://corymccartan.com/seine/reference/ei_wgt.md).

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(pop_urban, farm))
ei_ridge(spec)
#> An ecological inference model with 4 outcomes, 3 groups, and 1143 observations
#> Fit with penalty = 1.01453e-08

ei_ridge(pres_dem_hum + pres_rep_nix + pres_ind_wal + pres_abs ~
      vap_white + vap_black + vap_other | pop_urban + farm, data = elec_1968)
#> An ecological inference model with 4 outcomes, 3 groups, and 1143 observations
#> Fit with penalty = 1.01453e-08

# bounds inferred
all.equal(
  fitted(ei_ridge(spec, bounds = NULL)),
  fitted(ei_ridge(spec, bounds = 0:1))
)
#> [1] TRUE

# bounds enforced
min(fitted(ei_ridge(spec)))
#> [1] -0.1029557
min(fitted(ei_ridge(spec, bounds = 0:1)))
#> [1] 3.765682e-20
```
