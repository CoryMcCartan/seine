# Estimate Riesz representer for ecological inference

Fits a penalized Riesz regression for ecological inference, allowing for
overall estimates of conditional means using
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md).

## Usage

``` r
ei_riesz(x, ..., weights, penalty, scale = TRUE)

# S3 method for class 'formula'
ei_riesz(formula, data, total, weights, penalty, scale = TRUE, ...)

# S3 method for class 'ei_spec'
ei_riesz(x, weights, penalty, scale = TRUE, ...)

# S3 method for class 'data.frame'
ei_riesz(x, z, total, weights, penalty, scale = TRUE, ...)

# S3 method for class 'matrix'
ei_riesz(x, z, total, weights, penalty, scale = TRUE, ...)

# Default S3 method
ei_riesz(x, ...)
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

- penalty:

  The ridge penalty (a non-negative scalar), which must be specified.
  Recommended value is the same penalty used in
  [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md),
  which is stored in the `penalty` entry of the fitted model object.

- scale:

  If `TRUE`, scale covariates `z` to have unit variance.

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

- total:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A variable containing the total number of observations in each
  aggregate unit. For example, the column containing the total number of
  voters. Required by default.

- z:

  When `x` is a **data frame** or **matrix**, `w` are any covariates,
  specified as:

  - A **data frame** with numeric columns.

  - A **matrix**

  These are shifted to have mean zero. If `scale=TRUE` (the default),
  they are also scaled to have unit variance.

## Value

An `ei_riesz` object.

## Details

The regression is calculated using the singular value decomposition.

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

# Recommended: get ridge penalty from ei_ridge()
spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(pop_urban, farm))
m = ei_ridge(spec)

ei_riesz(spec, penalty = m$penalty)
#> A Riesz representer with 3 groups and 1143 observations
#> Fit with penalty = 1.01453e-08

rr = ei_riesz(~ vap_white + vap_black + vap_other | pop_urban + farm,
              data = elec_1968, total = pres_total, penalty = m$penalty)
summary(rr)
#> Second moment of representer:
#>    vap_white    vap_black    vap_other 
#>     4.461456    48.486776 31468.779992 
#> 
#> Second moment of representer (leave-one-out):
#>   vap_white   vap_black   vap_other 
#>     4.98174    51.34839 42521.67381 

# Examine the weights and check they have mean 1
head(weights(rr, group = "vap_black"))
#> [1]  3.3584139 -0.9872354  8.7084517  4.6428913 -4.6285576  6.0225661
colMeans(weights(rr))
#> vap_white vap_black vap_other 
#> 1.0000000 1.0000000 0.9999999 

mean(elec_1968$pres_ind_wal * weights(rr, "vap_white"))
#> [1] 0.3633165
```
