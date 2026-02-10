# Produce local ecological estimates

Projects predictions from a fitted regression model onto the accounting
constraint using a provided residual covariance matrix. This ensures
that each set of local estimates satisfies the accounting identity.
Local estimates may be truncated to variable bounds.

## Usage

``` r
ei_est_local(
  regr,
  data,
  total,
  b_cov,
  contrast = NULL,
  bounds = regr$blueprint$bounds,
  sum_one = NULL,
  conf_level = 0.95,
  regr_var = TRUE,
  unimodal = TRUE
)

# S3 method for class 'ei_est_local'
as.array(x, ...)
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

- total:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  A variable containing the total number of observations in each
  aggregate unit. For example, the column containing the total number of
  voters. Required if `data` is not an
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

- b_cov:

  A covariance matrix to use in projecting the local estimates `b` onto
  the accounting constraint, such as one estimated with
  [`ei_local_cov()`](https://corymccartan.com/seine/reference/ei_local_cov.md).
  When `b_cov` is a single number, it is used to specify the pairwise
  correlation between the local estimates across predictors, with the
  same correlation structure as the regression residuals assumed for `b`
  within each predictor. For example, set `b_cov=0` to assume no
  correlation in `b` across predictors, or set `b_cov=1` to assume
  perfect correlation in `b` across predictors, corresponding to a
  (local) neighborhood model. When there are multiple outcome variables
  and `b_cov` is a matrix with entries for each predictor, it will be
  applied identically to each outcome. Alternatively, a matrix with
  entries for each predictor-outcome combination may be provided, with
  entries in the order (Y1\|X1, Y1\|X2, ..., Y2\|X1, Y2\|X2, ...).

- contrast:

  If provided, a list containing entries `predictor` and `outcome`, each
  containing a contrast vector. If only one of `predictor` or `outcome`
  is provided, the contrast will be calculated for all levels of the
  other variable. For example `list(predictor = c(1, -1, 0))` will
  calculate the difference in each outcome between the first and second
  predictor groups; `list(outcome = c(1, -1))` will calculate the
  difference between the two outcomes for each predictor group; and
  `list(predictor = c(1, -1, 0), outcome = c(1, -1))` will calculate the
  difference in differences.

- bounds:

  A vector `c(min, max)` of bounds for the outcome, to which the local
  estimates will be truncated. In general, truncation will lead to
  violations of the accounting identity. If `bounds = NULL`, they will
  be inferred from the outcome variable: if it is contained within
  \\\[0, 1\]\\, for instance, then the bounds will be `c(0, 1)`. Setting
  `bounds = FALSE` forces unbounded estimates. The default uses the
  `bounds` attribute of `regr`, if available, or infers from the outcome
  variable otherwise. Note that bounds are currently not applied if
  `contrast` is provided.

- sum_one:

  If `TRUE`, the outcome variables are constrained to sum to one. Can
  only apply when `bounds` are enforced and there is more than one
  outcome variable. If `NULL`, infers `sum_one = TRUE` when the bounds
  are `c(0, 1)` the outcome variables sum to 1.

- conf_level:

  A numeric specifying the level for confidence intervals. If `FALSE`,
  no confidence intervals are calculated. For `regr` arguments from
  [`ei_wrap_model()`](https://corymccartan.com/seine/reference/ei_wrap_model.md),
  confidence intervals will not incorporate uncertainty in the
  prediction itself, just the residual. This will trigger a warning
  periodically.

- regr_var:

  If `TRUE`, incorporate uncertainty from the regression model when
  calculating confidence intervals. Only applies when `regr` is fitted
  with
  [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md),
  and requires that function be called with `vcov = TRUE`.

- unimodal:

  If `TRUE`, assume a unimodal residual distribution. Reduces width of
  confidence intervals by a factor of 2/3.

- x:

  An object of class `ei_est_local`

- ...:

  Additional arguments (ignored)

## Value

A data frame with estimates. The `.row` column in the output corresponds
to the observation index in the input. The `weight` column contains the
product of the predictor variable and total for each observation. Taking
a weighted average of the estimate against this column will produce a
global estimate. It has class `ei_est_local`.

## Details

Local estimates are produced jointly across outcome variables. When
bounds are applied, unless `sum_one = TRUE`, the estimates for each
observation may not satisfy logical constraints, including the
accounting identity.

Projections are done obliquely in accordance with `b_cov` via quadratic
programming. Occasionally, the quadratic program may be infeasible due
to the specific data, features of `b_cov`, or numerical errors. Various
relaxations of the accounting identity and `b_cov` are attempted in
these cases; indices where relaxations of `b_cov` were used are stored
in the `proj_relax` attribute of the output, and indices of infeasible
projections are stored in the `proj_misses` attribute.

## Methods (by generic)

- `as.array(ei_est_local)`: Format estimates an array with dimensions
  `<rows>*<predictors>*<outcomes>`. Does not work if the object has been
  sorted.

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm))

m = ei_ridge(spec)

ei_est_local(m, spec, b_cov = 0, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.99)
#> # A tibble: 13,716 × 8
#>     .row predictor outcome      weight estimate std.error conf.low conf.high
#>    <int> <chr>     <chr>         <dbl>    <dbl>     <dbl>    <dbl>     <dbl>
#>  1     1 vap_white pres_dem_hum  5877.   0.0507   0.0260   0          0.224 
#>  2     2 vap_white pres_dem_hum 16131.   0.0146   0.0137   0          0.106 
#>  3     3 vap_white pres_dem_hum  4872.   0.0105   0.0454   0          0.313 
#>  4     4 vap_white pres_dem_hum  3566.   0        0.0236   0          0.157 
#>  5     5 vap_white pres_dem_hum  8801.   0.0263   0.00161  0.0190     0.0371
#>  6     6 vap_white pres_dem_hum  1698.   0.110    0.0681   0          0.564 
#>  7     7 vap_white pres_dem_hum  4970.   0        0.0364   0          0.243 
#>  8     8 vap_white pres_dem_hum 22844.   0.0591   0.0152   0.00707    0.161 
#>  9     9 vap_white pres_dem_hum  7731.   0        0.0304   0          0.178 
#> 10    10 vap_white pres_dem_hum  5259.   0.0409   0.00651  0.00915    0.0843
#> # ℹ 13,706 more rows

b_cov = ei_local_cov(m, spec)
e_orth = ei_est_local(m, spec, b_cov = 0, bounds = c(0, 1), sum_one = TRUE)
e_nbhd = ei_est_local(m, spec, b_cov = 1, bounds = c(0, 1), sum_one = TRUE)
e_rcov = ei_est_local(m, spec, b_cov = b_cov, bounds = c(0, 1), sum_one = TRUE)
# average interval width
c(
    e_orth = mean(e_orth$conf.high - e_orth$conf.low),
    e_nbhd = mean(e_nbhd$conf.high - e_nbhd$conf.low),
    e_rcov = mean(e_rcov$conf.high - e_rcov$conf.low)
)
#>    e_orth    e_nbhd    e_rcov 
#> 0.3626324 0.3123206 0.4427580 
```
