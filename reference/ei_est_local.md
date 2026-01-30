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
  r_cov = NULL,
  bounds = regr$blueprint$bounds,
  sum_one = NULL,
  conf_level = FALSE,
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

- r_cov:

  A covariance matrix of the residuals to use in projecting the local
  estimates onto the accounting constraint, or a list of matrices, one
  for each outcome variable. Defaults to the identity matrix scaled by
  the residual variance of `regr`, corresponding to orthogonal
  projection. Set `r_cov=1` to use a degenerate covariance matrix
  corresponding to a (local) neighborhood model. When there are multiple
  outcome variables and `r_cov` is a matrix, it will be applied
  identically to each outcome.

- bounds:

  A vector `c(min, max)` of bounds for the outcome, to which the local
  estimates will be truncated. In general, truncation will lead to
  violations of the accounting identity. If `bounds = NULL`, they will
  be inferred from the outcome variable: if it is contained within
  \\\[0, 1\]\\, for instance, then the bounds will be `c(0, 1)`. Setting
  `bounds = FALSE` forces unbounded estimates. The default uses the
  `bounds` attribute of `regr`, if available, or infers from the outcome
  variable otherwise.

- sum_one:

  If `TRUE`, the outcome variables are constrained to sum to one. Can
  only apply when `bounds` are enforced and there is more than one
  outcome variable. The default `NULL` infers `sum_one = TRUE` when the
  bounds are `c(0, 1)` the outcome variables sum to 1.

- conf_level:

  A numeric specifying the level for confidence intervals. If `FALSE`
  (the default), no confidence intervals are calculated. For `regr`
  arguments from
  [`ei_wrap_model()`](https://corymccartan.com/seine/reference/ei_wrap_model.md),
  confidence intervals will not incorporate uncertainty in the
  prediction itself, just the residual. This will trigger a warning
  periodically.

- unimodal:

  If `TRUE`, assume a unimodal residual distribution. Improves width of
  confidence intervals by a factor of 4/9.

- x:

  An object of class `ei_est_local`

- ...:

  Additional arguments (ignored)

## Value

A data frame with estimates. The `.row` column in the output corresponds
to the observation index in the input. It has class `ei_est_local`,
supporting several methods.

## Details

Local estimates are produced independently for each outcome variable.
Truncation to bounds, if used, will in general lead to estimates that do
not satisfy the accounting identity.

## Methods (by generic)

- `as.array(ei_est_local)`: Format estimates an array with dimensions
  `<rows>*<predictors>*<outcomes>`. Does not work if the object has been
  sorted.

## Examples

``` r
if (FALSE) { # \dontrun{
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm))

m = ei_ridge(spec)

ei_est_local(m, spec, conf_level = 0.95)
suppressWarnings(ei_est_local(m, spec, bounds=c(0.01, 0.2)))
} # }
```
