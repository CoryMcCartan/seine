# Wrap another predictive model for use in `ei_est`

Stores additional data and attributes on a generic model class so that
it can be used as the `regr` argument to
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md). Given
the wide variety of model classes, there is no guarantee this function
will work. However, most model classes supporting a
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) method will work as
long as there is no transformation of the predictor variables as part of
the model formula or fitting.

## Usage

``` r
ei_wrap_model(x, data, predictors = NULL, outcome = NULL, ...)
```

## Arguments

- x:

  A model object, supporting
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) generics.

- data:

  A data frame or matrix containing the data used to fit the model, or
  an [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  object (recommended). If the latter, then the `predictors` and
  `outcome` arguments are ignored and need not be provided.

- predictors:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Predictor variables. This is the `x` variable in ecological regression
  that is of primary interest. For example, the columns containing the
  percentage of each racial group.

- outcome:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Outcome variables. This is the `y` variable in ecological regression
  that is of primary interest. For example, the columns containing the
  percentage of votes for each party.

- ...:

  Additional arguments passed to the
  [`predict()`](https://rdrr.io/r/stats/predict.html) method.

## Value

An `ei_wrapped` object, which has the information required to use the
provided `x` with
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md).

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal, pres_total,
               covariates = c(pop_urban, farm))

# Note: this is not a model recommended for valid ecological inference!
m = suppressWarnings(
    glm(pres_ind_wal ~ 0 + vap_white + vap_black + vap_other + pop_urban + farm,
        data = spec, family = "binomial")
)
m_wrap = ei_wrap_model(m, spec, type = "response")
print(m_wrap)
#> A wrapped <glm/lm> model with 1143 observations

ei_est(m_wrap, data = spec) # notice all estimates nonnegative
#> # A tibble: 3 × 6
#>   predictor outcome        estimate std.error conf.low conf.high
#>   <chr>     <chr>             <dbl>     <dbl>    <dbl>     <dbl>
#> 1 vap_white pres_ind_wal 0.339        0.0239    0.292     0.385 
#> 2 vap_black pres_ind_wal 0.741        0.0565    0.630     0.852 
#> 3 vap_other pres_ind_wal 0.00000298   0.00539  -0.0106    0.0106
```
