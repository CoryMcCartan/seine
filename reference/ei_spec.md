# Specify an ecological inference problem

Uses tidy-select syntax to specify outcomes, predictors, and covariates.
The result of this function can be passed directly into
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) or
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md), or
plotted with
[`plot()`](https://corymccartan.com/seine/reference/plot.ei_spec.md).

## Usage

``` r
ei_spec(
  data,
  predictors,
  outcome,
  total,
  covariates = NULL,
  preproc = NULL,
  strip = FALSE
)

# S3 method for class 'ei_spec'
weights(object, normalize = TRUE, ...)
```

## Arguments

- data:

  A data frame.

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

- total:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A variable containing the total number of observations in each
  aggregate unit. For example, the column containing the total number of
  voters. Required by default.

- covariates:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Covariates.

- preproc:

  An optional function which takes in a data frame of covariates and
  returns a modeling-ready numeric matrix of covariates. Useful to apply
  any preprocessing, such as a basis transformation, as part of the
  estimation process. Passed to
  [`rlang::as_function()`](https://rlang.r-lib.org/reference/as_function.html),
  and so supports `purrr`-style lambda functions. This function is
  called once when forming the `ei_spec` object so that the same
  processing is applied during all estimation steps. The function may
  also be re-used by other functions, like
  [`ei_bench()`](https://corymccartan.com/seine/reference/ei_bench.md)
  and
  [`ei_test_car()`](https://corymccartan.com/seine/reference/ei_test_car.md).
  The default is a call to
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) with the
  formula `~ 0 + .`, i.e., all covariates and no predictors.

- strip:

  Whether to strip common prefixes from column names within each group.
  For example, columns named `vap_white`, `vap_black`, and `vap_hisp`
  would be renamed `white`, `black` and `other` in the model and output.

- object:

  An ei_spec object.

- normalize:

  If `TRUE`, normalize the totals to have mean 1.

- ...:

  Additional arguments (ignored).

## Value

An `ei_spec` object, which is a data frame with additional attributes
recording `predictors`, `outcomes`, `total`, and `covariates`.

## Details

The function is lightweight and does not perform any checking of the
arguments, bounds, sum constraints, etc. All of these checks are
performed by functions that use `ei_spec` objects.

## Methods (by generic)

- `weights(ei_spec)`: Extract the totals from a specification

## Examples

``` r
data(elec_1968)
ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs, pres_total)
#> EI Specification
#> • Predictors: `vap_white`, `vap_black`, and `vap_other`
#> • Outcome: `pres_dem_hum`, `pres_rep_nix`, `pres_ind_wal`, and `pres_abs`
#> • Covariates: none
#> # A tibble: 1,143 × 7
#>   vap_white vap_black vap_other pres_dem_hum pres_rep_nix pres_ind_wal pres_abs
#>       <dbl>     <dbl>     <dbl>        <dbl>        <dbl>        <dbl>    <dbl>
#> 1     0.761    0.237   0.00173        0.199        0.0773        0.711  0.0122 
#> 2     0.860    0.137   0.00306        0.105        0.115         0.764  0.0161 
#> 3     0.610    0.389   0.000808       0.242        0.0489        0.687  0.0218 
#> 4     0.783    0.216   0.00106        0.141        0.0571        0.799  0.00290
#> 5     0.981    0.0181  0.000757       0.0375       0.222         0.727  0.0134 
#> # ℹ 1,138 more rows

# basis expansion
if (requireNamespace("bases", quietly = TRUE)) {
    spec = ei_spec(
        data = elec_1968,
        predictors = vap_white:vap_other,
        outcome = pres_dem_hum:pres_abs,
        total = pres_total,
        covariates = c(pop_city:pop_rural, farm:educ_coll, starts_with("inc_")),
        preproc = ~ bases::b_bart(.x, trees = 500)
    )
}
```
