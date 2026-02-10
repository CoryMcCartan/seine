# Test the coarsening at random (CAR) assumption

The CAR assumption, which is required for accurate estimates with
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md),
implies that the conditional expectation function (CEF) of the aggregate
outcomes will take a certain partially linear form. This function tests
that implication by comparing a fully nonparametric estimate of the CEF
to one in the partially linear form implied by CAR, and comparing their
goodness-of-fit. It is **very important** that `preproc` be used in the
`spec` to perform a rich basis transformation of the covariates and
predictors; a missing `preproc` will lead to a warning.

## Usage

``` r
ei_test_car(spec, weights, iter = 1000, use_chisq = FALSE)
```

## Arguments

- spec:

  An `ei_spec` object created with
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md).
  The object should use the `preproc` argument to specify a rich basis
  expansion of the covariates and predictors.

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

- iter:

  The number of permutations to use when estimating the null
  distribution. Ignored when `use_chisq = TRUE`.

- use_chisq:

  If `TRUE`, use the asymptotic chi-squared distribution for the Wald
  test statistic instead of conducting a permutation test. Only
  appropriate for larger sample sizes (Helwig 2022 recommends at least
  200 when a single predictor is used).

## Value

A tibble with one row per outcome variable and columns describing the
test results. The `p.value` column contains the p-values for the test.
P-values are not adjusted by default; pass them to
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) if desired.

## Details

The test is a Kennedy-Cade (1996) style permutation test on a Wald
statistic for the coefficients not included in the "reduced" model that
would be fit by
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md).
The test is carried out by fitting a regression on a fully
basis-expanded combination of covariates and predictors, and calculating
a Wald statistic for the

## References

Helwig, N. E. (2022). Robust Permutation Tests for Penalized Splines.
*Stats*, 5(3), 916-933.

Kennedy, P. E., & Cade, B. S. (1996). Randomization tests for multiple
regression. *Communications in Statistics-Simulation and Computation*,
25(4), 923-936.

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

## Examples

``` r
data(elec_1968)

# basis expansion: poly() with degree=2 not recommended in practice
preproc = if (requireNamespace("bases", quietly = TRUE)) {
    ~ bases::b_bart(.x, trees = 100)
} else {
    ~ poly(as.matrix(.x), degree=2, simple=TRUE)
}

spec = ei_spec(
    data = elec_1968,
    predictors = vap_white:vap_other,
    outcome = pres_dem_hum:pres_abs,
    total = pres_total,
    covariates = c(pop_city:pop_rural, farm:educ_coll, starts_with("inc_")),
    preproc = preproc
)

ei_test_car(spec, iter=19) # use a larger number in practice
#> # A tibble: 4 × 4
#>   outcome          W    df p.value
#>   <chr>        <dbl> <int>   <dbl>
#> 1 pres_dem_hum  300.   105    0.05
#> 2 pres_rep_nix  257.   105    0.05
#> 3 pres_ind_wal  357.   105    0.05
#> 4 pres_abs      206.   105    0.05
```
