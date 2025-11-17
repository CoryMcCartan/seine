# Benchmark sensitivity parameters from observed covariates

Produces a table of benchmark values for `c_outcome` and `c_predictor`
in [`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md)
for each covariate, following the methodology of Chernozhukov et al.
(2024).

## Usage

``` r
ei_bench(spec, preproc = NULL, subset = NULL)
```

## Arguments

- spec:

  An [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

- preproc:

  An optional function which takes in a `ei_spec` object (`spec` with
  one covariate removed) and returns a modified object that includes
  modified object. Useful to apply any preprocessing, such as a basis
  transformation, as part of the benchmarking process.

- subset:

  Passed on to
  [`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md).

## References

Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V.
(2024). *Long story short: Omitted variable bias in causal machine
learning* (No. w30302). National Bureau of Economic Research.

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
               total = pres_total, covariates = c(educ_elem, pop_urban, farm))

ei_bench(spec)
#> # A tibble: 9 × 7
#>   covariate predictor outcome      c_outcome c_predictor confounding est_chg
#>   <chr>     <chr>     <chr>            <dbl>       <dbl>       <dbl>   <dbl>
#> 1 educ_elem vap_white pres_ind_wal    0.0859      0.609       0.686   0.0272
#> 2 educ_elem vap_black pres_ind_wal    0.0859      1          -0.827  -0.148 
#> 3 educ_elem vap_other pres_ind_wal    0.0859      1           0.0880  0.422 
#> 4 pop_urban vap_white pres_ind_wal    0.0680      0.151      -0.940  -0.0195
#> 5 pop_urban vap_black pres_ind_wal    0.0680      0.104       0.863   0.0578
#> 6 pop_urban vap_other pres_ind_wal    0.0680      0.0560      0.245   0.322 
#> 7 farm      vap_white pres_ind_wal    0.125       0.275      -0.563  -0.0203
#> 8 farm      vap_black pres_ind_wal    0.125       0.102       0.420   0.0377
#> 9 farm      vap_other pres_ind_wal    0.125       0.289       0.471   1.72  

# preprocess to add all 2-way interactions
ei_bench(spec, preproc = function(s) {
    z_cols = match(attr(s, "ei_z"), names(s))
    s_out = s[-z_cols]
    z_new = model.matrix(~ .^2 - 1, data = s[z_cols])
    s_out = cbind(s_out, z_new)
    ei_spec(s_out, vap_white:vap_other, pres_ind_wal,
            total = attr(s, "ei_n"), covariates = colnames(z_new))
})
#> # A tibble: 9 × 7
#>   covariate predictor outcome      c_outcome c_predictor confounding est_chg
#>   <chr>     <chr>     <chr>            <dbl>       <dbl>       <dbl>   <dbl>
#> 1 educ_elem vap_white pres_ind_wal    0.193        0.797      0.422   0.0308
#> 2 educ_elem vap_black pres_ind_wal    0.193        1         -0.0493 -0.0187
#> 3 educ_elem vap_other pres_ind_wal    0.193        0.113     -0.549  -1.22  
#> 4 pop_urban vap_white pres_ind_wal    0.0860       0.143     -0.586  -0.0152
#> 5 pop_urban vap_black pres_ind_wal    0.0860       0.210      0.526   0.0681
#> 6 pop_urban vap_other pres_ind_wal    0.0860      -0.345     -1      -0.359 
#> 7 farm      vap_white pres_ind_wal    0.254        0.700     -0.294  -0.0237
#> 8 farm      vap_black pres_ind_wal    0.254        0.617      0.486   0.160 
#> 9 farm      vap_other pres_ind_wal    0.254        1          0.162   1.08  
```
