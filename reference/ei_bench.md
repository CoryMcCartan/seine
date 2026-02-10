# Benchmark sensitivity parameters from observed covariates

Produces a table of benchmark values for `c_outcome` and `c_predictor`
in [`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md)
for each covariate, following the methodology of Chernozhukov et al.
(2024).

## Usage

``` r
ei_bench(spec, subset = NULL, contrast = NULL)
```

## Arguments

- spec:

  An [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

- subset:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  An optional indexing vector describing the subset of units over which
  to calculate estimates.

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

# with preprocessed covariates
spec = ei_spec(
    data = elec_1968,
    predictors = vap_white:vap_other,
    outcome = pres_ind_wal,
    total = pres_total,
    covariates = c(educ_elem, pop_urban, farm),
    preproc = ~ model.matrix(~ .^2 - 1, data = .x)
)
ei_bench(spec)
#> # A tibble: 9 × 7
#>   covariate predictor outcome      c_outcome c_predictor confounding est_chg
#>   <chr>     <chr>     <chr>            <dbl>       <dbl>       <dbl>   <dbl>
#> 1 educ_elem vap_white pres_ind_wal    0.193        0.797      0.422   0.0308
#> 2 educ_elem vap_black pres_ind_wal    0.193        1         -0.0493 -0.0187
#> 3 educ_elem vap_other pres_ind_wal    0.193        0.113     -0.549  -1.22  
#> 4 pop_urban vap_white pres_ind_wal    0.0860       0.143     -0.586  -0.0152
#> 5 pop_urban vap_black pres_ind_wal    0.0860       0.210      0.526   0.0681
#> 6 pop_urban vap_other pres_ind_wal    0.0860       0         -1      -0.359 
#> 7 farm      vap_white pres_ind_wal    0.254        0.700     -0.294  -0.0237
#> 8 farm      vap_black pres_ind_wal    0.254        0.617      0.486   0.160 
#> 9 farm      vap_other pres_ind_wal    0.254        1          0.162   1.08  
ei_bench(spec, subset = pop_urban > 0.5)
#> # A tibble: 9 × 7
#>   covariate predictor outcome      c_outcome c_predictor confounding  est_chg
#>   <chr>     <chr>     <chr>            <dbl>       <dbl>       <dbl>    <dbl>
#> 1 educ_elem vap_white pres_ind_wal    0.193        0.797      0.181   0.0132 
#> 2 educ_elem vap_black pres_ind_wal    0.193        1         -0.128  -0.0486 
#> 3 educ_elem vap_other pres_ind_wal    0.193        0.113     -0.213  -0.472  
#> 4 pop_urban vap_white pres_ind_wal    0.0860       0.143      1       0.0493 
#> 5 pop_urban vap_black pres_ind_wal    0.0860       0.210     -1      -0.159  
#> 6 pop_urban vap_other pres_ind_wal    0.0860       0         -1      -0.982  
#> 7 farm      vap_white pres_ind_wal    0.254        0.700      0.0289  0.00233
#> 8 farm      vap_black pres_ind_wal    0.254        0.617      0.0813  0.0268 
#> 9 farm      vap_other pres_ind_wal    0.254        1          0.0823  0.545  

# with contrasts
spec = ei_spec(elec_1968, vap_white:vap_other, pres_rep_nix:pres_ind_wal,
               total = pres_total, covariates = c(educ_elem, pop_urban, farm))
ei_bench(spec, contrast = list(predictor = c(1, -1, 0)))
#> # A tibble: 6 × 7
#>   covariate predictor         outcome c_outcome c_predictor confounding  est_chg
#>   <chr>     <chr>             <chr>       <dbl>       <dbl>       <dbl>    <dbl>
#> 1 educ_elem vap_white - vap_… pres_r…   0.0337        1          -0.548 -0.0748 
#> 2 educ_elem vap_white - vap_… pres_i…   0.0852        1           0.855  0.175  
#> 3 pop_urban vap_white - vap_… pres_r…   0.00543       0.112       0.134  0.00327
#> 4 pop_urban vap_white - vap_… pres_i…   0.0672        0.112      -0.956 -0.0773 
#> 5 farm      vap_white - vap_… pres_r…   0.0255        0.141       0.217  0.0127 
#> 6 farm      vap_white - vap_… pres_i…   0.123         0.141      -0.478 -0.0580 
ei_bench(spec, contrast = list(outcome = c(1, -1)))
#> # A tibble: 9 × 7
#>   covariate predictor outcome          c_outcome c_predictor confounding est_chg
#>   <chr>     <chr>     <chr>                <dbl>       <dbl>       <dbl>   <dbl>
#> 1 educ_elem vap_white pres_rep_nix - …   0            0.601       -1     -0.0455
#> 2 educ_elem vap_black pres_rep_nix - …   0            1            1      0.204 
#> 3 educ_elem vap_other pres_rep_nix - …   0            1           -1     -0.125 
#> 4 pop_urban vap_white pres_rep_nix - …   0.00188      0.145        1      0.0197
#> 5 pop_urban vap_black pres_rep_nix - …   0.00188      0.104       -1     -0.0609
#> 6 pop_urban vap_other pres_rep_nix - …   0.00188      0.0365      -0.822 -0.544 
#> 7 farm      vap_white pres_rep_nix - …   0.00722      0.260        0.749  0.0238
#> 8 farm      vap_black pres_rep_nix - …   0.00722      0.102       -0.577 -0.0468
#> 9 farm      vap_other pres_rep_nix - …   0.00722      0.226       -0.691 -2.05  
```
