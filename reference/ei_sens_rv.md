# Robustness values for ecological inference

The robustness value is the minimum bound for both `c_outcome` and
`c_predictor` in
[`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md) such
that the bias bound is a certain value. For example, if the provided
bias bound is 0.5, meaning a bias of magnitude 0.5 would be considered
problematic, then the robustness value is the minimum amount of
confounding of outcome and predictor (more specifically, the Riesz
representer) that would lead to bias of magnitude 0.5.

## Usage

``` r
ei_sens_rv(est, bias_bound, confounding = 1)
```

## Arguments

- est:

  A set of estimates from
  [`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md) using
  both regression and Riesz representer.

- bias_bound:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A bias bound: an amount of bias which is considered substantial.
  Evaluated in the context of `est`, so that one can to refer to
  `std.error` and `estimate` as needed.

- confounding:

  The confounding parameter (\\\rho\\), which must be between 0 and 1
  (the adversarial worst-case).

## Value

A data frame of the same format as `est`, but with a new `rv` column
containing the robustness values.

## References

Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V.
(2024). *Long story short: Omitted variable bias in causal machine
learning* (No. w30302). National Bureau of Economic Research.

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
               total = pres_total, covariates = c(state, pop_urban, farm))
m = ei_ridge(spec)
rr = ei_riesz(spec, penalty = m$penalty)
est = ei_est(m, rr, spec)

ei_sens_rv(est, 0.1) # how much confounding for bias of 0.1
#> # A tibble: 3 × 5
#>   predictor outcome      estimate std.error      rv
#>   <chr>     <chr>           <dbl>     <dbl>   <dbl>
#> 1 vap_white pres_ind_wal    0.387    0.0246 0.331  
#> 2 vap_black pres_ind_wal    0.490    0.0505 0.108  
#> 3 vap_other pres_ind_wal   -1.19     0.528  0.00491
ei_sens_rv(est, 2 * std.error) # how much confounding for bias of 2 SE
#> # A tibble: 3 × 5
#>   predictor outcome      estimate std.error     rv
#>   <chr>     <chr>           <dbl>     <dbl>  <dbl>
#> 1 vap_white pres_ind_wal    0.387    0.0246 0.180 
#> 2 vap_black pres_ind_wal    0.490    0.0505 0.109 
#> 3 vap_other pres_ind_wal   -1.19     0.528  0.0506

# How much confounding to equalize all estimates (no polarization)
y_avg = weighted.mean(elec_1968$pres_ind_wal, elec_1968$pres_total)
ei_sens_rv(est, estimate - y_avg)
#> # A tibble: 3 × 5
#>   predictor outcome      estimate std.error     rv
#>   <chr>     <chr>           <dbl>     <dbl>  <dbl>
#> 1 vap_white pres_ind_wal    0.387    0.0246 0.164 
#> 2 vap_black pres_ind_wal    0.490    0.0505 0.154 
#> 3 vap_other pres_ind_wal   -1.19     0.528  0.0725

# Extract as matrix
as.matrix(ei_sens_rv(est, 0.2), "rv")
#>            outcome
#> predictor   pres_ind_wal
#>   vap_white  0.545120562
#>   vap_black  0.203654212
#>   vap_other  0.009796171

# Works for contrasts as well
est = ei_est(m, rr, spec, contrast = list(predictor=c(1, -1, 0)))
ei_sens_rv(est, estimate) # how much to eliminate disparity
#> # A tibble: 1 × 5
#>   predictor             outcome      estimate std.error     rv
#>   <chr>                 <chr>           <dbl>     <dbl>  <dbl>
#> 1 vap_white - vap_black pres_ind_wal   -0.102    0.0457 0.0935
```
