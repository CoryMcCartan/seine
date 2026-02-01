# Estimate ecological quantities

Produces estimates of overall conditional means from a fitted ecological
inference model or Riesz representer. If both a regression model and a
Riesz representer are provided, a debiased machine learning (DML)
estimate is produced.

## Usage

``` r
ei_est(
  regr = NULL,
  riesz = NULL,
  data,
  total,
  subset = NULL,
  contrast = NULL,
  outcome = NULL,
  conf_level = FALSE,
  use_student = TRUE
)

# S3 method for class 'ei_est'
as.matrix(x, which = "estimate", ...)

# S3 method for class 'ei_est'
vcov(object, ...)

# S3 method for class 'ei_est'
nobs(object, ...)
```

## Arguments

- regr:

  A fitted regression model, from
  [`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md),
  or another kind of regression model wrapped with
  [`ei_wrap_model()`](https://corymccartan.com/seine/reference/ei_wrap_model.md).
  If `riesz` is not provided and `regr` is an
  [`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md)
  object, then `riesz` will be set to the value of `regr` and `regr`
  will be set to `NULL`. This is so users can call this function as
  `ei_est(<riesz>, data = <data>)`.

- riesz:

  A fitted Riesz representer, from
  [`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md),
  or a matrix of Riesz weights

- data:

  The data frame, matrix, or
  [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md) object
  that was used to fit the regression or Riesz representer.

- total:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  A variable containing the total number of observations in each
  aggregate unit. For example, the column containing the total number of
  voters. Required if `data` is not an
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  object and `riesz` is not provided.

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

- outcome:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A vector or matrix of outcome variables. Only required if both `riesz`
  is provided alone (without `regr`) and `data` is not an
  [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md) object.

- conf_level:

  A numeric specifying the level for confidence intervals. If `FALSE`
  (the default), no confidence intervals are calculated. Standard errors
  are always returned.

- use_student:

  If `TRUE`, use construct confidence intervals from a Student-*t*
  distribution, which may improve coverage properties in small samples.

- x, object:

  An object of class `ei_est`

- which:

  Which column of `ei_est` to convert to a matrix. For example, pass
  `which="std.error"` to return standard errors instead of estimates.
  Partial matching supported.

- ...:

  Additional arguments (ignored)

## Value

A data frame with estimates. It has class `ei_est`, supporting several
methods, and two additional attributes: `vcov`, containing the estimated
covariance matrix for the estimates, and `n`, containing the number of
aggregate units used in estimation (the number of rows in `data`).

## Methods (by generic)

- `as.matrix(ei_est)`: Format estimates, standard errors, or other
  columns as a matrix.

- `vcov(ei_est)`: Extract full covariance matrix of estimates

- `nobs(ei_est)`: Extract number of units covered by estimates

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
rr = ei_riesz(spec, penalty = m$penalty)

ei_est(regr = m, data = spec, conf_level = 0.95) # Plug-in estimate
#> # A tibble: 12 × 6
#>    predictor outcome      estimate std.error  conf.low conf.high
#>    <chr>     <chr>           <dbl>     <dbl>     <dbl>     <dbl>
#>  1 vap_white pres_dem_hum  0.230    0.0222    0.186      0.273  
#>  2 vap_black pres_dem_hum  0.588    0.0457    0.498      0.678  
#>  3 vap_other pres_dem_hum  1.55     0.181     1.20       1.91   
#>  4 vap_white pres_rep_nix  0.382    0.0284    0.326      0.438  
#>  5 vap_black pres_rep_nix -0.0775   0.0165   -0.110     -0.0452 
#>  6 vap_other pres_rep_nix  0.494    0.0722    0.353      0.636  
#>  7 vap_white pres_ind_wal  0.387    0.0251    0.338      0.436  
#>  8 vap_black pres_ind_wal  0.488    0.0446    0.401      0.576  
#>  9 vap_other pres_ind_wal -1.05     0.131    -1.30      -0.789  
#> 10 vap_white pres_abs      0.00152  0.000251  0.00102    0.00201
#> 11 vap_black pres_abs      0.00109  0.000227  0.000641   0.00153
#> 12 vap_other pres_abs      0.00118  0.000470  0.000260   0.00210
ei_est(riesz = rr, data = spec) # Weighted (Riesz) estimate
#> # A tibble: 12 × 4
#>    predictor outcome      estimate std.error
#>    <chr>     <chr>           <dbl>     <dbl>
#>  1 vap_white pres_dem_hum  0.230    0.0224  
#>  2 vap_black pres_dem_hum  0.588    0.0948  
#>  3 vap_other pres_dem_hum  1.55     1.98    
#>  4 vap_white pres_rep_nix  0.382    0.0237  
#>  5 vap_black pres_rep_nix -0.0775   0.0695  
#>  6 vap_other pres_rep_nix  0.494    1.49    
#>  7 vap_white pres_ind_wal  0.387    0.0286  
#>  8 vap_black pres_ind_wal  0.488    0.115   
#>  9 vap_other pres_ind_wal -1.05     1.98    
#> 10 vap_white pres_abs      0.00152  0.000252
#> 11 vap_black pres_abs      0.00109  0.000840
#> 12 vap_other pres_abs      0.00118  0.0134  
est = ei_est(regr = m, riesz = rr, data = spec) # Double/debiased ML estimate
# Working with the output
as.matrix(est)
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal     pres_abs
#>   vap_white    0.2296810   0.38141344    0.3873821 0.0015234356
#>   vap_black    0.5856274  -0.07680138    0.4901047 0.0010692501
#>   vap_other    1.6369943   0.59796879   -1.2351021 0.0001390607
as.matrix(est, "std.error")
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal     pres_abs
#>   vap_white   0.02346882   0.03061363   0.02459607 0.0002801553
#>   vap_black   0.04997943   0.02670867   0.05056767 0.0005101201
#>   vap_other   0.47023983   0.33462954   0.55488728 0.0058597715
vcov(est)[1:4, 1:4]
#>                        vap_white:pres_dem_hum vap_black:pres_dem_hum
#> vap_white:pres_dem_hum           0.0005507853           0.0004144131
#> vap_black:pres_dem_hum           0.0004144131           0.0024979436
#> vap_other:pres_dem_hum           0.0021911962           0.0056012336
#> vap_white:pres_rep_nix           0.0006103737           0.0007945651
#>                        vap_other:pres_dem_hum vap_white:pres_rep_nix
#> vap_white:pres_dem_hum            0.002191196           0.0006103737
#> vap_black:pres_dem_hum            0.005601234           0.0007945651
#> vap_other:pres_dem_hum            0.221125499           0.0045939233
#> vap_white:pres_rep_nix            0.004593923           0.0009371941

# Contrasts
ei_est(regr = m, riesz = rr, data = spec, contrast = list(predictor = c(1, -1, 0)))
#> # A tibble: 4 × 4
#>   predictor             outcome       estimate std.error
#>   <chr>                 <chr>            <dbl>     <dbl>
#> 1 vap_white - vap_black pres_dem_hum -0.356     0.0471  
#> 2 vap_white - vap_black pres_rep_nix  0.458     0.0495  
#> 3 vap_white - vap_black pres_ind_wal -0.103     0.0458  
#> 4 vap_white - vap_black pres_abs      0.000454  0.000542
ei_est(regr = m, riesz = rr, data = spec,
       contrast = list(predictor = c(-1, 1, 0), outcome = c(1, -1, 0, 0)))
#> # A tibble: 1 × 4
#>   predictor             outcome                     estimate std.error
#>   <chr>                 <chr>                          <dbl>     <dbl>
#> 1 vap_black - vap_white pres_dem_hum - pres_rep_nix    0.814    0.0700

# Subsetting
est = ei_est(m, rr, data = spec, subset = (state == "Alabama"))
as.matrix(est)
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal    pres_abs
#>   vap_white  0.006649819   0.16292432    0.8158650 0.014560833
#>   vap_black  0.733211405  -0.04036192    0.2995285 0.007622046
#>   vap_other  1.803096491   0.16471413   -0.9702743 0.002463645
nobs(est)
#> [1] 67
```
