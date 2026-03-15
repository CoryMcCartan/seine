# Local estimates for ecological inference

Global ecological inference estimates describe average behavior across
all units in the data. Sometimes, however, estimates for each individual
unit (precinct, county, or other geography) are of interest: to map
spatial variation, to examine outliers, or to inform decisions at the
local level. **seine** provides two approaches for this.
[`ei_bounds()`](https://corymccartan.com/seine/reference/ei_bounds.md)
computes guaranteed-valid partial identification bounds for each unit,
relying only on the accounting identity and known limits on the outcome.
[`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md)
produces point estimates and confidence intervals under the Conditional
Average Representativeness (CAR) assumption and a model for within-unit
variation, using a fitted regression model. The former bounds are also
sometimes called the Duncan-Davis bounds.

## Setting up the analysis

We will use the `elec_1968` data included in the package, which contain
county-level election returns from Southern states in the 1968 U.S.
presidential election. We want to estimate the individual-level
association between race and presidential vote choice for each county.

``` r
library(seine)
data(elec_1968)
```

For the bounds, only the data are needed. For
[`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md),
we also need a fitted regression model, so we set up an `ei_spec` object
and call
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md).

``` r
spec = ei_spec(
    elec_1968,
    predictors = vap_white:vap_other,
    outcome = pres_dem_hum:pres_abs,
    total = pres_total,
    covariates = c(state, pop_city:pop_rural, farm:educ_coll, inc_00_03k:inc_25_99k),
    preproc = function(x) {
        x = model.matrix(~ 0 + ., x) # convert factors to dummies
        bases::b_bart(x, trees = 200)
    }
)

m = ei_ridge(spec)
```

See the main vignette
([`vignette("seine")`](https://corymccartan.com/seine/articles/seine.md))
for a full walkthrough of the estimation workflow.

## Partial identification bounds

The simplest local quantities to compute are partial identification
bounds: for each county and each predictor-outcome combination, the
range of values for the local estimand \\\beta\_{gj} = \Pr(Y = j \mid X
= x_g, G = g)\\ that is consistent with the observed aggregate data.
These bounds require no modeling assumptions beyond the accounting
identity and the specified bounds on the outcome.

[`ei_bounds()`](https://corymccartan.com/seine/reference/ei_bounds.md)
takes an `ei_spec` object (or a formula) and returns a data frame with
one row per unit per predictor-outcome combination.

``` r
bnds = ei_bounds(spec, bounds = c(0, 1))
head(bnds)
#> # A tibble: 6 × 6
#>    .row predictor outcome      weight    min    max
#>   <int> <chr>     <chr>         <dbl>  <dbl>  <dbl>
#> 1     1 vap_white pres_dem_hum  5877. 0      0.262 
#> 2     2 vap_white pres_dem_hum 16131. 0      0.122 
#> 3     3 vap_white pres_dem_hum  4872. 0      0.397 
#> 4     4 vap_white pres_dem_hum  3566. 0      0.180 
#> 5     5 vap_white pres_dem_hum  8801. 0.0190 0.0382
#> 6     6 vap_white pres_dem_hum  1698. 0      1
```

The `min` and `max` columns give the sharp bounds for each unit. The
`weight` column is the product of the predictor proportion and the
total, which is useful for aggregating. When the bounds are wide for
many units, this indicates a fundamentally difficult identification
problem, regardless of modeling choices.

To aggregate the bounds across all units to produce bounds on the global
estimands, pass `global = TRUE`.

``` r
ei_bounds(spec, bounds = c(0, 1), global = TRUE)
#> # A tibble: 12 × 4
#>    predictor outcome            min     max
#>    <chr>     <chr>            <dbl>   <dbl>
#>  1 vap_white pres_dem_hum 0.176     0.375  
#>  2 vap_black pres_dem_hum 0.00467   0.927  
#>  3 vap_other pres_dem_hum 0         1      
#>  4 vap_white pres_rep_nix 0.247     0.421  
#>  5 vap_black pres_rep_nix 0         0.808  
#>  6 vap_other pres_rep_nix 0         0.991  
#>  7 vap_white pres_ind_wal 0.213     0.415  
#>  8 vap_black pres_ind_wal 0.00874   0.946  
#>  9 vap_other pres_ind_wal 0         1.000  
#> 10 vap_white pres_abs     0.0000108 0.00179
#> 11 vap_black pres_abs     0         0.00844
#> 12 vap_other pres_abs     0         0.0910
```

The [`as.array()`](https://rdrr.io/r/base/array.html) method reformats
the output as a three-dimensional array with dimensions units by
predictors by outcomes, which can be convenient for further analysis.

``` r
bnds_arr = as.array(bnds)
dim(bnds_arr) # rows x predictors x outcomes
#> [1] 1143    3    4    2
```

Bounds on contrasts, such as the difference in vote share between White
and Black voters, can be computed directly by passing a `contrast`
argument.

``` r
ei_bounds(spec, bounds = c(0, 1), contrast = list(predictor = c(1, -1, 0)))
#> # A tibble: 4,572 × 5
#>     .row predictor             outcome         min    max
#>    <int> <chr>                 <chr>         <dbl>  <dbl>
#>  1     1 vap_white - vap_black pres_dem_hum -0.841 0.262 
#>  2     2 vap_white - vap_black pres_dem_hum -0.761 0.122 
#>  3     3 vap_white - vap_black pres_dem_hum -0.623 0.397 
#>  4     4 vap_white - vap_black pres_dem_hum -0.654 0.180 
#>  5     5 vap_white - vap_black pres_dem_hum -0.981 0.0382
#>  6     6 vap_white - vap_black pres_dem_hum -0.742 0.901 
#>  7     7 vap_white - vap_black pres_dem_hum -0.531 0.250 
#>  8     8 vap_white - vap_black pres_dem_hum -0.993 0.183 
#>  9     9 vap_white - vap_black pres_dem_hum -0.474 0.178 
#> 10    10 vap_white - vap_black pres_dem_hum -0.991 0.0855
#> # ℹ 4,562 more rows
```

## Estimating the local covariance

Point estimates and confidence intervals from
[`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md)
require a model for how the local estimands vary around their
conditional mean within each unit. This variation is captured by the
`b_cov` argument, a covariance matrix for the local estimands
\\\beta\_{gj}\\ around the regression prediction.

[`ei_local_cov()`](https://corymccartan.com/seine/reference/ei_local_cov.md)
estimates this covariance from the data, under the CAR assumption and an
additional homoskedasticity condition. Concretely, it fits a ridge
regression of the empirical second moments of the regression residuals
on the cross-products of the predictor variables, and uses the
polarization identity to recover the covariance for each
predictor-outcome pair. A small amount of shrinkage towards the global
residual covariance is applied, controlled by `prior_obs`.

``` r
b_cov = ei_local_cov(m, spec)
round(b_cov, 3)
#>                        vap_white:pres_dem_hum vap_black:pres_dem_hum
#> vap_white:pres_dem_hum                  0.006                 -0.008
#> vap_black:pres_dem_hum                 -0.008                  0.070
#> vap_other:pres_dem_hum                 -0.005                  0.058
#> vap_white:pres_rep_nix                 -0.003                 -0.002
#> vap_black:pres_rep_nix                  0.007                  0.006
#> vap_other:pres_rep_nix                  0.004                  0.010
#> vap_white:pres_ind_wal                 -0.003                  0.010
#> vap_black:pres_ind_wal                  0.001                 -0.076
#> vap_other:pres_ind_wal                  0.001                 -0.069
#> vap_white:pres_abs                      0.000                  0.000
#> vap_black:pres_abs                      0.000                  0.000
#> vap_other:pres_abs                      0.000                  0.000
#>                        vap_other:pres_dem_hum vap_white:pres_rep_nix
#> vap_white:pres_dem_hum                 -0.005                 -0.003
#> vap_black:pres_dem_hum                  0.058                 -0.002
#> vap_other:pres_dem_hum                  0.060                  0.014
#> vap_white:pres_rep_nix                  0.014                  0.031
#> vap_black:pres_rep_nix                 -0.001                 -0.022
#> vap_other:pres_rep_nix                 -0.020                 -0.056
#> vap_white:pres_ind_wal                 -0.009                 -0.027
#> vap_black:pres_ind_wal                 -0.058                  0.024
#> vap_other:pres_ind_wal                 -0.040                  0.042
#> vap_white:pres_abs                      0.000                  0.000
#> vap_black:pres_abs                      0.001                  0.000
#> vap_other:pres_abs                      0.000                  0.000
#>                        vap_black:pres_rep_nix vap_other:pres_rep_nix
#> vap_white:pres_dem_hum                  0.007                  0.004
#> vap_black:pres_dem_hum                  0.006                  0.010
#> vap_other:pres_dem_hum                 -0.001                 -0.020
#> vap_white:pres_rep_nix                 -0.022                 -0.056
#> vap_black:pres_rep_nix                  0.025                  0.041
#> vap_other:pres_rep_nix                  0.041                  0.104
#> vap_white:pres_ind_wal                  0.015                  0.052
#> vap_black:pres_ind_wal                 -0.031                 -0.052
#> vap_other:pres_ind_wal                 -0.040                 -0.084
#> vap_white:pres_abs                      0.000                  0.000
#> vap_black:pres_abs                      0.000                  0.000
#> vap_other:pres_abs                      0.000                  0.000
#>                        vap_white:pres_ind_wal vap_black:pres_ind_wal
#> vap_white:pres_dem_hum                 -0.003                  0.001
#> vap_black:pres_dem_hum                  0.010                 -0.076
#> vap_other:pres_dem_hum                 -0.009                 -0.058
#> vap_white:pres_rep_nix                 -0.027                  0.024
#> vap_black:pres_rep_nix                  0.015                 -0.031
#> vap_other:pres_rep_nix                  0.052                 -0.052
#> vap_white:pres_ind_wal                  0.030                 -0.025
#> vap_black:pres_ind_wal                 -0.025                  0.108
#> vap_other:pres_ind_wal                 -0.043                  0.110
#> vap_white:pres_abs                      0.000                  0.000
#> vap_black:pres_abs                      0.000                 -0.001
#> vap_other:pres_abs                      0.000                  0.000
#>                        vap_other:pres_ind_wal vap_white:pres_abs
#> vap_white:pres_dem_hum                  0.001                  0
#> vap_black:pres_dem_hum                 -0.069                  0
#> vap_other:pres_dem_hum                 -0.040                  0
#> vap_white:pres_rep_nix                  0.042                  0
#> vap_black:pres_rep_nix                 -0.040                  0
#> vap_other:pres_rep_nix                 -0.084                  0
#> vap_white:pres_ind_wal                 -0.043                  0
#> vap_black:pres_ind_wal                  0.110                  0
#> vap_other:pres_ind_wal                  0.125                  0
#> vap_white:pres_abs                      0.000                  0
#> vap_black:pres_abs                     -0.001                  0
#> vap_other:pres_abs                      0.000                  0
#>                        vap_black:pres_abs vap_other:pres_abs
#> vap_white:pres_dem_hum              0.000                  0
#> vap_black:pres_dem_hum              0.000                  0
#> vap_other:pres_dem_hum              0.001                  0
#> vap_white:pres_rep_nix              0.000                  0
#> vap_black:pres_rep_nix              0.000                  0
#> vap_other:pres_rep_nix              0.000                  0
#> vap_white:pres_ind_wal              0.000                  0
#> vap_black:pres_ind_wal             -0.001                  0
#> vap_other:pres_ind_wal             -0.001                  0
#> vap_white:pres_abs                  0.000                  0
#> vap_black:pres_abs                  0.000                  0
#> vap_other:pres_abs                  0.000                  0
```

The rows and columns are ordered by predictor within outcome, i.e.
(Y1\|X1, Y1\|X2, …, Y2\|X1, Y2\|X2, …). Large diagonal entries indicate
high within-unit variation in the local estimands, relative to the
regression predictions. Off-diagonal entries reflect the extent to
which, within a county, vote choices across racial groups tend to move
together.

## Local point estimates

With the regression model and a covariance structure in hand,
[`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md)
projects the regression predictions onto the accounting constraint to
produce local estimates that exactly satisfy the ecological identity
within each county. The function also computes valid confidence
intervals using Chebyshev’s inequality or, under a unimodality
assumption, the tighter unimodal bound.

The `b_cov` argument controls the assumed covariance structure. There
are three common choices:

- `b_cov = ei_local_cov(m, spec)`: the **data-estimated covariance**,
  which uses the structure estimated above. *This is the recommended
  choice*.
- `b_cov = 0`: the **orthogonal model**, which assumes no correlation in
  the local estimands across predictors within a unit. When the
  estimated `b_cov` is implausible, this can be an acceptable choice.
- `b_cov = 0.95` (or another value near 1): the **neighborhood model**,
  which assumes nearly perfect positive correlation across predictors.
  This is appropriate when counties are relatively homogeneous and the
  local estimands are expected to vary together. However, it leads to
  quite narrow confidence intervals, which may not cover the truth in
  practice.

We fit all three and compare the average width of the resulting
confidence intervals.

``` r
e_rcov = ei_est_local(m, spec, b_cov = b_cov, bounds = c(0, 1), sum_one = TRUE)
e_orth = ei_est_local(m, spec, b_cov = 0,    bounds = c(0, 1), sum_one = TRUE)
e_nbhd = ei_est_local(m, spec, b_cov = 0.95, bounds = c(0, 1), sum_one = TRUE)

c(
    estimated = mean(e_rcov$conf.high - e_rcov$conf.low),
    orthogonal = mean(e_orth$conf.high - e_orth$conf.low),
    neighborhood = mean(e_nbhd$conf.high - e_nbhd$conf.low)
)
#>    estimated   orthogonal neighborhood 
#>    0.4165165    0.2876576    0.2334284
```

The orthogonal model tends to produce the widest intervals. The
neighborhood model can be narrower when the predictor proportions within
a county are nearly uniform, because the accounting constraint then
binds tightly. The estimated covariance typically lies between the two
extremes, and as mentioned above is the recommended default.

Setting `sum_one = TRUE` enforces the constraint that the local vote
shares across candidates sum to one within each racial group, which is
the correct restriction for these data. The `bounds = c(0, 1)` argument
truncates estimates to the unit interval and caps the standard errors
implied by the width of the bounds.

To examine results for a specific county, we can filter the output.

``` r
subset(e_rcov, .row == 1)
#> # A tibble: 12 × 8
#>     .row predictor outcome      weight estimate std.error conf.low conf.high
#>    <int> <chr>     <chr>         <dbl>    <dbl>     <dbl>    <dbl>     <dbl>
#>  1     1 vap_white pres_dem_hum 5877.   0.111     0.0684   0          0.262 
#>  2     1 vap_black pres_dem_hum 1831.   0.478     0.218    0          0.841 
#>  3     1 vap_other pres_dem_hum   13.3  0.902     0.244    0.176      1     
#>  4     1 vap_white pres_rep_nix 5877.   0.101     0.0293   0.0139     0.102 
#>  5     1 vap_black pres_rep_nix 1831.   0         0.0941   0          0.281 
#>  6     1 vap_other pres_rep_nix   13.3  0.0964    0.233    0          0.791 
#>  7     1 vap_white pres_ind_wal 5877.   0.775     0.0627   0.620      0.934 
#>  8     1 vap_black pres_ind_wal 1831.   0.512     0.200    0          1     
#>  9     1 vap_other pres_ind_wal   13.3  0         0.289    0          0.861 
#> 10     1 vap_white pres_abs     5877.   0.0128    0.00119  0.00927    0.0160
#> 11     1 vap_black pres_abs     1831.   0.0103    0.00385  0          0.0218
#> 12     1 vap_other pres_abs       13.3  0.00123   0.00954  0          0.0297
```

The [`as.array()`](https://rdrr.io/r/base/array.html) method provides a
convenient view of the point estimates as a three-dimensional array.

``` r
head(as.array(e_rcov)[, , "pres_dem_hum"])
#>        vap_white vap_black vap_other
#> [1,] 0.111049420 0.4776159 0.9024011
#> [2,] 0.053312727 0.4074466 0.9000749
#> [3,] 0.003871053 0.6155130 0.8955696
#> [4,] 0.041992105 0.4968093 0.8936855
#> [5,] 0.029622310 0.4278036 0.8956958
#> [6,] 0.108993598 0.6715685 0.8941690
```

Finally, taking a weighted average of the local estimates against the
`weight` column reproduces the global estimate from
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md).

``` r
grps = split(e_rcov, ~ predictor + outcome)
sapply(grps, function(g) weighted.mean(g$estimate, g$weight))
#>     vap_black.pres_abs     vap_other.pres_abs     vap_white.pres_abs 
#>           1.520710e-03           7.274520e-04           1.552457e-03 
#> vap_black.pres_dem_hum vap_other.pres_dem_hum vap_white.pres_dem_hum 
#>           5.166252e-01           9.086030e-01           2.625707e-01 
#> vap_black.pres_ind_wal vap_other.pres_ind_wal vap_white.pres_ind_wal 
#>           4.691351e-01           1.498855e-10           3.182293e-01 
#> vap_black.pres_rep_nix vap_other.pres_rep_nix vap_white.pres_rep_nix 
#>           1.271900e-02           9.066955e-02           4.176476e-01
```

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V.
(2024). Long story short: Omitted variable bias in causal machine
learning (No. w30302). *National Bureau of Economic Research.*

This vignette was originally produced by a large language model, and
then reviewed and edited by the package authors.
