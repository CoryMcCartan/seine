# Conduct a sensitivity analysis for estimated ecological quantities

Relates confounding of an omitted variable with predictor or outcome to
bias in ecological estimates, using the nonparametric sensitivity
analysis of Chernozhukov et al. (2024).

## Usage

``` r
ei_sens(
  est,
  c_outcome = seq(0, 1, 0.01)^2,
  c_predictor = seq(0, 1, 0.01)^2,
  bias_bound = NULL,
  confounding = 1,
  expand_ci = TRUE
)
```

## Arguments

- est:

  A set of estimates from
  [`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md) using
  both regression and Riesz representer.

- c_outcome:

  The (nonparametric) partial \\R^2\\ of the omitted variables with the
  outcome variables. Must be between 0 and 1. Can be a vector, in which
  case all combinations of values with `c_predictor` are used.

- c_predictor:

  How much variation latent variables create in the Riesz representer,
  i.e. \\1-R^2\\ of the true Riesz representer on the estimated one
  without the omitted variable. Must be between 0 and 1. Can be a
  vector, in which case all combinations of values with `c_outcome` are
  used.

- bias_bound:

  If provided, overrides `c_predictor` and finds values of `c_predictor`
  that correspond to (the absolute value of) the provided amount of
  bias.

- confounding:

  The confounding parameter (\\\rho\\), which must be between 0 and 1
  (the adversarial worst-case).

- expand_ci:

  If `TRUE` and confidence intervals are present in `est`, expand the
  width of the intervals in each direction by the calculated bias bound.

## Value

A data frame of the same format as `est`, but with additional columns:
`c_outcome` and `c_predictor`, matching all combinations of those
arguments, and `bias_bound`, containing the bound on the amount of bias.
The data frame has additional class `ei_sens`, which supports a
[`plot.ei_sens()`](https://corymccartan.com/seine/reference/plot.ei_sens.md)
method.

## Details

The parameter `c_predictor` equals \\1 - R^2\_{\alpha\sim\alpha_s}\\,
where \\\alpha\\ is the true Riesz representer and \\\alpha_s\\ is the
Riesz representer with the observed covariates. The RR can be
equivalently expressed as \$\$ \alpha(n, \bar x_j, z, a) = u(n, \bar
x_j)\partial\_{\bar x_j} \log f(\bar x_j\mid z, a), \$\$ where \\u(n,
\bar x_j)\\ is the weighting term, \\A\\ is the unobserved confounder,
and \\f\\ is the conditional density. The corresponding `c_predictor` is
then \$\$ 1 - R^2\_{\alpha\sim\alpha_s} = 1 - \\ \frac{\mathbb{E}\[u(N,
\bar X_j)(\partial\_{\bar x_j} \log f(\bar x_j\mid Z))^2\]}{
\mathbb{E}\[u(N, \bar X_j)(\partial\_{\bar x_j} \log f(\bar x_j\mid Z,
A))^2\]}. \$\$

The bounds here are plug-in estimates and do not incorporate sampling
uncertainty. As such, they may fail to cover the true value in finite
samples, even under large enough sensitivity parameters; see Section 5
of Chernozhukov et al. (2024).

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

ei_sens(est, c_outcome=0.2)
#> # A tibble: 303 × 7
#>    predictor outcome      estimate std.error c_outcome c_predictor bias_bound
#>    <chr>     <chr>           <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
#>  1 vap_white pres_ind_wal    0.387    0.0246       0.2      0         0      
#>  2 vap_black pres_ind_wal    0.490    0.0505       0.2      0         0      
#>  3 vap_other pres_ind_wal   -1.19     0.528        0.2      0         0      
#>  4 vap_white pres_ind_wal    0.387    0.0246       0.2      0.0001    0.00111
#>  5 vap_black pres_ind_wal    0.490    0.0505       0.2      0.0001    0.00392
#>  6 vap_other pres_ind_wal   -1.19     0.528        0.2      0.0001    0.0909 
#>  7 vap_white pres_ind_wal    0.387    0.0246       0.2      0.0004    0.00221
#>  8 vap_black pres_ind_wal    0.490    0.0505       0.2      0.0004    0.00784
#>  9 vap_other pres_ind_wal   -1.19     0.528        0.2      0.0004    0.182  
#> 10 vap_white pres_ind_wal    0.387    0.0246       0.2      0.0009    0.00332
#> # ℹ 293 more rows

# How much variation would the regression residual need to explain of
# Riesz representer to cause bias of 0.4?
ei_sens(est, c_outcome=1, bias_bound=0.4)
#> # A tibble: 3 × 7
#>   predictor outcome      estimate std.error c_outcome c_predictor bias_bound
#>   <chr>     <chr>           <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
#> 1 vap_white pres_ind_wal    0.387    0.0246         1    0.723           0.4
#> 2 vap_black pres_ind_wal    0.490    0.0505         1    0.172           0.4
#> 3 vap_other pres_ind_wal   -1.19     0.528          1    0.000388        0.4

# Update confidence intervals and extract as matrix
est = ei_est(m, rr, spec, conf_level=0.95)
sens = ei_sens(est, c_outcome=0.5, c_predictor=0.2)
as.matrix(sens, "conf.high")
#>            outcome
#> predictor   pres_ind_wal
#>   vap_white    0.5229302
#>   vap_black    0.8985835
#>   vap_other    7.0320672

# Works for contrasts as well
est = ei_est(m, rr, spec, contrast = list(predictor=c(1, -1, 0)))
ei_sens(est, c_outcome=0.5, c_predictor=0.5)
#> # A tibble: 1 × 7
#>   predictor          outcome estimate std.error c_outcome c_predictor bias_bound
#>   <chr>              <chr>      <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
#> 1 vap_white - vap_b… pres_i…   -0.102    0.0457       0.5         0.5      0.737
```
