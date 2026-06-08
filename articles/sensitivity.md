# Sensitivity analysis for ecological inference

This vignette demonstrates the sensitivity analysis tools in **seine**,
using the `elec_1968` data on county-level voting in Southern states in
the 1968 U.S. presidential election. Sensitivity analysis is essential
for ecological inference (EI) because all EI methods rely on an
untestable identifying assumption—here, Conditional Average
Representativeness, or CAR—that is unlikely to hold exactly in practice.
The tools in **seine** are based on a nonparametric sensitivity
framework developed by Chernozhukov et al. (2024).

## Setting up the analysis

We begin by loading the 1968 election data, and defining an `ei_spec`
object that records the outcome, predictor, covariate, and total-count
columns, following the setup from
[`vignette("seine")`](https://corymccartan.com/seine/articles/seine.md).
We use a BART basis expansion for nonparametric covariate adjustment,
which we strongly recommend to avoid dependence on linearity
assumptions.

``` r

library(seine)
data(elec_1968)

spec = ei_spec(
    elec_1968,
    predictors = vap_white:vap_other,
    outcome = pres_dem_hum:pres_abs,
    total = pres_total,
    covariates = c(state, pop_city:pop_rural, farm:educ_coll, inc_00_03k:inc_25_99k),
    preproc = function(x) {
        x = model.matrix(~ 0 + ., x) # convert factors to dummies
        bases::b_bart(x, trees = 250)
    }
)
```

We fit the regression model with
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) and
the Riesz representer with
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md),
then combine them with
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md) to
estimate vote choice by race using double machine learning (DML). We
focus on the *contrast* between White and Black voters, which is a
direct measure of racially polarized voting. See the main vignette
([`vignette("seine")`](https://corymccartan.com/seine/articles/seine.md))
for a full walkthrough of this estimation workflow.

``` r

m = ei_ridge(spec)
rr = ei_riesz(spec, penalty = m$penalty)

est = ei_est(m, rr, spec, contrast = list(predictor = c(1, -1, 0)), conf_level = FALSE)
print(est)
#> # A tibble: 4 × 4
#>   predictor             outcome       estimate std.error
#>   <chr>                 <chr>            <dbl>     <dbl>
#> 1 vap_white - vap_black pres_dem_hum -0.350     0.0484  
#> 2 vap_white - vap_black pres_rep_nix  0.489     0.0485  
#> 3 vap_white - vap_black pres_ind_wal -0.138     0.0413  
#> 4 vap_white - vap_black pres_abs     -0.000327  0.000975
```

While we normally do not recommend setting `conf_level = FALSE` to
suppress confidence intervals, here we do, so that the output can more
easily fit on the screen. If confidence intervals are present in `est`,
they will be adjusted by the sensitivity analysis below.

## Sensitivity analysis

The estimates above rest on the CAR assumption: that, conditional on the
observed covariates, individual vote choice is independent of the
individual’s race. In practice, this assumption is unlikely to hold
exactly, as there may be unobserved confounders. **seine** provides a
number of tools to evaluate how sensitive the results are to violations
of that assumption.

The sensitivity framework considers the relationship between an
unobserved confounding variable and (i) the outcome and (ii) the Riesz
representer, measured in terms of partial \\R^2\\ values (`c_outcome`
and `c_predictor`, respectively). Stronger relationships indicate more
confounding and therefore more potential bias in the original estimates.

The [`ei_sens()`](https://corymccartan.com/seine/reference/ei_sens.md)
function provides a simple interface to this framework. Users provide
values for the sensitivity parameters, and a bound on the absolute bias
is returned. In the following example, we investigate the effect of an
omitted confounder that explains 50% of the residual variation in the
outcome and 20% of the variation in the Riesz representer.

``` r

ei_sens(est, c_outcome = 0.5, c_predictor = 0.2)
#> # A tibble: 4 × 7
#>   predictor          outcome estimate std.error c_outcome c_predictor bias_bound
#>   <chr>              <chr>      <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
#> 1 vap_white - vap_b… pres_d… -3.50e-1  0.0484         0.5         0.2     0.322 
#> 2 vap_white - vap_b… pres_r…  4.89e-1  0.0485         0.5         0.2     0.381 
#> 3 vap_white - vap_b… pres_i… -1.38e-1  0.0413         0.5         0.2     0.418 
#> 4 vap_white - vap_b… pres_a… -3.27e-4  0.000975       0.5         0.2     0.0143
```

We can also work backwards and ask what one of the sensitivity
parameters would have to be in order to produce a certain amount of
bias. For example, if we assumed a worst-case scenario where the
confounder explains the entire outcome (`c_outcome = 1`), we can ask how
strongly that confounder would need to be related to the Riesz
representer to produce a bias of up to 5pp.

``` r

ei_sens(est, c_outcome = 1, bias_bound = 0.05)
#> # A tibble: 4 × 7
#>   predictor          outcome estimate std.error c_outcome c_predictor bias_bound
#>   <chr>              <chr>      <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
#> 1 vap_white - vap_b… pres_d… -3.50e-1  0.0484           1     0.00300       0.05
#> 2 vap_white - vap_b… pres_r…  4.89e-1  0.0485           1     0.00215       0.05
#> 3 vap_white - vap_b… pres_i… -1.38e-1  0.0413           1     0.00179       0.05
#> 4 vap_white - vap_b… pres_a… -3.27e-4  0.000975         1     0.603         0.05
```

For all of the outcomes except `pres_abs`, whose estimate is much
smaller than 0.05, the answer is not very much!

### Benchmarking

The `c_outcome` parameter is relatively easy to understand, but
`c_predictor` is more difficult to interpret (though see the methodology
paper for more discussion). To help understand plausible values of these
parameters, we can conduct a **benchmarking analysis** that treats each
of our *observed* covariates in turn as a hypothetical *unobserved*
confounder, and calculates the implied values of the sensitivity
parameters.

``` r

bench = ei_bench(spec, contrast = list(predictor = c(1, -1, 0)))
#> ⠙ ETA:22s  Benchmarking state [1/13]
#> ⠹ ETA:19s  Benchmarking pop_urban [3/13]
#> ⠸ ETA:17s  Benchmarking pop_rural [4/13]
#> ⠼ ETA:13s  Benchmarking nonfarm [6/13]
#> ⠴ ETA:11s  Benchmarking educ_elem [7/13]
#> ⠦ ETA: 8s  Benchmarking educ_coll [9/13]
#> ⠧ ETA: 6s  Benchmarking inc_00_03k [10/13]
#> ⠇ ETA: 2s  Benchmarking inc_08_25k [12/13]
#> ⠇ ETA: 0s  Benchmarking inc_25_99k [13/13]

subset(bench, outcome == "pres_rep_nix")
#> # A tibble: 13 × 7
#>    covariate  predictor       outcome c_outcome c_predictor confounding  est_chg
#>    <chr>      <chr>           <chr>       <dbl>       <dbl>       <dbl>    <dbl>
#>  1 state      vap_white - va… pres_r…   0.197        0.277      -0.0737 -2.95e-2
#>  2 pop_city   vap_white - va… pres_r…   0.0129       0.633      -0.0204 -2.81e-3
#>  3 pop_urban  vap_white - va… pres_r…   0.00584      0           1       9.20e-3
#>  4 pop_rural  vap_white - va… pres_r…   0            0.111       1       1.25e-2
#>  5 farm       vap_white - va… pres_r…   0.0236       0.116       0.209   2.01e-2
#>  6 nonfarm    vap_white - va… pres_r…   0.0334       0.305       0.108   1.85e-2
#>  7 educ_elem  vap_white - va… pres_r…   0.00979      0.0473     -0.0144 -5.87e-4
#>  8 educ_hsch  vap_white - va… pres_r…   0.0323       0.154      -0.0897 -1.14e-2
#>  9 educ_coll  vap_white - va… pres_r…   0.0370       0.147       0.271   3.62e-2
#> 10 inc_00_03k vap_white - va… pres_r…   0.0135       0.150       0.336   2.74e-2
#> 11 inc_03_08k vap_white - va… pres_r…   0            0           1       1.49e-2
#> 12 inc_08_25k vap_white - va… pres_r…   0.0239       0.0684      0.0843  6.39e-3
#> 13 inc_25_99k vap_white - va… pres_r…   0.00445      0.0174      1       1.82e-2
```

The table above shows the benchmark values for each covariate for the
racially polarized Nixon vote estimand. The `confounding` column is an
additional component of the sensitivity analysis that is discussed in
the paper; the default value is 1, which is a conservative worst-case
bound. The benchmark values here show that `state` is far and away the
strongest observed confounder, whose inclusion changes the estimate by
29pp. If the unobserved confounders were as strong as `state`, we might
expect a significant amount of bias, as we will see next.

### Bias contour plot

Rather than perform this sensitivity analysis on a single set of
sensitivity parameters, we can run it across all combinations of
parameter values, and visualize the results on a **bias contour plot.**
We can further overlay the benchmarking values to help interpret the
results.

``` r

sens = ei_sens(est) # the default evaluates on a grid of parameters
plot(sens, "pres_rep_nix", bench = bench, bounds = c(-1, 1))
```

![Bias contour plot for the racially polarized Nixon
vote](sensitivity_files/figure-html/unnamed-chunk-7-1.png)

    #> Warning in graphics::par(oldpar): graphical parameter "cin" cannot be set
    #> Warning in graphics::par(oldpar): graphical parameter "cra" cannot be set
    #> Warning in graphics::par(oldpar): graphical parameter "csi" cannot be set
    #> Warning in graphics::par(oldpar): graphical parameter "cxy" cannot be set
    #> Warning in graphics::par(oldpar): graphical parameter "din" cannot be set
    #> Warning in graphics::par(oldpar): graphical parameter "page" cannot be set

The contour lines indicate how much bias could result from an unobserved
confounder with the specified sensitivity parameters. The blue dashed
contours correspond to bias of 1, 2, and 3 standard errors. This is a
helpful value to compare against, because bias of that size corresponds
to a predictable drop in coverage rates of confidence intervals. For
example, bias of 1 standard error means that a confidence interval with
95% nominal coverage will actually have coverage of only around 80%.

The red asterisks indicate the benchmark values for each covariate. Most
are clustered in the lower-left corner and can’t be distinguished. In
contrast, the benchmark for `state` shows that an unobserved confounder
of that strength could lead to bias of around 30pp, which is substantial
compared to the estimate itself, which is 49pp.

### Robustness value

Finally, it can be helpful to summarize the sensitivity analysis by a
single number. The
[`ei_sens_rv()`](https://corymccartan.com/seine/reference/ei_sens_rv.md)
function calculates the **robustness value**, which measures the minimum
strength of an unobserved confounder that would lead to a bias of a
given amount. Here, we consider bias sufficient to eliminate any
evidence of racially polarized voting, i.e., bias equal to the estimated
difference between White and Black voters.

``` r

ei_sens_rv(est, bias_bound = estimate)
#> # A tibble: 4 × 5
#>   predictor             outcome       estimate std.error      rv
#>   <chr>                 <chr>            <dbl>     <dbl>   <dbl>
#> 1 vap_white - vap_black pres_dem_hum -0.350     0.0484   0.317  
#> 2 vap_white - vap_black pres_rep_nix  0.489     0.0485   0.362  
#> 3 vap_white - vap_black pres_ind_wal -0.138     0.0413   0.110  
#> 4 vap_white - vap_black pres_abs     -0.000327  0.000975 0.00803
```

The robustness value (one for each predictor/outcome combination) is
relatively small for Wallace’s vote share, indicating low robustness
(high sensitivity). In particular, it is far smaller than the amount of
confounding benchmarked by the observed `state` variable. For Humphrey
and Nixon’s vote shares, however, the robustness values are larger,
indicating more confidence in the finding of racially polarized voting
for those candidates.

As with any single-number summary, it is important to consider
sensitivity beyond the single value, by using the contour plot and the
benchmarking analysis.

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V.
(2024). Long story short: Omitted variable bias in causal machine
learning (No. w30302). *National Bureau of Economic Research.*
