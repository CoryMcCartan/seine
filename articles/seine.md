# Analyzing aggregate data

Ecological inference (EI) is the statistical problem of learning
individual-level associations from aggregate-level data. EI commonly
arises when two datasets are joined using a shared geographic
identifier, and when individual data are not released for privacy
reasons. To take some recent examples from the *New York Times*:

- [Estimating COVID vaccine uptake by political
  beliefs](https://www.nytimes.com/interactive/2021/04/17/us/vaccine-hesitancy-politics.html)
- [Understanding which demographics supported a progressive mayoral
  candidate](https://www.nytimes.com/interactive/2025/06/25/nyregion/nyc-mayor-election-results-map-mamdani-cuomo.html)
- [Evaluating the differential impact of tariffs on political
  partisans](https://www.nytimes.com/interactive/2025/03/15/business/economy/tariffs-trump-maps-voters.html)

EI is also used in public health and epidemiology, and is widely applied
in litigation under the federal Voting Rights Act of 1965 (VRA) to
establish the presence of racially polarized voting.

## Preparing data

As an example of an ecological analysis, we will use the `elec_1968`
data included in the package. The data contain county-level election
returns from Southern states in the 1968 U.S. presidential election
along with a number of covariates taken from the 1970 U.S. census. The
counties here are the **aggregation units**; in other analyses, states,
precincts, or cities might be the aggregation units.

``` r
library(seine)
data(elec_1968)

print(elec_1968)
#> # A tibble: 1,143 × 41
#>    fips  state  abbr  region division county    pop pop_city pop_urban pop_rural
#>    <chr> <chr>  <chr> <chr>  <chr>    <chr>   <dbl>    <dbl>     <dbl>     <dbl>
#>  1 01001 Alaba… AL    South  East So… Autau…  24460        0     0.536     0.464
#>  2 01003 Alaba… AL    South  East So… Baldw…  59382        0     0.266     0.734
#>  3 01005 Alaba… AL    South  East So… Barbo…  22543        0     0.404     0.596
#>  4 01007 Alaba… AL    South  East So… Bibb …  13812        0     0         1    
#>  5 01009 Alaba… AL    South  East So… Bloun…  26853        0     0.163     0.837
#>  6 01011 Alaba… AL    South  East So… Bullo…  11824        0     0.366     0.634
#>  7 01013 Alaba… AL    South  East So… Butle…  22007        0     0.365     0.635
#>  8 01015 Alaba… AL    South  East So… Calho… 103092        0     0.641     0.359
#>  9 01017 Alaba… AL    South  East So… Chamb…  36356        0     0.437     0.563
#> 10 01019 Alaba… AL    South  East So… Chero…  15606        0     0         1    
#> # ℹ 1,133 more rows
#> # ℹ 31 more variables: pop_white <dbl>, pop_black <dbl>, pop_aian <dbl>,
#> #   pop_asian <dbl>, pop_hisp <dbl>, vap <dbl>, vap_white <dbl>,
#> #   vap_black <dbl>, vap_other <dbl>, farm <dbl>, nonfarm <dbl>,
#> #   educ_elem <dbl>, educ_hsch <dbl>, educ_coll <dbl>, cvap <dbl>,
#> #   cvap_white <dbl>, cvap_black <dbl>, cvap_other <dbl>, inc_00_03k <dbl>,
#> #   inc_03_08k <dbl>, inc_08_25k <dbl>, inc_25_99k <dbl>, pres_dem_hum <dbl>, …
```

We are interested in estimating the individual-level association between
race and presidential vote choice. The **outcome** variables are the
proportion of votes cast for each candidate: `pres_dem_hum`,
`pres_rep_nix`, `pres_ind_wal`, and `pres_abs`, where the latter are
abstentions and ballots cast for other candidates. The **predictor**
variables are the proportions of the voting-age population in each
racial group: `vap_white`, `vap_black`, and `vap_other`. The data also
contain a number of **covariates**, such as education and income, which
we discuss below.

Ideally, these would be the proportion of each racial group within the
population that actually cast a ballot for President. Since those
proportions are unobserved, they would have to be estimated using a
first stage of ecological inference with an outcome variable measuring
turnout. Alternatively, one could include non-voters as another category
of outcome variable, so that both outcome and predictor variables are
proportions relative to the total voting-age population. For
demonstration purposes, we will ignore this issue and proceed as if
turnout were uniform across racial groups in every county.

These data have already been cleaned. Often, outcomes and predictors are
measured as counts, or may have been rounded, so that they do not sum to
exactly 1. **seine** provides the
[`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
function to assist in preprocessing. To see this in action, suppose we
wanted to set up the turnout problem mentioned in the previous
paragraph. The
[`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
function would let us create a new turnout proportion variable from our
existing data.

``` r
elec_1968_turn = ei_proportions(elec_1968, turnout = pres_total,
                                .total = vap, clamp = 0.01)

subset(elec_1968_turn, select = c(fips, state, county, turnout, vap, .other))
#> # A tibble: 1,143 × 6
#>    fips  state   county          turnout   vap .other
#>    <chr> <chr>   <chr>             <dbl> <dbl>  <dbl>
#>  1 01001 Alabama Autauga County    0.606 12744  0.394
#>  2 01003 Alabama Baldwin County    0.568 33012  0.432
#>  3 01005 Alabama Barbour County    0.645 12370  0.355
#>  4 01007 Alabama Bibb County       0.601  7575  0.399
#>  5 01009 Alabama Blount County     0.566 15856  0.434
#>  6 01011 Alabama Bullock County    0.721  6019  0.279
#>  7 01013 Alabama Butler County     0.593 12341  0.407
#>  8 01015 Alabama Calhoun County    0.484 55547  0.516
#>  9 01017 Alabama Chambers County   0.504 21117  0.496
#> 10 01019 Alabama Cherokee County   0.614  9215  0.386
#> # ℹ 1,133 more rows
```

The function normalizes `pres_total` by `vap` and stores the result in a
column labeled `turnout`. It also stores the remaining proportion (i.e.,
the non-voters) in the `.other` column, by default. In this data, there
is one county which had higher turnout than 1970 VAP. The `clamp = 0.01`
argument tells
[`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
to allow that kind of excess up to 1% of the total, and round those
proportions down to 1. Any proportions in excess of 1.01 would throw an
error. You can read about other functionality and customization of
[`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
in the function’s documentation.

## Avoiding the ecological fallacy

The core challenge of ecological inference is that only marginal
proportions are observed (racial groups, candidate vote shares), but we
are interested in joint data (candidate vote shares *within* each racial
group). The key to overcoming this challenge is assuming some kind of
homogeneity across aggregation units. Enough homogeneity means that
information can be shared across aggregation units to estimate the
missing joint proportions.

More precisely, a researcher needs to believe that **coarsening at
random** (CAR) holds in order to conduct EI. Coarsening at random means
that unobserved joint data of interest are mean-independent of the
predictors and the number of people in each aggregation unit, given
covariates.[¹](#fn1)

In these data, CAR means that once we know a set of covariate values for
a county, such as its education and age, learning about the racial
composition of the county does not change our beliefs about the
candidate preference within each racial group or the total turnout.

For example, take the three counties shown below, which have been
selected by a clustering algorithm to be similar on the observed
covariates: urbanity, agriculture, education, and income.

| state    | county              | pop_urban |   farm | educ_elem | educ_hsch | educ_coll | inc_00_03k | inc_03_08k | inc_08_25k | inc_25_99k |
|:---------|:--------------------|----------:|-------:|----------:|----------:|----------:|-----------:|-----------:|-----------:|-----------:|
| Virginia | Charles City County |         0 | 0.0409 |    0.5401 |    0.3669 |    0.0930 |     0.1778 |     0.4486 |     0.3570 |     0.0166 |
| Virginia | Greene County       |         0 | 0.0737 |    0.5827 |    0.3517 |    0.0656 |     0.1736 |     0.3832 |     0.4368 |     0.0064 |
| Virginia | Louisa County       |         0 | 0.0620 |    0.5300 |    0.3964 |    0.0736 |     0.1823 |     0.4200 |     0.3814 |     0.0163 |

CAR means that the preference for, e.g., George Wallace among White
voters in these counties is roughly the same and is unrelated to the
fact that the demographics are quite different between the counties:

| state    | county              | vap_white | vap_black | vap_other | pres_total |
|:---------|:--------------------|----------:|----------:|----------:|-----------:|
| Virginia | Charles City County |    0.2138 |    0.7000 |    0.0862 |       1960 |
| Virginia | Greene County       |    0.9073 |    0.0917 |    0.0010 |       1549 |
| Virginia | Louisa County       |    0.6669 |    0.3318 |    0.0014 |       3964 |

If we believed that in the majority-Black Charles City County, racial
resentment might increase the preference for Wallace compared to the
heavily majority-White Greene County, then CAR would be violated.

For now, we will proceed under the CAR assumption, though there are
serious reasons to doubt its applicability in these data. Later, we’ll
discuss how to conduct a [sensitivity analysis](#sensitivity-analysis)
to evaluate how possible violations might affect our conclusions.

## Ecological estimation

Once we’ve evaluated the CAR assumption, we can proceed with estimation.
**seine** implements double/debiased machine learning (DML), which means
we fit two models before combining them for a final estimate:

1.  A **regression** model of the outcome variables on the predictor
    variables and covariates
2.  A **Riesz representer** model, which yields a special set of
    “weights” that can be used in estimation.

By carefully combining the fitted regression and Riesz representer, we
can reduce the sensitivity to biases in each component.

### Setup

**seine** provides both a formula interface and a tidy interface through
a new [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
object. We recommend the
[`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
approach for most analyses, since it dovetails well with the other
estimation and sensitivity functions. We will demonstrate both
approaches here, however.

To create an EI *specification*, we call
[`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md) and
use `tidyselect` syntax to specify the outcome, predictors, covariates,
and the column with the total number of people in each aggregation unit.
The function returns an `ei_spec` object, which is just a data frame
with some additional metadata about these variables.

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

print(spec)
#> EI Specification
#> • Predictors: `vap_white`, `vap_black`, and `vap_other`
#> • Outcome: `pres_dem_hum`, `pres_rep_nix`, `pres_ind_wal`, and `pres_abs`
#> • Covariates (465 after preprocessing):`state`, `pop_city`, `pop_urban`, `pop_rural`, `farm`, `nonfarm`, `educ_elem`, `educ_hsch`, `educ_coll`, `inc_00_03k`, `inc_03_08k`, `inc_08_25k`, and `inc_25_99k`
#> # A tibble: 1,143 × 20
#>   vap_white vap_black vap_other pres_dem_hum pres_rep_nix pres_ind_wal pres_abs
#>       <dbl>     <dbl>     <dbl>        <dbl>        <dbl>        <dbl>    <dbl>
#> 1     0.761    0.237   0.00173        0.199        0.0773        0.711  0.0122 
#> 2     0.860    0.137   0.00306        0.105        0.115         0.764  0.0161 
#> 3     0.610    0.389   0.000808       0.242        0.0489        0.687  0.0218 
#> 4     0.783    0.216   0.00106        0.141        0.0571        0.799  0.00290
#> 5     0.981    0.0181  0.000757       0.0375       0.222         0.727  0.0134 
#> # ℹ 1,138 more rows
#> # ℹ 13 more variables: state <chr>, pop_city <dbl>, pop_urban <dbl>,
#> #   pop_rural <dbl>, farm <dbl>, nonfarm <dbl>, educ_elem <dbl>,
#> #   educ_hsch <dbl>, educ_coll <dbl>, inc_00_03k <dbl>, inc_03_08k <dbl>,
#> #   inc_08_25k <dbl>, inc_25_99k <dbl>
```

The only other argument to
[`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md) is
`preproc`, which describes preprocessing done to the covariates before
model fitting. This argument powers the nonparametric estimation in
**seine**: by using various basis expansions in `preproc`, flexible and
assumption-lean models can be fit. We *strongly recommend* using a
nonparametric basis expansion, because otherwise the EI estimates are
dependent on the covariates entering the regression model linearly.

Here, we are using `b_bart()` from the
[`bases`](https://cran.r-project.org/package=bases) package, which
produces a basis expansion that allows for approximately fitting a
Bayesian Additive Regression Trees (BART) model. Other options include
`b_tpsob()`, `b_rff()`, and `b_inter()`, or functions from the `splines`
package.

### Fitting the regression

Any machine learning method can be used to fit the regression model.
However, due to the aggregation process that led to our data, there is
certain structure in the regression function that can be leveraged for
improved estimation. We recommend using
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) to
fit the regression model, because it will automatically use this
structure, and automatically determine the ridge penalty using a
closed-form expression for the leave-one-out errors.

Using the tidy interface, fitting the regression is as simple as calling
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) on
the `ei_spec` object:

``` r
m = ei_ridge(spec)

print(m)
#> An ecological inference model with 4 outcomes, 3 groups, and 1143 observations
#> Fit with penalty = 226.687
```

We can see that
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) has
automatically selected a small ridge penalty. By default, all covariates
are centered and scaled to have unit variance. This is generally
appropriate when penalizing all coefficients equally, as is done by
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md).
But in some cases it may not be appropriate, and this behavior can be
suppressed by providing `scale = FALSE`.

Alternatively, we could use the formula interface, which would also let
us specify our own interaction terms; here, we interact `state` with all
other variables. Nonparametric basis expansions like
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) and
[`bases::b_tpsob()`](http://corymccartan.com/bases/reference/b_tpsob.md)
can also be used in the formula interface. Formulas in **seine** require
the user to separate the predictors and covariates by a vertical bar.

``` r
m_form = ei_ridge(
    cbind(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs) ~
        vap_white + vap_black + vap_other |
        state * (pop_urban + pop_rural + farm + educ_hsch + educ_coll +
                     inc_03_08k + inc_08_25k + inc_25_99k),
    data = elec_1968, total = pres_total
)

print(m_form)
#> An ecological inference model with 4 outcomes, 3 groups, and 1143 observations
#> Fit with penalty = 5.22562
```

The [`summary()`](https://rdrr.io/r/base/summary.html) method of fitted
regression objects shows summary statistics for fitted values, which can
help diagnose misspecification, and shows the \\R^2\\ values for each
outcome variable. Here, racial demographics and covariates explain a
substantial amount of the total variation in vote shares. The fitted
values are almost all between 0 and 1, but the presence of some negative
predictions indicates there is at least some model misspecification.

``` r
summary(m)
#> Fitted values:
#>   pres_dem_hum        pres_rep_nix        pres_ind_wal        pres_abs         
#>  Min.   :-0.005007   Min.   :-0.006611   Min.   :0.03058   Min.   :-0.0022343  
#>  1st Qu.: 0.237510   1st Qu.: 0.199913   1st Qu.:0.26450   1st Qu.:-0.0001501  
#>  Median : 0.291894   Median : 0.308484   Median :0.38012   Median : 0.0001260  
#>  Mean   : 0.301626   Mean   : 0.296998   Mean   :0.40004   Mean   : 0.0013324  
#>  3rd Qu.: 0.382241   3rd Qu.: 0.386513   3rd Qu.:0.53370   3rd Qu.: 0.0006569  
#>  Max.   : 0.692182   Max.   : 0.714572   Max.   :0.82686   Max.   : 0.0294563  
#> 
#> R-squared by outcome:
#> pres_dem_hum pres_rep_nix pres_ind_wal     pres_abs 
#>    0.7210448    0.7432126    0.8156087    0.6104691
```

### Fitting the Riesz representer

The Riesz representer is less familiar, but no less easy to fit. Using
the tidy interface, we simply pass the `ei_spec` object to
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md).
Unlike
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md),
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md)
requires a penalty to be specified. A good default is to use the same
penalty as was used in the regression.

``` r
rr = ei_riesz(spec, penalty = m$penalty)
```

We could also use the formula interface. It is critical to provide
exactly the same formula and data to both
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md) and
[`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md)
(though the Riesz representer does not use the outcome variable); the
tidy interface obviates the need to worry about this.

``` r
rr_form = ei_riesz(
    ~ vap_white + vap_black + vap_other |
        state * (pop_urban + pop_rural + farm + educ_hsch + educ_coll +
                     inc_03_08k + inc_08_25k + inc_25_99k),
    data = elec_1968, total = pres_total, penalty = m_form$penalty
)
```

As with the regression model, the
[`summary()`](https://rdrr.io/r/base/summary.html) function provides
useful information for evaluating the Riesz representer.

``` r
summary(rr)
#> Second moment of representer:
#>   vap_white   vap_black   vap_other 
#>    10.62269   149.29675 10025.02281 
#> 
#> Second moment of representer (leave-one-out):
#>   vap_white   vap_black   vap_other 
#>    13.60462   194.83194 39604.31314
```

Large second moments of the Riesz representer are indicative of a more
difficult EI problem, likely due to limited variation in the predictor,
given covariates. Here we see that there is very little information for
the `other` group, and the representer is highly variable. Comparing the
in-sample and leave-one-out second moments can also help identify cases
of possible overfitting, where a higher penalty may be useful.

### DML for ecological estimates

With the regression function and Riesz representer now fitted, we are
ready to combine them to estimate our quantities of interest: vote
choice by race. This is accomplished with the
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md)
function, which takes in both fitted models and the original `ei_spec`
object, and returns a tidy data frame of estimates. The `conf_level`
argument is optional and produces confidence intervals of the specified
width from the asymptotic Normal approximation.

``` r
est = ei_est(m, rr, spec, conf_level = 0.95)
print(est)
#> # A tibble: 12 × 6
#>    predictor outcome       estimate std.error  conf.low conf.high
#>    <chr>     <chr>            <dbl>     <dbl>     <dbl>     <dbl>
#>  1 vap_white pres_dem_hum  0.234     0.0252    0.185      0.284  
#>  2 vap_black pres_dem_hum  0.608     0.0564    0.497      0.718  
#>  3 vap_other pres_dem_hum  1.13      0.246     0.649      1.61   
#>  4 vap_white pres_rep_nix  0.422     0.0358    0.352      0.492  
#>  5 vap_black pres_rep_nix -0.0651    0.0246   -0.113     -0.0168 
#>  6 vap_other pres_rep_nix -0.000759  0.142    -0.279      0.277  
#>  7 vap_white pres_ind_wal  0.342     0.0182    0.306      0.378  
#>  8 vap_black pres_ind_wal  0.457     0.0427    0.373      0.540  
#>  9 vap_other pres_ind_wal -0.130     0.252    -0.625      0.365  
#> 10 vap_white pres_abs      0.00169   0.000342  0.00102    0.00236
#> 11 vap_black pres_abs      0.00101   0.000853 -0.000663   0.00268
#> 12 vap_other pres_abs     -0.000129  0.00196  -0.00398    0.00372
```

The same call works with the formula interface.

``` r
est_form = ei_est(m_form, rr_form, elec_1968)
```

Often, a particular *contrast* is of interest, such as the difference in
vote shares between two groups. The `contrast=` argument to
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md) allows
for estimating these contrasts directly, with proper uncertainty
quantification. Here, we estimate the difference in vote shares between
predictor group 1 (White voters) and predictor group 2 (Black voters).
This is a measure of racially polarized voting.

``` r
est_c = ei_est(m, rr, spec, contrast = list(predictor = c(1, -1, 0)), conf_level = 0.95)
print(est_c)
#> # A tibble: 4 × 6
#>   predictor             outcome       estimate std.error conf.low conf.high
#>   <chr>                 <chr>            <dbl>     <dbl>    <dbl>     <dbl>
#> 1 vap_white - vap_black pres_dem_hum -0.373      0.0538  -0.479    -0.268  
#> 2 vap_white - vap_black pres_rep_nix  0.487      0.0497   0.390     0.585  
#> 3 vap_white - vap_black pres_ind_wal -0.115      0.0402  -0.194    -0.0358 
#> 4 vap_white - vap_black pres_abs      0.000683   0.00100 -0.00128   0.00265
```

Occasionally, it is helpful to examine the estimates in a different
format. The [`as.matrix()`](https://rdrr.io/r/base/matrix.html) method
works on `ei_est` objects and can be used on any column of the object,
such as the estimate or standard error. The full (asymptotic) covariance
matrix of all estimates is also accessible via
[`vcov()`](https://rdrr.io/r/stats/vcov.html).

``` r
as.matrix(est)
#>            outcome
#> predictor   pres_dem_hum  pres_rep_nix pres_ind_wal      pres_abs
#>   vap_white    0.2344005  0.4220673541    0.3418389  0.0016933258
#>   vap_black    0.6075045 -0.0650607957    0.4565462  0.0010100695
#>   vap_other    1.1309934 -0.0007588442   -0.1301076 -0.0001289748

as.matrix(est, which = "conf.low")
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal      pres_abs
#>   vap_white    0.1850358    0.3517833    0.3060911  0.0010229815
#>   vap_black    0.4967846   -0.1132949    0.3728618 -0.0006630292
#>   vap_other    0.6492101   -0.2787199   -0.6253274 -0.0039820101
```

Sometimes, estimates within a set of geographies are of interest. The
`subset=` argument to
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md) allows
for producing estimates in these smaller areas.

``` r
as.matrix(ei_est(m, rr, spec, subset = pop_city >= 0.9))
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal      pres_abs
#>   vap_white    0.2957165   0.48695283    0.2161180  0.0012127546
#>   vap_black    0.6392747  -0.05398564    0.4144446  0.0002661182
#>   vap_other    1.1385898  -0.01494923   -0.1234434 -0.0001998880
as.matrix(ei_est(m, rr, spec, subset = state == "Mississippi"))
#>            outcome
#> predictor   pres_dem_hum pres_rep_nix pres_ind_wal      pres_abs
#>   vap_white   0.02160305  0.177981594    0.8005845 -1.691338e-04
#>   vap_black   0.70403019  0.015038086    0.2803618  5.699620e-04
#>   vap_other   1.12329112  0.008639923   -0.1318805 -5.260328e-05
```

Finally,
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md)
actually also works with a regression model alone, or a Riesz
representer alone. However, these estimates are not debiased, and may
have higher error. They generally have improperly calibrated confidence
intervals.

``` r
# Not recommended
est_m = ei_est(regr = m, data = spec)
est_rr = ei_est(riesz = rr, data = spec)

sd(est_m$estimate - est_rr$estimate) # estimates (here) are close
#> [1] 7.846052e-06
sd(est_m$std.error - est_rr$std.error) # standard errors are very different
#> [1] 0.4412622
```

## Local estimates

Sometimes, it is of interest to produce estimates that are even more
fine-grained than what is possible with the `subset` argument to
[`ei_est()`](https://corymccartan.com/seine/reference/ei_est.md):
estimates for a single precinct or geography. **seine** provides two
functions for this purpose:
[`ei_bounds()`](https://corymccartan.com/seine/reference/ei_bounds.md),
which produces guaranteed-valid partial identification bounds for each
geography, and
[`ei_est_local()`](https://corymccartan.com/seine/reference/ei_est_local.md),
which produces point estimates and confidence intervals for each
geography under CAR and a few more assumptions. The former bounds are
also sometimes referred to as the Duncan–Davis bounds. See
[`vignette("local")`](https://corymccartan.com/seine/articles/local.md)
for a full walkthrough of this functionality.

## Sensitivity analysis

The entire analysis so far has rested on the critical CAR assumption. In
practice, no such independence assumption ever holds exactly. Thus, it
is important to evaluate how sensitive the results are to violations of
that identifying assumption. **seine** provides a suite of tools for
this purpose; see
[`vignette("sensitivity")`](https://corymccartan.com/seine/articles/sensitivity.md)
for a full walkthrough.

## References

McCartan, C., & Kuriwaki, S. (2025+). Identification and semiparametric
estimation of conditional means from aggregate data. Working paper
[arXiv:2509.20194](https://arxiv.org/abs/2509.20194).

------------------------------------------------------------------------

1.  A slightly weaker assumption is possible; see the methodology paper
    (McCartan and Kuriwaki 2025) for details.
