
<!-- README.md is generated from README.Rmd. Please edit that file -->

# <span style="color: #425682"><b>seine</b></span>: Semiparametric Ecological Inference <a href="https://corymccartan.com/seine/"><img src="man/figures/logo.svg" align="right" height="144" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/CoryMcCartan/seine/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CoryMcCartan/seine/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Ecological inference (EI) is the statistical problem of learning
individual-level associations from aggregate-level data. Without certain
*identifying assumptions* and proper *estimation methods*, researchers
can easily draw incorrect conclusions from aggregate data, finding
patterns where none exist, missing important individual-level patterns,
or even concluding an effect runs in the opposite direction than it
actually does. This is known as the *ecological fallacy*.

The **seine** package allows researchers to perform modern ecological
inference quickly, accurately, and transparently.

- **Double/debiased machine learning** allows for controlling for
  confounding covariates, which increases the plausibility of
  identifying assumptions. Machine learning can be used to estimate the
  key regression model and avoid strong parametric assumptions made by
  existing EI methods.
- **Sensitivity analysis** and **benchmarking** let researchers
  understand how violations of their assumptions will affect results.
- A **tidy interface** makes the package modular and easy to use, and
  works well with pipe-based workflows.
- **Minimal dependencies** and efficient estimation routines keep
  everything **fast** and lightweight.

## Installation

You can install the development version of **seine** with

``` r
remotes::install_github("CoryMcCartan/seine")
```

## Example

The package contains county-level data from the 1968 U.S. presidential
election for southern states, where Hubert Humphrey (Democrat) ran
against Richard Nixon (Republican) and George Wallace (independent;
segregationist).

In this data, one might be interested in estimating vote choice by race.
To do so, **seine** provides both a formula interface and a tidy
interface through a new `ei_spec()` object. Here, we’ll demonstrate the
latter.

``` r
library(seine)
data(elec_1968)

spec = ei_spec(
    elec_1968, 
    predictors = vap_white:vap_other,
    outcome = pres_dem_hum:pres_ind_wal, 
    total = pres_total,
    covariates = c(state, pop_city:pop_rural, farm:educ_coll, 
                   inc_00_03k:inc_25_99k)
)

m = ei_ridge(spec)
rr = ei_riesz(spec, penalty = m$penalty)

ei_est(regr = m, riesz = rr, data = spec, conf_level = 0.95)
#> # A tibble: 9 × 6
#>   predictor outcome      estimate std.error conf.low conf.high
#>   <chr>     <chr>           <dbl>     <dbl>    <dbl>     <dbl>
#> 1 vap_white pres_dem_hum   0.225     0.0241   0.178     0.273 
#> 2 vap_black pres_dem_hum   0.584     0.0601   0.467     0.702 
#> 3 vap_other pres_dem_hum   2.92      0.744    1.46      4.38  
#> 4 vap_white pres_rep_nix   0.435     0.0365   0.363     0.506 
#> 5 vap_black pres_rep_nix  -0.0242    0.0367  -0.0963    0.0478
#> 6 vap_other pres_rep_nix  -4.75      0.991   -6.69     -2.80  
#> 7 vap_white pres_ind_wal   0.339     0.0197   0.300     0.377 
#> 8 vap_black pres_ind_wal   0.437     0.0436   0.351     0.522 
#> 9 vap_other pres_ind_wal   2.84      0.840    1.20      4.49
```

This workflow is explained in more detail on the [Get
Started](https://corymccartan.com/seine/articles/seine.html) page, along
with demonstrations of data preprocessing and sensitivity analysis.

## Name

In addition to being an acronym for “semiparametric ecological
inference,” **seine** (pronounced SAYN) refers to a type of fishing net,
which is used to catch fish like tuna in aggregate. A net also is
visually similar to the tomography lines which are commonly used to
visualize the interaction of the latent data and the accounting identity
in 2×2 ecological inference.
