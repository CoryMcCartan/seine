
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **seine**: Semiparametric Ecological Inference <a href="https://corymccartan.com/seine/"><img src="man/figures/logo.svg" align="right" height="144" /></a>

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
  identifying assumptions.
- **Sensitivity analysis** and **benchmarking** let researchers
  understand how violations of their assumptions will affect results.
- A **tidy interface** makes the package modular and easy to use, and
  works well with pipe-based workflows.
- **Minimal dependencies** and efficient estimation routines keep
  everything **fast** and lightweight.

<h2>

<span style="text-color: red">This software is experimental</span>
</h2>

**seine** and the underlying statistical methods are under active
development. The package is not currently recommended for non-expert
use.

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
    covariates = c(farm:educ_coll, inc_00_03k:inc_25_99k)
)

m = ei_ridge(spec)
rr = ei_riesz(spec, penalty = m$penalty)

ei_est(regr = m, riesz = rr, data = spec, conf_level = 0.95)
#> # A tibble: 9 × 6
#>   predictor outcome estimate std.error conf.low conf.high
#>   <chr>     <chr>      <dbl>     <dbl>    <dbl>     <dbl>
#> 1 white     dem_hum   0.291     0.0229   0.246     0.336 
#> 2 black     dem_hum   0.285     0.0433   0.200     0.369 
#> 3 other     dem_hum   3.80      1.14     1.57      6.03  
#> 4 white     rep_nix   0.436     0.0360   0.365     0.506 
#> 5 black     rep_nix  -0.0290    0.0313  -0.0903    0.0322
#> 6 other     rep_nix  -3.86      1.16    -6.14     -1.58  
#> 7 white     ind_wal   0.272     0.0195   0.234     0.310 
#> 8 black     ind_wal   0.740     0.0632   0.616     0.864 
#> 9 other     ind_wal   1.01      1.52    -1.97      4.00
```

This workflow is explained in more detail in the [package
vignette](./articles/seine.html), along with demonstrations of data
preprocessing and sensitivity analysis.
