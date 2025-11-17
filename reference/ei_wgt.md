# Estimation weighting models

Several built-in helper functions to generate estimation weights from a
vector of unit totals, or an existing
[`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
object.

## Usage

``` r
ei_wgt_unif(x)

ei_wgt_prop(x)

ei_wgt_sqrt(x)
```

## Arguments

- x:

  A numeric vector of unit totals, or an existing
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

## Value

A numeric vector of estimation weights with the same number of
observations as `x`. These will have mean 1.

## Functions

- `ei_wgt_unif()`: Uniform weights across units with any population.
  Appropriate if the unit-level variance is constant, i.e.,
  homosekdastic.

- `ei_wgt_prop()`: Weights proportional to the totals. Appropriate if
  the unit-level variance is inversely proportional to the number of
  observations.

- `ei_wgt_sqrt()`: Weights proportional to the square root of the
  totals. Appropriate if the unit-level variance is inversely
  proportional to the square root of the number of observations.

## Examples

``` r
data(elec_1968)

ei_wgt_unif(head(elec_1968$pres_total))
#> [1] 1 1 1 1 1 1

spec = ei_spec(head(elec_1968), predictors = vap_white:vap_other,
               outcome = pres_ind_wal, total = pres_total)
ei_wgt_prop(spec)
#> [1] 0.8851989 2.1516032 0.9151221 0.5219934 1.0283945 0.4976879
ei_wgt_sqrt(spec)
#> [1] 0.9722265 1.5157519 0.9885225 0.7465854 1.0479170 0.7289967
```
