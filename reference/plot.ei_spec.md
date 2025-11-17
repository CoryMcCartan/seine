# Plot an EI specification

Useful for quickly visualizing scatterplots of outcome versus predictor
variables.

## Usage

``` r
# S3 method for class 'ei_spec'
plot(x, ..., pch = 16, cex = 0.2)
```

## Arguments

- x:

  An [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

- ...:

  Additional arguments passed to
  [`pairs()`](https://rdrr.io/r/graphics/pairs.html).

- pch, cex:

  As in [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

## Examples

``` r
data(elec_1968)
spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs, pres_total)
plot(spec)

```
