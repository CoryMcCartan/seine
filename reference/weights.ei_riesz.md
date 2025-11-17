# Extract Riesz representer weights

Extracts a single set of Riesz representer weights from an `ei_riesz`
object, for a selected group.

## Usage

``` r
# S3 method for class 'ei_riesz'
weights(object, group = TRUE, loo = FALSE, ...)
```

## Arguments

- object:

  An
  [`ei_riesz()`](https://corymccartan.com/seine/reference/ei_riesz.md)
  object.

- group:

  The group for which to extract the weights, as a numeric index or a
  character column name. The special (default) value `TRUE` will return
  a matrix of weights, with each column corresponding to a group.

- loo:

  If `TRUE`, return the leave-one-out weights

- ...:

  Additional arguments (ignored)

## Value

A numeric vector of weights
