# Convert counts to proportions

Divides counts in specified columns by a specified total or by their
sum, possibly storing the result in a new column. Also checks for
out-of-bounds and missing values, and can create a column containing the
remainder amount so that all proportions sum to 1.

## Usage

``` r
ei_proportions(data, ..., .total = ".total", .other = ".other", clamp = 0.001)
```

## Arguments

- data:

  A data frame.

- ...:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Columns to convert to proportions. If one column is completely
  missing, it will be imputed as 1 minus the total.

- .total:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Column to use as the total. If this column does not exist, it will be
  created as the sum of the selected columns. Supports renaming syntax.

- .other:

  \<[`tidy-select`](https://tidyselect.r-lib.org/reference/language.html)\>
  Column to store the remainder, so that the selected columns plus
  `.other` sum to 1. If the selected columns do sum to 1 after
  normalization, this argument will not be used; otherwise it will be
  created or overwritten. The calculation of `.other` is performed
  *after* clamping (see below).

- clamp:

  Proportions that are `clamp` below 0 or above 1 will be rounded to 0
  and 1, respectively. Values outside `clamp` will throw an error. Set
  `clamp=0` to disable or `clamp=Inf` to allow for out-of-bounds
  proportions (not recommended).

## Value

A modified data frame. Unselected columns are unmodified.

## Examples

``` r
data(elec_1968)
# Make a data frame with counts
d_unnorm = with(head(elec_1968, 10), data.frame(
  vap = vap,
  vap_white = vap * vap_white,
  vap_black = vap * vap_black,
  vap_other = vap * vap_other
))

ei_proportions(d_unnorm, vap_white:vap_black, .total=vap) # `.other` column created
#>      vap vap_white  vap_black vap_other       .other
#> 1  12744 0.7611425 0.23713120        22 0.0017263026
#> 2  33012 0.8595662 0.13737429       101 0.0030594935
#> 3  12370 0.6104285 0.38876314        10 0.0008084074
#> 4   7575 0.7832343 0.21570957         8 0.0010561056
#> 5  15856 0.9811428 0.01810040        12 0.0007568113
#> 6   6019 0.3912610 0.60873899         0 0.0000000000
#> 7  12341 0.6791184 0.32015234         9 0.0007292764
#> 8  55547 0.8501269 0.14738870       138 0.0024843826
#> 9  21117 0.7267604 0.27252924        15 0.0007103282
#> 10  9215 0.9290288 0.07010309         8 0.0008681498
ei_proportions(d_unnorm, vap_white:vap_other) # no total provided
#>      vap vap_white  vap_black    vap_other .total
#> 1  12744 0.7611425 0.23713120 0.0017263026  12744
#> 2  33012 0.8595662 0.13737429 0.0030594935  33012
#> 3  12370 0.6104285 0.38876314 0.0008084074  12370
#> 4   7575 0.7832343 0.21570957 0.0010561056   7575
#> 5  15856 0.9811428 0.01810040 0.0007568113  15856
#> 6   6019 0.3912610 0.60873899 0.0000000000   6019
#> 7  12341 0.6791184 0.32015234 0.0007292764  12341
#> 8  55547 0.8501269 0.14738870 0.0024843826  55547
#> 9  21117 0.7267604 0.27252924 0.0007103282  21117
#> 10  9215 0.9290288 0.07010309 0.0008681498   9215
# renaming allowed
ei_proportions(d_unnorm, white=vap_white, black=vap_black,
               .total=c(total=vap), .other="vap_other")
#>    total     white      black    vap_other
#> 1  12744 0.7611425 0.23713120 0.0017263026
#> 2  33012 0.8595662 0.13737429 0.0030594935
#> 3  12370 0.6104285 0.38876314 0.0008084074
#> 4   7575 0.7832343 0.21570957 0.0010561056
#> 5  15856 0.9811428 0.01810040 0.0007568113
#> 6   6019 0.3912610 0.60873899 0.0000000000
#> 7  12341 0.6791184 0.32015234 0.0007292764
#> 8  55547 0.8501269 0.14738870 0.0024843826
#> 9  21117 0.7267604 0.27252924 0.0007103282
#> 10  9215 0.9290288 0.07010309 0.0008681498
```
