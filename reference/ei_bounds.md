# Compute partial identification bounds on local ecological quantities

For each observation, computes the minimum and maximum value of each
local estimand that is consistent with the accounting identity and the
bounds on the outcome. Optionally aggregate the computed bounds across
units.

## Usage

``` r
ei_bounds(
  x,
  ...,
  total,
  contrast = NULL,
  bounds = c(0, 1),
  sum_one = NULL,
  global = FALSE
)

# S3 method for class 'ei_spec'
ei_bounds(
  x,
  total,
  contrast = NULL,
  bounds = c(0, 1),
  sum_one = NULL,
  global = FALSE,
  ...
)

# S3 method for class 'formula'
ei_bounds(
  formula,
  data,
  total,
  contrast = NULL,
  bounds = c(0, 1),
  sum_one = NULL,
  global = FALSE,
  ...
)

# S3 method for class 'data.frame'
ei_bounds(
  x,
  y,
  total,
  contrast = NULL,
  bounds = c(0, 1),
  sum_one = NULL,
  global = FALSE,
  ...
)

# S3 method for class 'matrix'
ei_bounds(
  x,
  y,
  total,
  contrast = NULL,
  bounds = c(0, 1),
  sum_one = NULL,
  global = FALSE,
  ...
)

# Default S3 method
ei_bounds(x, ...)

# S3 method for class 'ei_bounds'
as.array(x, ...)
```

## Arguments

- x:

  An object of class `ei_bounds`

- ...:

  Additional arguments (ignored)

- total:

  \<[`data-masking`](https://rlang.r-lib.org/reference/args_data_masking.html)\>
  A variable containing the total number of observations in each
  aggregate unit. For example, the column containing the total number of
  voters. Required for computing weights unless `x` is an
  [`ei_spec()`](https://corymccartan.com/seine/reference/ei_spec.md)
  object.

- contrast:

  If provided, a list containing entries `predictor` and `outcome`, each
  containing a contrast vector. If only one of `predictor` or `outcome`
  is provided, the contrast will be calculated for all levels of the
  other variable. For example `list(predictor = c(1, -1, 0))` will
  calculate the difference in each outcome between the first and second
  predictor groups; `list(outcome = c(1, -1))` will calculate the
  difference between the two outcomes for each predictor group; and
  `list(predictor = c(1, -1, 0), outcome = c(1, -1))` will calculate the
  difference in differences.

- bounds:

  A vector `c(min, max)` of bounds for the outcome. If `bounds = NULL`,
  they will be inferred from the outcome variable: if it is contained
  within \\\[0, 1\]\\, for instance, then the bounds will be `c(0, 1)`.
  The default `bounds = FALSE` uses an unbounded outcome.

- sum_one:

  If `TRUE`, the outcome variables are constrained to sum to one within
  each predictor group. Can only apply when `bounds` are enforced and
  there is more than one outcome variable. If `NULL`, infers
  `sum_one = TRUE` when the bounds are `c(0, 1)` and the outcome
  variables sum to 1.

- global:

  If `TRUE`, aggregate the bounds across units to produce bounds on the
  global estimands.

- formula:

  A formula such as `y ~ x0 + x1` specifying the outcome `y` and the
  predictors of interest `x`. The predictors should form a partition,
  that is, `x0 + x1 = 1` for each observation. Users can be include more
  than two predictors as well, e.g.
  `pct_white + pct_black + pct_hisp + pct_other`. If there are just two
  predictors, it is acceptable to only include one in the formula; the
  other will be formed as 1 minus the provided predictor.

- data:

  When a **formula** is used, `data` is a **data frame** containing both
  the predictors and the outcome.

- y:

  When `x` is a **data frame** or **matrix**, `y` is the outcome
  specified as:

  - A **data frame** with numeric columns.

  - A **matrix**

  - A numeric **vector**.

  When the outcome is a proportion, you can use
  [`ei_proportions()`](https://corymccartan.com/seine/reference/ei_proportions.md)
  to assist in preparing it.

## Value

A data frame with bounds. The `.row` column in the output corresponds to
the observation index in the input. The `min` and `max` columns contain
the minimum and maximum values for each local estimand. The `weight`
column contains the product of the predictor variable and total for each
observation, where applicable. Taking a weighted average of the bounds
against this column will produce global bounds. It has class
`ei_bounds`.

## Methods (by generic)

- `as.array(ei_bounds)`: Format bounds as an array with dimensions
  `<rows>*<predictors>*<outcomes>*2`. Does not work if the object has
  been sorted.

## Examples

``` r
data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm))

ei_bounds(spec, bounds = c(0, 1))
#> # A tibble: 13,716 × 6
#>     .row predictor outcome      weight     min    max
#>    <int> <chr>     <chr>         <dbl>   <dbl>  <dbl>
#>  1     1 vap_white pres_dem_hum  5877. 0       0.262 
#>  2     2 vap_white pres_dem_hum 16131. 0       0.122 
#>  3     3 vap_white pres_dem_hum  4872. 0       0.397 
#>  4     4 vap_white pres_dem_hum  3566. 0       0.180 
#>  5     5 vap_white pres_dem_hum  8801. 0.0190  0.0382
#>  6     6 vap_white pres_dem_hum  1698. 0       1     
#>  7     7 vap_white pres_dem_hum  4970. 0       0.250 
#>  8     8 vap_white pres_dem_hum 22844. 0.00707 0.183 
#>  9     9 vap_white pres_dem_hum  7731. 0       0.178 
#> 10    10 vap_white pres_dem_hum  5259. 0.00915 0.0855
#> # ℹ 13,706 more rows
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

# Infer bounds
ei_bounds(pres_ind_wal ~ vap_white, data = elec_1968, total = pres_total, bounds = NULL)
#> # A tibble: 2,286 × 6
#>     .row predictor outcome      weight   min   max
#>    <int> <chr>     <chr>         <dbl> <dbl> <dbl>
#>  1     1 vap_white pres_ind_wal  5877. 0.620 0.934
#>  2     2 vap_white pres_ind_wal 16131. 0.726 0.889
#>  3     3 vap_white pres_ind_wal  4872. 0.487 1    
#>  4     4 vap_white pres_ind_wal  3566. 0.743 1    
#>  5     5 vap_white pres_ind_wal  8801. 0.721 0.741
#>  6     6 vap_white pres_ind_wal  1698. 0     1    
#>  7     7 vap_white pres_ind_wal  4970. 0.644 1    
#>  8     8 vap_white pres_ind_wal 22844. 0.668 0.844
#>  9     9 vap_white pres_ind_wal  7731. 0.640 1    
#> 10    10 vap_white pres_ind_wal  5259. 0.828 0.905
#> # ℹ 2,276 more rows

# With contrast
ei_bounds(
    spec,
    bounds = c(0, 1),
    contrast = list(predictor = c(1, -1, 0), outcome = c(1, -1, 0, 0))
)
#> # A tibble: 1,143 × 5
#>     .row predictor             outcome                        min   max
#>    <int> <chr>                 <chr>                        <dbl> <dbl>
#>  1     1 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.942 0.588
#>  2     2 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.895 0.958
#>  3     3 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.704 0.523
#>  4     4 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.727 0.445
#>  5     5 vap_white - vap_black pres_dem_hum - pres_rep_nix -1.21  0.831
#>  6     6 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.854 0.973
#>  7     7 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.633 0.466
#>  8     8 vap_white - vap_black pres_dem_hum - pres_rep_nix -1.13  0.953
#>  9     9 vap_white - vap_black pres_dem_hum - pres_rep_nix -0.614 0.551
#> 10    10 vap_white - vap_black pres_dem_hum - pres_rep_nix -1.06  0.950
#> # ℹ 1,133 more rows

# manually aggregate min/max
# easier with dplyr:
# summarize(across(min:max, ~ weighted.mean(.x, weight)), .by=c(predictor, outcome))
grp_units = split(ei_bounds(spec, bounds = c(0, 1)), ~ predictor + outcome)
do.call(rbind, lapply(grp_units, function(b) {
    tibble::tibble(
        predictor = b$predictor[1],
        outcome = b$outcome[1],
        min = weighted.mean(b$min, b$weight),
        max = weighted.mean(b$max, b$weight)
    )
}))
#> # A tibble: 12 × 4
#>    predictor outcome            min     max
#>  * <chr>     <chr>            <dbl>   <dbl>
#>  1 vap_black pres_abs     0         0.00844
#>  2 vap_other pres_abs     0         0.0910 
#>  3 vap_white pres_abs     0.0000108 0.00179
#>  4 vap_black pres_dem_hum 0.00467   0.927  
#>  5 vap_other pres_dem_hum 0         1      
#>  6 vap_white pres_dem_hum 0.176     0.375  
#>  7 vap_black pres_ind_wal 0.00874   0.946  
#>  8 vap_other pres_ind_wal 0         1.000  
#>  9 vap_white pres_ind_wal 0.213     0.415  
#> 10 vap_black pres_rep_nix 0         0.808  
#> 11 vap_other pres_rep_nix 0         0.991  
#> 12 vap_white pres_rep_nix 0.247     0.421  
```
