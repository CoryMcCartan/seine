# Methods for [`ei_ridge`](https://corymccartan.com/seine/reference/ei_ridge.md) models

Models fitted with
[`ei_ridge()`](https://corymccartan.com/seine/reference/ei_ridge.md)
support various generic methods.

## Usage

``` r
# S3 method for class 'ei_ridge'
predict(object, new_data, type = "numeric", ...)

# S3 method for class 'ei_ridge'
fitted(object, ...)

# S3 method for class 'ei_ridge'
residuals(object, ...)

# S3 method for class 'ei_ridge'
vcov(object, ...)

# S3 method for class 'ei_ridge'
summary(object, ...)

# S3 method for class 'ei_ridge'
weights(object, normalize = TRUE, ...)
```

## Arguments

- object:

  A fitted
  [ei_ridge](https://corymccartan.com/seine/reference/ei_ridge.md) model

- new_data:

  A data frame, matrix, or
  [ei_spec](https://corymccartan.com/seine/reference/ei_spec.md) of new
  predictors.

- type:

  The type of predictions to generate; only `"numeric"` is supported.

- ...:

  Additional arguments (ignored)

- normalize:

  If `TRUE`, normalize the weights to have mean 1.

## Functions

- `predict(ei_ridge)`: Predict from an `ei_ridge` model.

- `fitted(ei_ridge)`: Extract fitted values.

- `residuals(ei_ridge)`: Extract residuals.

- `vcov(ei_ridge)`: Extract unscaled covariance of coefficient
  estimates. Covariance estimate is not currently
  heteroskedasticity-robust. Multiply by `sigma2` from the fitted model
  to get the covariance matrix for a particular outcome variable.

- `summary(ei_ridge)`: Summarize the model's fitted values and \\R^2\\.

- `weights(ei_ridge)`: Extract estimation weights from a fitted model.
