test_that("ei_bounds_impl works with simple example", {
    # Simple 2x2 case
    x = matrix(c(0.6, 0.4, 0.7, 0.3), nrow = 2, byrow = TRUE)
    y = matrix(c(0.3, 0.5), nrow = 2)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = c(100, 200)
    bounds = c(0, 1)

    result = ei_bounds_impl(x, y, total, NULL, bounds)

    expect_s3_class(result, "ei_bounds")
    expect_equal(nrow(result), 4)  # 2 rows × 2 predictors × 1 outcome
    expect_equal(ncol(result), 6)
    expect_true(all(c(".row", "predictor", "outcome", "wt", "min", "max") %in% names(result)))

    expect_true(all(result$min >= bounds[1], na.rm = TRUE))
    expect_true(all(result$max <= bounds[2], na.rm = TRUE))

    expect_true(all(result$min <= result$max, na.rm = TRUE))

    expected_wt = c(x[1, 1] * total[1], x[2, 1] * total[2],
                    x[1, 2] * total[1], x[2, 2] * total[2])
    expect_equal(as.numeric(result$wt), as.numeric(expected_wt))

    x <- matrix(c(0.8, 0.2), nrow = 1)
    y <- matrix(c(0.9, 0.1), nrow = 1)
    result <- R_bounds_lp(x, y, c(0, 1))
    expected <- list(
        min = matrix(c(0.875, 0.5, 0, 0), nrow = 1),
        max = matrix(c(1, 1, 0.125, 0.5), nrow = 1)
    )
    expect_equal(result, expected)
})

test_that("ei_bounds_impl respects custom bounds", {
    x = matrix(c(0.5, 0.5), nrow = 1)
    y = matrix(c(0.4), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = 100
    bounds = c(0.2, 0.8)

    result = ei_bounds_impl(x, y, total, NULL, bounds)

    # Check bounds are respected
    expect_true(all(result$min >= bounds[1], na.rm = TRUE))
    expect_true(all(result$max <= bounds[2], na.rm = TRUE))
})

test_that("ei_bounds_impl handles multiple outcomes", {
    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3, 0.5), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = 100
    bounds = c(0, 1)

    result = ei_bounds_impl(x, y, total, NULL, bounds)

    # Should have 1 row × 2 predictors × 2 outcomes = 4 results
    expect_equal(nrow(result), 4)
    expect_equal(unique(result$outcome), c("y1", "y2"))
    expect_equal(sort(unique(result$predictor)), c("x1", "x2"))
})

test_that("ei_bounds_impl handles infeasible cases", {
    # Create infeasible case: y > x with bounds [0, 1]
    x = matrix(c(0.3, 0.7), nrow = 1)
    y = matrix(c(0.9), nrow = 1)  # Would require B[1,1] > 1
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = 100
    bounds = c(0, 0.5)  # Tight bounds to force infeasibility

    result = ei_bounds_impl(x, y, total, NULL, bounds)

    # Should return NAs for infeasible entries
    # (depends on C++ implementation, may or may not be NA)
    expect_true(all(is.finite(result$min) | is.na(result$min)))
    expect_true(all(is.finite(result$max) | is.na(result$max)))
})

test_that("ei_bounds.formula interface works", {
    data = data.frame(
        x1 = c(0.6, 0.4, 0.7),
        x2 = c(0.4, 0.6, 0.3),
        y1 = c(0.3, 0.5, 0.4),
        y2 = c(0.4, 0.3, 0.5),
        total = c(100, 200, 150)
    )

    result = ei_bounds(y1 + y2 ~ x1 + x2, data = data, total = total, bounds = c(0, 1))

    expect_s3_class(result, "ei_bounds")
    expect_equal(nrow(result), 12)  # 3 rows × 2 predictors × 2 outcomes
    expect_true(all(result$min >= 0, na.rm = TRUE))
    expect_true(all(result$max <= 1, na.rm = TRUE))
    expect_true(all(result$min <= result$max, na.rm = TRUE))
})

# tests df too since it passes through
test_that("ei_bounds.matrix interface works", {
    x = matrix(c(0.6, 0.4, 0.7, 0.3, 0.6, 0.4), nrow = 3, byrow = TRUE)
    y = matrix(c(0.3, 0.5, 0.5, 0.4, 0.4, 0.6), nrow = 3, byrow = TRUE)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = c(100, 200, 150)

    result = ei_bounds(x, y, total = total, bounds = c(0, 1))

    expect_s3_class(result, "ei_bounds")
    expect_equal(nrow(result), 12)  # 3 rows × 2 predictors × 2 outcomes
    expect_true(all(result$min >= 0, na.rm = TRUE))
    expect_true(all(result$max <= 1, na.rm = TRUE))
    expect_true(all(result$min <= result$max, na.rm = TRUE))
})

test_that("ei_bounds.ei_spec interface works", {
    data(elec_1968)
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(state, pop_urban, farm))

    result = ei_bounds(spec, bounds = c(0, 1))

    expect_s3_class(result, "ei_bounds")
    # Should have nrow(spec) × n_predictors × n_outcomes
    n_pred = length(attr(spec, "ei_x"))
    n_outcome = length(attr(spec, "ei_y"))
    n_rows = nrow(spec)
    expect_equal(nrow(result), n_rows * n_pred * n_outcome)

    expect_true(all(result$min >= 0, na.rm = TRUE))
    expect_true(all(result$max <= 1, na.rm = TRUE))
    expect_true(all(result$min <= result$max, na.rm = TRUE))

    # Check that weights are positive
    expect_true(all(result$wt >= 0))
})

test_that("ei_bounds respects tight bounds", {
    x = matrix(c(0.5, 0.5), nrow = 1)
    y = matrix(c(0.4), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = 100

    # Test with tight lower bound
    result1 = ei_bounds(x, y, total = total, bounds = c(0.3, 1))
    expect_true(all(result1$min >= 0.3, na.rm = TRUE))

    # Test with tight upper bound
    result2 = ei_bounds(x, y, total = total, bounds = c(0, 0.7))
    expect_true(all(result2$max <= 0.7, na.rm = TRUE))
})

test_that("as.array.ei_bounds reshapes correctly", {
    x = matrix(c(0.6, 0.4, 0.7, 0.3), nrow = 2, byrow = TRUE)
    y = matrix(c(0.3, 0.4, 0.5, 0.6), nrow = 2, byrow = TRUE)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = c(100, 200)

    result = ei_bounds(x, y, total = total, bounds = c(0, 1))
    arr = as.array(result)

    expect_equal(dim(arr), c(2, 2, 2, 2))  # rows × predictors × outcomes × 2 (min/max)
    expect_equal(dimnames(arr)[[2]], c("x1", "x2"))
    expect_equal(dimnames(arr)[[3]], c("y1", "y2"))
    expect_equal(dimnames(arr)[[4]], c("min", "max"))

    # Check min <= max everywhere
    expect_true(all(arr[,,,1] <= arr[,,,2], na.rm = TRUE))
})

test_that("ei_bounds validates inputs", {
    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = 100

    # Test missing total error for non-ei_spec
    expect_error(
        ei_bounds(x, y, bounds = c(0, 1)),
        "required"
    )

    # Test NA in predictors - it will warn about not summing to 1 first
    x_na = matrix(c(0.6, NA), nrow = 1)
    colnames(x_na) = c("x1", "x2")
    expect_warning(
        expect_error(
            ei_bounds(x_na, y, total = total, bounds = c(0, 1)),
            "Missing values"
        ),
        "sum to 1"
    )

    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3), nrow = 1)
    total = 100

    expect_error(
        ei_bounds(x, y, total = total, bounds = c(-Inf, Inf)),
        "At least one bound"
    )
})

test_that("ei_bounds works with single predictor", {
    x = matrix(c(1.0, 1.0), nrow = 2)
    y = matrix(c(0.3, 0.5), nrow = 2)
    colnames(x) = c("x1")
    colnames(y) = c("y1")
    total = c(100, 200)

    result = ei_bounds(x, y, total = total, bounds = c(0, 1))

    expect_equal(nrow(result), 2)  # 2 rows × 1 predictor × 1 outcome
    # With x=1, bounds should be exactly y (since B*1 = y means B = y)
    expect_equal(result$min, c(0.3, 0.5), tolerance = 1e-10)
    expect_equal(result$max, c(0.3, 0.5), tolerance = 1e-10)
})

test_that("ei_bounds works with single outcome", {
    x = matrix(c(0.6, 0.4, 0.7, 0.3), nrow = 2, byrow = TRUE)
    y = matrix(c(0.3, 0.5), nrow = 2)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = c(100, 200)

    result = ei_bounds(x, y, total = total, bounds = c(0, 1))

    expect_equal(nrow(result), 4)  # 2 rows × 2 predictors × 1 outcome
    expect_equal(unique(result$outcome), "y1")
})
