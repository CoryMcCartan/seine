test_that("ei_bounds_bridge works with simple example", {
    # Simple 2x2 case
    x = matrix(c(0.6, 0.4, 0.7, 0.3), nrow = 2, byrow = TRUE)
    y = matrix(c(0.3, 0.5), nrow = 2)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = c(100, 200)
    bounds = c(0, 1)

    result = ei_bounds_bridge(x, y, total, NULL, bounds, NULL)

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

    # Test calling R_bounds_lp directly (old C++ API with 3 arguments)
    x_test <- matrix(c(0.8, 0.2), nrow = 1)
    y_test <- matrix(c(0.9, 0.1), nrow = 1)
    result <- R_bounds_lp(x_test, y_test, c(0, 1))
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
    bounds = c(0.2, 0.8)

    # Create identity contrast matrix for no contrast
    contrast_mat = diag(2)  # 2 = 1 outcome * 2 predictors

    result = ei_bounds_impl(x, y, contrast_mat, bounds, FALSE)

    # Check bounds are respected
    expect_true(all(result$min >= bounds[1], na.rm = TRUE))
    expect_true(all(result$max <= bounds[2], na.rm = TRUE))
})

test_that("ei_bounds_bridge handles multiple outcomes", {
    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3, 0.5), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = 100
    bounds = c(0, 1)

    result = ei_bounds_bridge(x, y, total, NULL, bounds, NULL)

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
    bounds = c(0, 0.5)  # Tight bounds to force infeasibility

    # Create identity contrast matrix for no contrast
    contrast_mat = diag(2)  # 2 = 1 outcome * 2 predictors

    result = ei_bounds_impl(x, y, contrast_mat, bounds, FALSE)

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

    expect_true(all(result$min >= 0 - 1e-8, na.rm = TRUE))
    expect_true(all(result$max <= 1 + 1e-8, na.rm = TRUE))
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

test_that("global bounds work", {
    x = matrix(c(0.6, 0.4, 0.7, 0.3), nrow = 2, byrow = TRUE)
    y = matrix(c(0.3, 0.5), nrow = 2)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = c(100, 200)

    result = ei_bounds(x, y, total = total, bounds = c(0, 1), global = TRUE)

    expect_equal(nrow(result), 2)  # 2 predictors × 1 outcome
})

test_that("contrast matching exact linear combination yields tight bounds", {
    # When the contrast matches the x proportions that define the constraint,
    # the bounds should be exactly equal to y (tight bounds)
    # For a 2-predictor problem with x=c(0.4, 0.6) and y=0.8,
    # the constraint is 0.4*B[1] + 0.6*B[2] = 0.8
    # So the contrast c(0.4, 0.6) should yield bounds of exactly [0.8, 0.8]

    x = matrix(c(0.4, 0.6), nrow = 1)
    y = matrix(c(0.8), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1")
    total = 100
    contrast = list(predictor = c(0.4, 0.6))

    result = ei_bounds(x, y, total = total, bounds = c(0, 1), contrast = contrast)

    # Should have 1 row (for the contrast, not 2 rows for individual predictors)
    expect_equal(nrow(result), 1)
    # Bounds should be tight at y=0.8
    expect_equal(result$min, 0.8, tolerance = 1e-10)
    expect_equal(result$max, 0.8, tolerance = 1e-10)
})

test_that("contrast with sum_one constraint yields tight bounds for c(1,1,1)", {
    # With 3 outcomes and a contrast of c(1, 1, 1), if we have a sum_one
    # constraint (entries sum to 1 across outcomes within predictors),
    # the contrast should yield tight bounds of [1, 1]
    # This is because B[i,1] + B[i,2] + B[i,3] = 1 for each predictor i,
    # so summing across all 3 outcomes must equal 1

    x = matrix(c(0.5, 0.3, 0.2), nrow = 1)
    y = matrix(c(0.2, 0.5, 0.3), nrow = 1)
    colnames(x) = c("x1", "x2", "x3")
    colnames(y) = c("y1", "y2", "y3")
    total = 100
    contrast = list(outcome = c(1, 1, 1))

    result = ei_bounds(x, y, total = total, bounds = c(0, 1),
                       contrast = contrast, sum_one = TRUE)

    # Should have 1 row for each predictor (3 total)
    expect_equal(nrow(result), 3)
    # All bounds should be tight at 1
    expect_equal(result$min, c(1, 1, 1), tolerance = 1e-10)
    expect_equal(result$max, c(1, 1, 1), tolerance = 1e-10)

    # near-sum-to-1 also works
    data(elec_1968)
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                total = pres_total, covariates = c(state, pop_urban, farm))
    result = ei_bounds(spec, bounds = c(0, 1), global = TRUE, contrast = list(outcome = c(1, 1, 1, 0)))
    expect_true(all(result$max <= 1 + 1e-8))
    expect_true(all(result$min >= 0.9))
})

test_that("contrast with difference c(1, -1) computes bounds on differences", {
    # Test a contrast that computes y1 - y2
    # If y1 = 0.3, y2 = 0.5, x = c(0.6, 0.4)
    # The constraints are:
    #   0.6*B[1,1] + 0.4*B[2,1] = 0.3
    #   0.6*B[1,2] + 0.4*B[2,2] = 0.5
    # The contrast c(1, -1) gives bounds on B[i,1] - B[i,2]
    # With bounds [0, 1], the difference can range from -1 to 1

    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3, 0.5), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = 100
    contrast = list(outcome = c(1, -1))

    result = ei_bounds(x, y, total = total, bounds = c(0, 1), contrast = contrast)

    # Should have 2 rows (one per predictor)
    expect_equal(nrow(result), 2)
    # Min should be at least -1, max at most 1
    expect_true(all(result$min >= -1))
    expect_true(all(result$max <= 1))
    expect_true(all(result$min <= result$max))
})

test_that("trivial contrast c(1, 0, 0) matches original bounds for first outcome", {
    # A contrast like c(1, 0, 0) should give the same bounds as the
    # original bounds for the first outcome (y1) only

    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3, 0.5, 0.2), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2", "y3")
    total = 100
    contrast = list(outcome = c(1, 0, 0))

    result_contrast = ei_bounds(x, y, total = total, bounds = c(0, 1), contrast = contrast)
    result_original = ei_bounds(x, y, total = total, bounds = c(0, 1))

    # Filter original results to just y1
    result_y1 = result_original[result_original$outcome == "y1", ]

    # Should have same number of rows (2 predictors)
    expect_equal(nrow(result_contrast), 2)
    expect_equal(nrow(result_y1), 2)

    # Bounds should match
    expect_equal(result_contrast$min, result_y1$min, tolerance = 1e-10)
    expect_equal(result_contrast$max, result_y1$max, tolerance = 1e-10)
})

test_that("contrast with multiple outcomes respects the linear combination", {
    # Test a non-trivial contrast with 2 outcomes
    # If y1 = 0.3, y2 = 0.5, x = c(0.6, 0.4)
    # The constraints are:
    #   0.6*B[1,1] + 0.4*B[2,1] = 0.3
    #   0.6*B[1,2] + 0.4*B[2,2] = 0.5
    # A contrast of c(2, 1) should give bounds on 2*B[i,1] + 1*B[i,2]

    x = matrix(c(0.6, 0.4), nrow = 1)
    y = matrix(c(0.3, 0.5), nrow = 1)
    colnames(x) = c("x1", "x2")
    colnames(y) = c("y1", "y2")
    total = 100
    contrast = list(outcome = c(2, 1))  # 2*y1 + 1*y2

    result = ei_bounds(x, y, total = total, bounds = c(0, 1), contrast = contrast)

    # Should have 2 rows (one per predictor)
    expect_equal(nrow(result), 2)
    # The bounds should respect the linear combination
    expect_true(all(result$min <= result$max))
    # Minimum possible: 2*0 + 1*0 = 0
    expect_true(all(result$min >= 0))
    # Maximum possible: 2*1 + 1*1 = 3
    expect_true(all(result$max <= 3))
})
