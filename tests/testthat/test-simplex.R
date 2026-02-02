test_that("simplex works with <= constraints only", {
    # Maximize 3x1 + 2x2
    # Subject to: x1 + x2 <= 4
    #             2x1 + x2 <= 5
    #             x1, x2 >= 0
    result <- simplex(
        a = c(3, 2),
        A1 = rbind(c(1, 1), c(2, 1)),
        b1 = c(4, 5),
        maxi = TRUE
    )
    
    expect_equal(result$solved, 1)
    expect_equal(result$soln, c(x1 = 1, x2 = 3))
    expect_equal(result$value, 9)
})

test_that("simplex works with >= constraints", {
    # Test case from explore/simplex.R
    enj <- c(200, 6000, 3000, -200)
    fat <- c(800, 6000, 1000, 400)
    vitx <- c(50, 3, 150, 100)
    vity <- c(10, 10, 75, 100)
    vitz <- c(150, 35, 75, 5)
    
    result <- simplex(
        a = enj,
        A1 = fat,
        b1 = 13800,
        A2 = rbind(vitx, vity, vitz),
        b2 = c(600, 300, 550),
        maxi = TRUE
    )
    
    expect_equal(result$solved, 1)
    expect_equal(as.numeric(result$soln), c(0, 0, 13.8, 0), tolerance = 1e-10)
    expect_equal(result$value, 41400, tolerance = 1e-10)
})

test_that("simplex works with minimization", {
    # Minimize 3x1 + 2x2
    # Subject to: x1 + x2 >= 4
    #             2x1 + x2 >= 5
    #             x1, x2 >= 0
    result <- simplex(
        a = c(3, 2),
        A2 = rbind(c(1, 1), c(2, 1)),
        b2 = c(4, 5),
        maxi = FALSE
    )
    
    expect_equal(result$solved, 1)
    expect_equal(as.numeric(result$soln), c(1, 3), tolerance = 1e-10)
    expect_equal(result$value, 9, tolerance = 1e-10)
})

test_that("simplex handles scalar constraints", {
    # Single constraint case
    result <- simplex(
        a = c(2, 3),
        A1 = c(1, 1),
        b1 = 5,
        maxi = TRUE
    )
    
    expect_equal(result$solved, 1)
    expect_true(result$value <= 15 + 1e-10) # Should be 15 at (0,5)
})

test_that("simplex returns correct structure", {
    result <- simplex(
        a = c(3, 2),
        A1 = rbind(c(1, 1), c(2, 1)),
        b1 = c(4, 5),
        maxi = TRUE
    )
    
    expect_s3_class(result, "simplex")
    expect_named(result, c("soln", "solved", "value", "maxi", "slack", "obj"))
    expect_length(result$soln, 2)
    expect_length(result$slack, 2)
    expect_true(result$maxi)
})

test_that("simplex handles slack and surplus variables correctly", {
    result <- simplex(
        a = c(1, 1),
        A1 = c(1, 0),
        b1 = 5,
        A2 = c(0, 1),
        b2 = 2,
        maxi = TRUE
    )
    
    expect_equal(result$solved, 1)
    expect_true("slack" %in% names(result))
    expect_true("surplus" %in% names(result))
})
