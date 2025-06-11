test_that("ridge regression methods agree", {
    y = mtcars$mpg
    x = model.matrix(y ~ ., mtcars)
    w = sqrt(y) / mean(sqrt(y))

    penalties = c(0, 0.1, 1, 10)
    udv = svd(x)
    udvw = svd(x * sqrt(w))
    for (lambda in penalties) {
        fit_naive = ridge_naive(x, y, rep(1, length(y)), penalty=lambda)
        fit_svd = ridge_svd(udv, y, rep(1, length(y)), penalty=lambda)

        expect_equal(fit_naive$coef, fit_svd$coef, ignore_attr=TRUE)
        expect_equal(fit_naive$fitted, fit_svd$fitted, ignore_attr=TRUE)

        fit_naive = ridge_naive(x, y, w, penalty=lambda)
        fit_svd = ridge_svd(udvw, y, sqrt(w), penalty=lambda)

        expect_equal(fit_naive$coef, fit_svd$coef, ignore_attr=TRUE)
        expect_equal(fit_naive$fitted, fit_svd$fitted, ignore_attr=TRUE)
    }
})

test_that("Riesz regression methods agree", {
    y = mtcars$mpg
    tot = y^2 / mean(y^2)
    w = sqrt(y) / mean(sqrt(y))
    z = model.matrix(y ~ ., mtcars)
    x = rbeta(length(y), 1, 3)
    x = cbind(x, 1 - x)
    xz = row_kronecker(x, z)

    penalties = c(0.0001, 0.1, 1, 10)
    udv = svd(xz * sqrt(w))
    for (lambda in penalties) {
        fit_naive = riesz_naive(xz, ncol(z), tot, w, group=1, penalty=lambda)
        fit_svd = riesz_svd(xz, udv, ncol(z), tot, w, sqrt(w), group=1, penalty=lambda)

        expect_equal(fit_naive$alpha, fit_svd$alpha, ignore_attr=TRUE)
        expect_equal(fit_naive$loo, fit_svd$loo, ignore_attr=TRUE)
    }
})

test_that("leave-one-out shortcut is correct for ridge", {
    y = mtcars$mpg
    x = model.matrix(y ~ ., mtcars)
    w = sqrt(y) / mean(sqrt(y))
    udv = svd(x * sqrt(w))
    pen = 0.2

    loo_act = numeric(length(y))
    for (i in seq_along(y)) {
        loo_act[i] = y[i] - x[i, ] %*% ridge_naive(x[-i, ], y[-i], w[-i], penalty=pen)$coef
    }

    fit_naive = ridge_naive(x, y, w, penalty=pen)
    fit_svd = ridge_svd(udv, y, sqrt(w), penalty=pen)
    hat_naive = ridge_hat_naive(x, w, penalty=pen)
    hat_svd = ridge_hat_svd(udv, pen)

    expect_equal(fit_naive$fitted, fit_svd$fitted, ignore_attr=TRUE)
    expect_equal(hat_naive, hat_svd, ignore_attr=TRUE)

    loo_hat = c(y - fit_svd$fitted) / (1 - hat_svd)

    expect_equal(loo_act, loo_hat)
})

test_that("leave-one-out shortcut is correct for Riesz regression", {
    skip("Need to fix LOO weights")
    y = mtcars$mpg
    w = sqrt(y) / mean(sqrt(y))
    z = model.matrix(y ~ ., mtcars)
    x = rbeta(length(y), 1, 3)
    x = cbind(x, 1 - x)
    xz = row_kronecker(x, z, int_scale = 10)
    pen = 0.2

    loo_act = numeric(length(y))
    for (i in seq_along(y)) {
        wr = w
        wr[i] = 0
        loo_act[i] = riesz_naive(xz, ncol(z), wr, penalty=pen)$alpha[i]
        # loo_act[i] = y[i] - x[i, ] %*% riesz_naive(xz[-i, ], ncol(z), w[-i], penalty=pen)$coef
    }

    fit_naive = riesz_naive(xz, ncol(z), w, penalty=pen)
    plot(fit_naive$loo, loo_act)
    expect_equal(fit_naive$loo, loo_act)
})
