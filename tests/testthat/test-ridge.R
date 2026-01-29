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
    y = mtcars$mpg
    tot = sqrt(y) / mean(sqrt(y))
    # w = rep(1, length(y))
    w = 1/tot
    z = model.matrix(y ~ ., mtcars)
    x = rbeta(length(y), 1, 3)
    x = cbind(x, 1 - x)
    xz = row_kronecker(x, z, int_scale = 10)
    pen = 0.2

    loo_act = numeric(length(y))
    for (i in seq_along(y)) {
        wr = w
        tr = tot
        wr[i] = 0
        tr[i] = 0
        loo_act[i] = riesz_naive(xz, ncol(z), tr, wr, penalty=pen)$alpha[i]
    }

    fit_naive = riesz_naive(xz, ncol(z), tot, w, penalty=pen)
    # plot(fit_naive$loo, loo_act)
    # mean(abs(fit_naive$loo - loo_act)) / mean(abs(loo_act))
    expect_equal(fit_naive$loo, loo_act, tolerance = 0.2)
})

test_that("ridge constraints work", {
    d = elec_1968
    form = pres_dem_hum + pres_rep_nix + pres_ind_wal + pres_abs ~ vap_white  |
        pop_urban + pop_rural + farm + educ_elem + educ_hsch + educ_coll +
        inc_00_03k + inc_03_08k + inc_08_25k + inc_25_99k + log(pop) + pres_turn

    m = ei_ridge(form, data=elec_1968)
    m01 = ei_ridge(form, data=elec_1968, bounds=0:1, sum_one=FALSE)
    m01s = ei_ridge(form, data=elec_1968, bounds=c(0, 1), sum_one=TRUE)
    m01def = ei_ridge(form, data=elec_1968, bounds=NULL, sum_one=NULL)

    expect_true(min(fitted(m)) < 0)
    expect_true(min(fitted(m01)) > -.Machine$double.eps)
    expect_true(all(ei_est(m01, data=elec_1968, total=pres_total)$estimate > 0))
    expect_true(all(ei_est(m01, data=elec_1968, total=pres_total)$estimate < 1))

    expect_true(min(fitted(m01s)) > -.Machine$double.eps)
    expect_true(all(ei_est(m01s, data=elec_1968, total=pres_total)$estimate > 0))
    expect_true(all(ei_est(m01s, data=elec_1968, total=pres_total)$estimate < 1))

    tots = rowSums(as.matrix(ei_est(m01s, data=elec_1968, total=pres_total)))
    expect_true(all.equal(tots, c(vap_white=1, .other=1)))
    expect_identical(m01def, m01s) # check defaults infer correctly
})
