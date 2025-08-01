# fn() takes `k` arguments; test gradient in `i` at value `x`
test_gr = function(fn, x, tol=1e-10) {
    fn = rlang::as_function(fn)
    k = length(x)
    sapply(seq_len(k), function(i) {
        eps = rep(0, k)
        eps[i] = tol
        (fn(x + eps) - fn(x)) / eps[i]
    })
}

test_that("constraining Jacobian is correct", {
    for (i in 1:20) {
        upars = rnorm(5 + 2*rpois(1, 2))
        pars = constr_pars(upars)

        jac_calc = constr_jac(pars)
        jac_act = rbind(
            test_gr(~ pars_to_L(constr_pars(.))[c(1, 2, 4)], upars),
            test_gr(~ constr_pars(.)$rho, upars),
            test_gr(~ constr_pars(.)$beta, upars)
        )
        expect_equal(jac_calc, jac_act, tolerance = 1e-5)
    }
})

test_that("TMVN inference is correct", {
    b_true = c(0.5, 0.9)
    zscores = replicate(20, {
        spec = ei_synthetic(n=500, b_loc=b_true)
        m = ei_tmvn(spec)
        (m$est$beta[, 1] - b_true) / sqrt(diag(m$vcov)[4:5])
    })

    expect_true(all(abs(zscores) < 5))
})
