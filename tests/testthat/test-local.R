test_that("local estimates satisfy constraints", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))

    m = ei_ridge(spec)

    ests <- ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE, conf_level = 0.95)
    expect_true(all(ests$estimate > -1e-12))
    expect_true(all(ests$conf.low > -1e-12))
    expect_true(all(ests$conf.high > -1e-12))
    expect_true(all(ests$estimate < 1 + 1e-12))
    expect_true(all(ests$conf.low < 1 + 1e-12))
    expect_true(all(ests$conf.high < 1 + 1e-12))

    ea = as.array(ests)
    expect_identical(dim(ea), c(nrow(spec), 3L, 4L))
    expect_lt(mean(abs(rowSums(ea, dims = 2) - 1)), 5e-16)
})


test_that("Local oblique projection is calculated correctly", {
    n_y = 2
    n_x = 3
    H = diag(n_y) %x% local_basis(c(0.0, 0.3, 0.7))
    r_cov = (diag(n_y)*1.2 - 0.2) %x% (diag(n_x)*0.5 + 0.5)
    R = chol(r_cov)

    Pi = oblique_proj(H, R) %*% chol2inv(R) # get back to proj matrix Pi for which we have tests
    expect_equal(max(abs(Pi %*% Pi - Pi)), 0) # test idempotency
    expect_equal(max(abs(Pi %*% H - H)), 0) # test range is fixed
    expect_equal(sum(diag(Pi)), ncol(H)) # test rank
    # test self-adjoint wrt Sigma
    Inv_Sigma <- solve(t(R) %*% R)
    W_Pi <- Inv_Sigma %*% Pi
    expect_equal(max(abs(W_Pi - t(W_Pi))), 0, tolerance = 1e-8)
})

test_that("local confidence intervals are wider with regression variance", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))

    m = ei_ridge(spec, vcov = TRUE)

    e1 = ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = FALSE)
    e2 = ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = TRUE)
    w1 = e1$conf.high - e1$conf.low
    w2 = e2$conf.high - e2$conf.low
    expect_true(all(w2 >= w1))
})

test_that("local confidence intervals are very narrow with neighborhood model", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))

    m = ei_ridge(spec, vcov = TRUE)

    e1 = ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = FALSE)
    e2 = ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = TRUE)
    w1 = e1$conf.high - e1$conf.low
    w2 = e2$conf.high - e2$conf.low
    expect_true(all(w2 >= w1))
})


test_that("local confidence intervals work with contrasts", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))

    m = ei_ridge(spec, vcov = TRUE)

    ests = ei_est_local(
        m,
        spec,
        contrast = list(predictor = c(1, -1, 0), outcome = c(1, -1, 0, 0)),
        bounds = c(0, 1),
        sum_one = TRUE,
        conf_level = 0.95,
        regr_var = T
    )
    expect_equal(nrow(ests), nrow(spec))
    expect_lt(mean(ests$estimate), 0)
})
