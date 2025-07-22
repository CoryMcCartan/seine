test_that("Estimation methods agree when there is no penalization", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))
    m = ei_ridge(spec, penalty=0)#, weights=elec_1968$pres_total)
    rr = ei_riesz(spec, penalty=m$penalty)#, weights=elec_1968$pres_total)

    est_p = ei_est(m, data=spec)
    est_w = ei_est(rr, data=spec)
    est_d = ei_est(m, rr, data=spec)

    expect_equal(est_p$estimate, est_w$estimate)
    expect_equal(est_p$estimate, est_d$estimate)
    expect_true(all(est_d$std.error > est_p$std.error))
})

test_that("Estimation methods work with character predictors", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total, covariates = c(state, pop_urban))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty=m$penalty)

    expect_no_error({
        est_p = ei_est(m, data=spec, subset = state=="Louisiana")
        est_w = ei_est(rr, data=spec, subset = state=="Louisiana")
        est_d = ei_est(m, rr, data=spec, subset = state=="Louisiana")
    })
})

test_that("Shrinkage has the expected effect", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total, covariates = c(pop_urban, farm))
    m0 = ei_ridge(spec, penalty=0)
    rr0 = ei_riesz(spec, penalty=m0$penalty)
    m = ei_ridge(spec, penalty=5)
    rr = ei_riesz(spec, penalty=m$penalty)

    est_p0 = ei_est(m0, data=spec)
    est_w0 = ei_est(rr0, data=spec)
    est_d0 = ei_est(m0, rr, data=spec)
    est_p = ei_est(m, data=spec)
    est_w = ei_est(rr, data=spec)
    est_d = ei_est(m, rr, data=spec)

    expect_lt(sd(est_p$estimate), sd(est_p0$estimate))
    expect_lt(sd(est_w$estimate), sd(est_w0$estimate))
    expect_lt(sd(est_d$estimate), sd(est_d0$estimate))

    expect_lt(mean(est_p$std.error), mean(est_p0$std.error))
    expect_lt(mean(est_w$std.error), mean(est_w0$std.error))
    # expect_lt(mean(est_d$std.error), mean(est_d0$std.error))
})

test_that("Estimation methods work with no predictors", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total)
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty=m$penalty)

    expect_no_error({
        est_p = ei_est(m, data=spec)
        est_w = ei_est(rr, data=spec)
        est_d = ei_est(m, rr, data=spec)
    })
})
