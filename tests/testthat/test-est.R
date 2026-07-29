test_that("Estimation methods agree when there is no penalization", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))
    m = ei_ridge(spec, penalty=0)
    rr = ei_riesz(spec, penalty=m$penalty)

    est_p = ei_est(m, data=spec)
    est_w = ei_est(rr, data=spec)
    est_d = ei_est(m, rr, data=spec)

    expect_equal(est_p$estimate, est_w$estimate)
    expect_equal(est_p$estimate, est_d$estimate)
    expect_true(all(est_d$std.error > est_p$std.error))
})

test_that("Estimation methods work with single predictor", {
    m = ei_ridge(pres_ind_wal ~ vap_black | farm, elec_1968, penalty=0)
    rr = ei_riesz(pres_ind_wal ~ vap_black | farm, elec_1968,
                  penalty=m$penalty, total=pres_total)

    est_p = ei_est(m, data=elec_1968, total=pres_total)
    est_w = ei_est(rr, data=elec_1968, total=pres_total, outcome=pres_ind_wal)
    est_d = ei_est(m, rr, data=elec_1968, total=pres_total)

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
    # expect_lt(mean(est_d$std.error), mean(est_d0$std.error)) # not always
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

test_that("Contrasts work with no predictors", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total)
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty=m$penalty)

    est0 = subset(ei_est(m, rr, spec), outcome == "pres_dem_hum")
    estc = ei_est(
        m, rr, spec,
        contrast = list(predictor = c(1, -1, 0), outcome = c(1, 0, 0, 0, 0))
    )

    pt0 = est0$estimate[1] - est0$estimate[2]
    se0 = sqrt(sum(est0$std.error[1:2]^2))
    expect_equal(estc$estimate, pt0)
    expect_lte(estc$std.error, se0)
})

test_that("TMLE estimates respect bounds [0,1] and sum to 1", {
    # Use full outcome set (pres_dem_hum:pres_abs) so outcomes sum to 1
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm))
    m  = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)

    est = ei_est(regr = m, riesz = rr, data = spec, bounds = 0:1, tmle = TRUE)
    mat = as.matrix(est)

    # estimates lie in [0, 1] within 1e-3 (targeting zeroes out the EIF residual)
    expect_true(all(mat >= -1e-3))
    expect_true(all(mat <= 1 + 1e-3))
    # estimates sum to exactly 1 across outcomes for each predictor group
    # (exact because y rows and Q* rows both sum to 1, so wtd sums to 0)
    row_sums = rowSums(mat)
    expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("TMLE requires an ei_riesz object", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total)
    m = ei_ridge(spec)
    expect_error(
        suppressWarnings(ei_est(regr = m, data = spec, bounds = 0:1, tmle = TRUE)),
        class = "rlang_error"
    )
})

test_that("TMLE targeting solves the empirical score equations", {
    y = rbind(
        c(.70, .20, .10), c(.15, .65, .20), c(.30, .25, .45),
        c(.40, .35, .25), c(.55, .10, .35), c(.20, .50, .30)
    )
    q = rbind(
        c(.40, .35, .25), c(.35, .35, .30), c(.25, .40, .35),
        c(.30, .40, .30), c(.45, .25, .30), c(.30, .30, .40)
    )
    h = cbind(seq(-1, 1, length.out = nrow(y)), c(1, -1, 1, -1, 1, -1))

    out = tmle_target(
        y,
        list(yhat = q, preds = list(group_1 = q)),
        h,
        list(group_1 = h),
        0:1,
        TRUE
    )
    score = crossprod(h, y - out$yhat) / nrow(y)

    expect_equal(score, matrix(0, ncol(h), ncol(y)), tolerance = 1e-5)
    expect_equal(rowSums(out$yhat), rep(1, nrow(y)), tolerance = 1e-12)
})

test_that("Componentwise TMLE respects finite and one-sided bounds", {
    y01 = rbind(
        c(.70, .20), c(.15, .65), c(.30, .25),
        c(.40, .35), c(.55, .10), c(.20, .50)
    )
    q01 = rbind(
        c(.40, .35), c(.35, .35), c(.25, .40),
        c(.30, .40), c(.45, .25), c(.30, .30)
    )
    h = cbind(seq(-1, 1, length.out = nrow(y01)), c(1, -1, 1, -1, 1, -1))
    rm_cf = list(group_1 = h, group_2 = h)

    check_target = function(y, q, bounds) {
        rl = list(preds = list(group_1 = q, group_2 = q), yhat = q)
        out = tmle_target(y, rl, h, rm_cf, bounds, FALSE)
        expect_equal(crossprod(h, y - out$yhat), matrix(0, 2, 2), tolerance = 1e-5)
        expect_true(all(out$yhat > bounds[1] | !is.finite(bounds[1])))
        expect_true(all(out$yhat < bounds[2] | !is.finite(bounds[2])))
        expect_true(all(out$preds[[1]] > bounds[1] | !is.finite(bounds[1])))
        expect_true(all(out$preds[[1]] < bounds[2] | !is.finite(bounds[2])))
    }

    check_target(2 + 3 * y01, 2 + 3 * q01, c(2, 5))
    check_target(2 + y01, 2 + q01, c(2, Inf))
    check_target(1 - y01, 1 - q01, c(-Inf, 1))
})

test_that("TMLE estimates respect non-unit componentwise bounds", {
    data = elec_1968
    data$pres_dem_hum = 2 + 3 * data$pres_dem_hum
    spec = ei_spec(data, vap_white:vap_other, pres_dem_hum, total = pres_total,
                   covariates = c(pop_urban, farm))
    m = suppressWarnings(ei_ridge(spec, bounds = c(2, 5), sum_one = FALSE))
    rr = ei_riesz(spec, penalty = m$penalty)

    est = ei_est(m, rr, spec, bounds = c(2, 5), sum_one = FALSE, tmle = TRUE)
    expect_true(all(est$estimate > 2))
    expect_true(all(est$estimate < 5))
})

test_that("TMLE respects outcome contrasts", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm))
    m = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)

    est = ei_est(regr = m, riesz = rr, data = spec, bounds = 0:1, tmle = TRUE)
    est_contrast = ei_est(
        regr = m, riesz = rr, data = spec, bounds = 0:1, tmle = TRUE,
        contrast = list(outcome = c(1, -1, 0, 0))
    )

    expect_equal(est_contrast$estimate,
                 unname(drop(as.matrix(est) %*% c(1, -1, 0, 0))),
                 tolerance = 1e-8)
})

test_that("TMLE works for a subset of observations", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, state))
    m = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)

    expect_no_warning({
        est = ei_est(regr = m, riesz = rr, data = spec,
                     subset = state == "Louisiana", bounds = 0:1, tmle = TRUE)
    })

    mat = as.matrix(est)
    expect_true(all(is.finite(mat)))
    expect_true(all(mat >= -1e-3))
    expect_true(all(mat <= 1 + 1e-3))
    expect_equal(unname(rowSums(mat)), rep(1, nrow(mat)), tolerance = 1e-6)
})

test_that("Counterfactual Riesz representers recover normalized fitted weights", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total, covariates = c(pop_urban, farm))
    rr = ei_riesz(spec, penalty = 0.01)
    rm_cf = tmle_riesz_cf(rr, spec)

    data = hardhat::forge(spec, rr$blueprint)$predictors
    idx_x = match(rr$blueprint$ei_x, colnames(data))
    x = pull_x(data, idx_x)
    rm_cf_array = array(unlist(rm_cf), dim = c(nrow(x), ncol(x), ncol(x)))
    rm = rowSums(
        sweep(rm_cf_array, c(1, 3), x, `*`),
        dims = 2
    )

    expect_named(rm_cf, colnames(x))
    expect_equal(
        unname(rm),
        unname(scale_cols(rr$weights, 1 / colMeans(rr$weights))),
        tolerance = 1e-8
    )
})

test_that("Simplex projection produces nonnegative unit-sum rows", {
    x = rbind(
        c(0.2, 0.3, 0.5),
        c(-1, 2, 2),
        c(-1, -1, -1),
        c(2, -1, 0)
    )
    projected = tmle_project_simplex(x)

    expect_equal(
        projected,
        rbind(c(0.2, 0.3, 0.5), c(0, 0.5, 0.5), c(1/3, 1/3, 1/3), c(1, 0, 0)),
        tolerance = 1e-7
    )
    expect_true(all(projected >= 0))
    expect_equal(rowSums(projected), rep(1, nrow(projected)), tolerance = 1e-7)
})

test_that("TMLE accepts componentwise bounds and rejects unsupported simplex bounds", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec, penalty = 0.01)
    rr = ei_riesz(spec, penalty = m$penalty)

    local_mocked_bindings(
        tmle_target = function(y, rl, rm, rm_cf, bounds, sum_one) rl,
        .package = "seine"
    )

    expect_no_error(
        ei_est(m, rr, spec, bounds = 0:1, sum_one = TRUE, tmle = TRUE)
    )
    expect_no_error(
        ei_est(m, rr, spec, bounds = c(0, 1), sum_one = NULL, tmle = TRUE)
    )
    expect_no_error(
        ei_est(m, rr, spec, bounds = c(0, 1), sum_one = FALSE, tmle = TRUE)
    )
    expect_error(tmle_fluctuation(c(-1, 1), TRUE), "sum_one = TRUE")
    expect_error(
        ei_est(m, data = spec, bounds = c(0, 1), sum_one = TRUE, tmle = TRUE),
        "requires an <ei_riesz> object"
    )
    expect_warning(
        ei_est(m, rr, spec, bounds = c(0, 1), sum_one = TRUE, tmle = FALSE),
        "only used when `tmle = TRUE`"
    )
})
