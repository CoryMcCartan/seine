test_that("TMLE estimates respect bounds [0,1] and sum to 1", {
    # Use full outcome set (pres_dem_hum:pres_abs) so outcomes sum to 1
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm))
    m  = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)

    m_targeted = ei_target(m, rr, spec)
    expect_s3_class(m_targeted, "ei_targeted")
    expect_s3_class(m_targeted, "ei_wrapped")
    expect_equal(fitted(m_targeted), m_targeted$yhat)
    expect_equal(residuals(m_targeted), m_targeted$y - m_targeted$yhat)
    expect_equal(vcov(m_targeted), vcov(m))
    expect_equal(weights(m_targeted), weights(m))
    expect_no_error(capture.output(summary(m_targeted)))
    est = ei_est(regr = m_targeted, riesz = rr, data = spec)
    mat = as.matrix(est)

    # estimates lie in [0, 1] within 1e-3 (targeting zeroes out the EIF residual)
    expect_true(all(mat >= -1e-3))
    expect_true(all(mat <= 1 + 1e-3))
    # estimates sum to exactly 1 across outcomes for each predictor group
    # (exact because y rows and Q* rows both sum to 1, so wtd sums to 0)
    row_sums = rowSums(mat)
    expect_true(all(abs(row_sums - 1) < 1e-6))
    expect_true(all(est$conf.low >= 0))
    expect_true(all(est$conf.high <= 1))
})

test_that("Targeting requires an ei_riesz object", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total)
    m = ei_ridge(spec)
    expect_error(
        ei_target(m, data = spec, bounds = 0:1),
        class = "rlang_error"
    )
})

test_that("Unbounded targeting solves the empirical score equations", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total)
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)

    m_targeted = ei_target(m, rr, spec)
    rm = scale_cols(rr$weights, 1 / colMeans(rr$weights))
    expect_s3_class(m_targeted, "ei_targeted")
    expect_equal(dim(m_targeted$eps), c(ncol(rr$weights), ncol(m_targeted$y)))
    expect_equal(crossprod(rm, m_targeted$y - m_targeted$yhat),
                 matrix(0, ncol(rm), ncol(m_targeted$y)),
                 tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("Targeted regressions support sensitivity analysis", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec, bounds = 0:1)
    rr = ei_riesz(spec, penalty = m$penalty)
    m_targeted = suppressWarnings(ei_target(m, rr, spec, bounds = 0:1))
    est = ei_est(m_targeted, rr, spec)

    expect_no_error(ei_sens(est, c_outcome = 0.2, c_predictor = 0.2))
})

test_that("Targeted regressions support local regression variation", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE, vcov = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)
    m_targeted = ei_target(m, rr, spec)

    e0 = ei_est_local(m_targeted, spec, b_cov = 0, regr_var = FALSE)
    e1 = ei_est_local(m_targeted, spec, b_cov = 0, regr_var = TRUE)
    expect_true(all(e1$std.error >= e0$std.error))
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

    m_targeted = ei_target(m, rr, spec, bounds = c(2, 5), sum_one = FALSE)
    est = ei_est(m_targeted, rr, spec)
    expect_true(all(est$estimate > 2))
    expect_true(all(est$estimate < 5))
})

test_that("TMLE respects outcome contrasts", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm))
    m = suppressWarnings(ei_ridge(spec, bounds = 0:1, sum_one = TRUE))
    rr = ei_riesz(spec, penalty = m$penalty)

    m_targeted = ei_target(m, rr, spec, bounds = 0:1, sum_one = TRUE)
    est = ei_est(regr = m_targeted, riesz = rr, data = spec)
    est_contrast = ei_est(
        regr = m_targeted, riesz = rr, data = spec,
        contrast = list(outcome = c(1, -1, 0, 0))
    )

    expect_equal(est_contrast$estimate,
                 unname(drop(as.matrix(est) %*% c(1, -1, 0, 0))),
                 tolerance = 1e-8)
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

    expect_no_error(
        ei_target(m, rr, spec, bounds = 0:1, sum_one = TRUE)
    )
    expect_no_error(
        ei_target(m, rr, spec, bounds = c(0, 1), sum_one = NULL)
    )
    expect_no_error(
        ei_target(m, rr, spec, bounds = c(0, 1), sum_one = FALSE)
    )
    expect_error(tmle_fluctuation(c(-1, 1), TRUE), "sum_one = TRUE")
    expect_error(
        ei_target(m, data = spec, bounds = c(0, 1), sum_one = TRUE),
        "must be an <ei_riesz> object"
    )
})
