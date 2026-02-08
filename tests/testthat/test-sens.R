test_that("ei_sens works with default parameters", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(state, pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    sens = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3)

    expect_s3_class(sens, "ei_sens")
    expect_s3_class(sens, "ei_est")
    expect_true("c_outcome" %in% names(sens))
    expect_true("c_predictor" %in% names(sens))
    expect_true("bias_bound" %in% names(sens))
    expect_equal(sens$c_outcome, rep(0.5, nrow(sens)))
    expect_equal(sens$c_predictor, rep(0.3, nrow(sens)))
    expect_true(all(sens$bias_bound >= 0))
})

test_that("ei_sens works with multiple c_outcome and c_predictor values", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    c_out = c(0.1, 0.5)
    c_pred = c(0.2, 0.4, 0.6)
    sens = ei_sens(est, c_outcome = c_out, c_predictor = c_pred)

    n_combos = length(c_out) * length(c_pred)
    n_predictors = length(unique(est$predictor))
    expect_equal(nrow(sens), n_combos * n_predictors)

    expect_true(all(sens$c_outcome %in% c_out))
    expect_true(all(sens$c_predictor %in% c_pred))
})

test_that("ei_sens bias_bound parameter overrides c_predictor", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    sens = ei_sens(est, c_outcome = 0.5, bias_bound = 0.1)

    expect_true("bias_bound" %in% names(sens))
    expect_true("c_predictor" %in% names(sens))
    expect_equal(unique(sens$bias_bound), 0.1)
    expect_true(all(sens$c_predictor >= 0 & sens$c_predictor <= 1))
})

test_that("ei_sens confounding parameter affects bias bounds", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    sens1 = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, confounding = 1)
    sens0.5 = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, confounding = 0.5)

    expect_true(all(sens1$bias_bound >= sens0.5$bias_bound))
})

test_that("ei_sens expand_ci updates confidence intervals", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec, conf_level = 0.95)

    sens_expand = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, expand_ci = TRUE)
    sens_no_expand = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, expand_ci = FALSE)

    width_expand = sens_expand$conf.high - sens_expand$conf.low
    width_no_expand = sens_no_expand$conf.high - sens_no_expand$conf.low
    expect_true(all(width_expand > width_no_expand))
})

test_that("ei_sens validation", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    est_no_dml = ei_est(m, data = spec)

    expect_error(
        ei_sens(est_no_dml, c_outcome = 0.5, c_predictor = 0.3),
        "must be fit with DML"
    )

    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    expect_error(
        ei_sens(est, c_outcome = -0.1, c_predictor = 0.5),
        "must be between 0 and 1"
    )
    expect_error(
        ei_sens(est, c_outcome = 0.5, c_predictor = 1.5),
        "must be between 0 and 1"
    )
    expect_error(
        ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, confounding = -0.1),
        "must be between 0 and 1"
    )
    expect_error(
        ei_sens(est, c_outcome = 0.5, c_predictor = 0.3, confounding = c(0.5, 0.8)),
        "must be a single value"
    )
})

test_that("ei_sens preserves attributes", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    sens = ei_sens(est, c_outcome = 0.5, c_predictor = 0.3)

    expect_false(is.null(attr(sens, "sens_s")))
    expect_false(is.null(attr(sens, "bounds_inf")))
})

test_that("ei_sens works with contrasts", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec, contrast = list(predictor = c(1, -1, 0)))

    expect_no_error({
        sens = ei_sens(est, c_outcome = 0.5, c_predictor = 0.5)
    })
    expect_s3_class(sens, "ei_sens")
})

test_that("ei_sens_rv calculates robustness values", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    rv = ei_sens_rv(est, 0.1)

    expect_true("rv" %in% names(rv))
    expect_true(all(rv$rv >= 0))
    expect_true(all(rv$rv <= 1))
})

test_that("ei_sens_rv works with bias_bound as expression", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    rv = ei_sens_rv(est, 2 * std.error)

    expect_true("rv" %in% names(rv))
    expect_equal(nrow(rv), nrow(est))
})

test_that("ei_sens_rv works with confounding parameter", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    rv1 = ei_sens_rv(est, 0.1, confounding = 1)
    rv0.5 = ei_sens_rv(est, 0.1, confounding = 0.5)

    expect_true(all(rv1$rv <= rv0.5$rv))
})

test_that("ei_sens_rv validates bias_bound", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    expect_error(
        ei_sens_rv(est),
        "bias_bound.*required"
    )
    expect_error(
        ei_sens_rv(est, "invalid"),
        "must be numeric"
    )
})

test_that("ei_sens_rv preserves original data frame structure", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec)

    rv = ei_sens_rv(est, 0.1)

    expect_equal(nrow(rv), nrow(est))
    expect_true(all(names(est) %in% names(rv)))
})

test_that("ei_sens_rv works with contrasts", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))
    m = ei_ridge(spec)
    rr = ei_riesz(spec, penalty = m$penalty)
    est = ei_est(m, rr, spec, contrast = list(predictor = c(1, -1, 0)))

    expect_no_error({
        rv = ei_sens_rv(est, estimate)
    })
    expect_true("rv" %in% names(rv))
})

test_that("ei_bench produces benchmark values", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm, educ_elem))

    bench = ei_bench(spec)

    expect_s3_class(bench, "ei_bench")
    expect_true(all(
        c("covariate", "predictor", "outcome", "c_outcome", "c_predictor", "confounding") %in%
            names(bench)
    ))
})

test_that("ei_bench works with subset parameter", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    expect_no_error({
        bench = ei_bench(spec, subset = pop_urban > 0.5)
    })
    expect_s3_class(bench, "ei_bench")
})

test_that("ei_bench creates one row per covariate per predictor-outcome pair", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec)

    n_covariates = 2
    n_predictors = 3
    n_outcomes = 1
    expected_rows = n_covariates * n_predictors * n_outcomes

    expect_equal(nrow(bench), expected_rows)
})

test_that("ei_bench benchmark values are in valid ranges", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec)

    expect_true(all(bench$c_outcome >= 0 & bench$c_outcome <= 1))
    expect_true(all(bench$c_predictor >= 0 & bench$c_predictor <= 1))
    expect_true(all(bench$confounding >= -1 & bench$confounding <= 1))
})

test_that("ei_bench handles multiple outcomes", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec)

    n_outcomes = length(unique(bench$outcome))
    expect_true(n_outcomes > 1)
})


test_that("ei_bench leave-one-out logic produces different values per covariate", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec)

    covs = unique(bench$covariate)
    expect_equal(length(covs), 2)

    bench_urban = bench[bench$covariate == "pop_urban", ]
    bench_farm = bench[bench$covariate == "farm", ]
    expect_false(all(bench_urban$c_outcome == bench_farm$c_outcome))
})

test_that("ei_bench works with predictor contrast", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec, contrast = list(predictor = c(1, -1, 0)))

    expect_s3_class(bench, "ei_bench")
    expect_true(nrow(bench) > 0)
    expect_equal(length(unique(bench$predictor)), 1)
    expect_true(all(bench$c_outcome >= 0 & bench$c_outcome <= 1))
    expect_true(all(bench$c_predictor >= 0 & bench$c_predictor <= 1))
})

test_that("ei_bench works with outcome contrast", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec, contrast = list(outcome = c(1, -1, 0, 0)))

    expect_s3_class(bench, "ei_bench")
    expect_true(nrow(bench) > 0)
    expect_equal(length(unique(bench$outcome)), 1)
    expect_true(all(bench$c_outcome >= 0 & bench$c_outcome <= 1))
    expect_true(all(bench$c_predictor >= 0 & bench$c_predictor <= 1))
})

test_that("ei_bench works with predictor-outcome contrast", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec, contrast = list(predictor = c(1, -1, 0), outcome = c(1, -1, 0, 0)))

    expect_s3_class(bench, "ei_bench")
    expect_true(nrow(bench) > 0)
    expect_equal(length(unique(bench$predictor)), 1)
    expect_equal(length(unique(bench$outcome)), 1)
    expect_true(all(bench$c_outcome >= 0 & bench$c_outcome <= 1))
    expect_true(all(bench$c_predictor >= 0 & bench$c_predictor <= 1))
})

test_that("ei_bench with contrast has correct dimensions", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec, contrast = list(predictor = c(1, -1, 0)))

    n_covariates = 2
    n_contrasts = 1
    n_outcomes = 1
    expected_rows = n_covariates * n_contrasts * n_outcomes

    expect_equal(nrow(bench), expected_rows)
})

test_that("ei_bench with contrast values remain valid", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total, covariates = c(pop_urban, farm))

    bench = ei_bench(spec, contrast = list(predictor = c(1, -1, 0)))

    expect_true(all(!is.na(bench$c_outcome)))
    expect_true(all(!is.na(bench$c_predictor)))
    expect_true(all(!is.na(bench$confounding)))
    expect_true(all(is.finite(bench$c_outcome)))
    expect_true(all(is.finite(bench$c_predictor)))
})
