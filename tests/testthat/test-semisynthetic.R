data(elec_1968)

# A single-outcome spec and a multi-outcome spec used across tests
spec1 = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
                total = pres_total, covariates = c(pop_urban, farm))
spec4 = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                total = pres_total, covariates = c(pop_urban, farm))

# Pre-fit objects shared across tests (fitted once to speed up the suite)
m1   = ei_ridge(spec1)
rr1  = ei_riesz(spec1, penalty = m1$penalty)
est1 = ei_est(m1, rr1, spec1)
bc1  = ei_local_cov(m1, spec1)

m4   = ei_ridge(spec4)
rr4  = ei_riesz(spec4, penalty = m4$penalty)
est4 = ei_est(m4, rr4, spec4)
bc4  = ei_local_cov(m4, spec4)

# --------------------------------------------------------------------------
# Structure tests
# --------------------------------------------------------------------------

test_that("ei_semisynthetic returns a correctly formatted ei_spec", {
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    expect_s3_class(ss, "ei_spec")
    expect_true(is.data.frame(ss))

    expected_cols = c(attr(spec1, "ei_x"), attr(spec1, "ei_y"),
                      attr(spec1, "ei_z"), ".confounder")
    expect_setequal(names(ss), expected_cols)
    expect_true(".confounder" %in% names(ss))

    expect_false(".confounder" %in% attr(ss, "ei_z"))
    expect_equal(attr(ss, "ei_z"), attr(spec1, "ei_z"))

    expect_equal(nrow(ss), nrow(spec1))
})


test_that("ei_semisynthetic works with multiple outcomes", {
    expect_no_error({
        ss = ei_semisynthetic(m4, rr4, spec4, est = est4, b_cov = bc4,
                              c_outcome = 0.2, c_predictor = 0.1)
    })
    expect_s3_class(ss, "ei_spec")
    expect_setequal(attr(ss, "ei_y"), attr(spec4, "ei_y"))
})

# --------------------------------------------------------------------------
# Predictor and covariate preservation
# --------------------------------------------------------------------------

test_that("predictor columns are unchanged", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    x_nm = attr(spec1, "ei_x")
    expect_equal(as.matrix(ss[, x_nm]), as.matrix(spec1[, x_nm]))
})

test_that("covariate columns are unchanged", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    z_nm = attr(spec1, "ei_z")
    expect_equal(as.matrix(ss[, z_nm]), as.matrix(spec1[, z_nm]))
})

test_that("totals are unchanged", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    expect_equal(attr(ss, "ei_n"), attr(spec1, "ei_n"))
})

test_that("outcome columns are replaced with new simulated values", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    y_nm = attr(spec1, "ei_y")
    expect_false(identical(ss[[y_nm]], spec1[[y_nm]]))
})

# --------------------------------------------------------------------------
# Confounder statistical properties
# --------------------------------------------------------------------------

tol = 0.001  # tolerance for statistical property checks

# Helper: reconstruct z_riesz exactly as the implementation does
make_z_riesz = function(spec, rr) {
    z_proc = as.matrix(attr(spec, "ei_z_proc"))
    scale_cols(shift_cols(z_proc, rr$z_shift), rr$z_scale)
}

test_that("actual c_predictor matches requested for specific predictor", {
    set.seed(1968)
    target = 0.15
    pred = "vap_white"
    ss = ei_semisynthetic(m1, rr1, spec1, predictor = pred, outcome = "pres_ind_wal",
                          est = est1, b_cov = bc1,
                          c_outcome = 0.1, c_predictor = target)
    a = ss$.confounder

    z_riesz = make_z_riesz(spec1, rr1)
    x = as.matrix(spec1[attr(spec1, "ei_x")])
    total = attr(spec1, "ei_n")
    j = match(pred, attr(spec1, "ei_x"))

    z_aug = cbind(z_riesz, a)
    fit_aug = ei_riesz_impl(x, z_aug, total, rep(1, nrow(spec1)), rr1$penalty)

    actual = 1 - rr1$nu2[j] / fit_aug$nu2[j]
    expect_lt(abs(actual - target), tol)
})

test_that("c_predictor holds for specific predictor in multi-outcome spec", {
    set.seed(1968)
    target = 0.12
    pred = "vap_black"
    ss = ei_semisynthetic(m4, rr4, spec4, predictor = pred, outcome = "pres_dem_hum",
                          est = est4, b_cov = bc4,
                          c_outcome = 0.1, c_predictor = target)
    a = ss$.confounder

    z_riesz = make_z_riesz(spec4, rr4)
    x = as.matrix(spec4[attr(spec4, "ei_x")]); x = x / rowSums(x)
    total = attr(spec4, "ei_n")
    j = match(pred, attr(spec4, "ei_x"))

    z_aug = cbind(z_riesz, a)
    fit_aug = ei_riesz_impl(x, z_aug, total, rep(1, nrow(spec4)), rr4$penalty)

    actual = 1 - rr1$nu2[j] / fit_aug$nu2[j]
    expect_lt(abs(actual - target), tol)
})

test_that("actual c_outcome matches requested for specific outcome", {
    set.seed(1968)
    target = 0.2
    out = "pres_ind_wal"
    ss = ei_semisynthetic(m1, rr1, spec1, predictor = "vap_white", outcome = out,
                          est = est1, b_cov = bc1,
                          c_outcome = target, c_predictor = 0.1)
    a = ss$.confounder
    resid_y = as.matrix(ss[out]) - m1$fitted[, out, drop = FALSE]
    r2 = cor(as.vector(resid_y), a)^2
    expect_lt(abs(r2 - target), tol)
})

test_that("c_outcome holds for specific outcome in multi-outcome spec; others unaffected", {
    set.seed(1968)
    target = 0.15
    out = "pres_rep_nix"
    ss = ei_semisynthetic(m4, rr4, spec4, predictor = "vap_white", outcome = out,
                          est = est4, b_cov = bc4,
                          c_outcome = target, c_predictor = 0.1)
    a = ss$.confounder

    resid_target = as.matrix(ss[out]) - m4$fitted[, out, drop = FALSE]
    r2 = cor(as.vector(resid_target), a)^2
    expect_lt(abs(r2 - target), tol)

    # Other outcomes should NOT be shifted by the confounder
    out_other = "pres_dem_hum"
    resid_other = as.matrix(ss[out_other]) - m4$fitted[, out_other, drop = FALSE]
    r2_other = cor(as.vector(resid_other), a)^2
    expect_lt(r2_other, tol)
})

test_that("confounding = 0 produces no outcome shift", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, predictor = "vap_white", outcome = "pres_ind_wal",
                          est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.15, confounding = 0)
    # With confounding = 0, lambda_eff = 0, so eta_new = eta (no confounder shift)
    expect_equal(attr(ss, "eta"), attr(ss, "b_loc"), ignore_attr = TRUE)
})

test_that("confounding matches squared correlation of fitted and RR changes", {
    set.seed(1968)
    target = 0.35
    pred = "vap_white"
    out = "pres_ind_wal"
    ss = ei_semisynthetic(m1, rr1, spec1, predictor = pred, outcome = out,
                          est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.15, confounding = target)
    a = ss$.confounder

    z_riesz = make_z_riesz(spec1, rr1)
    x = as.matrix(spec1[attr(spec1, "ei_x")]); x = x / rowSums(x)
    total = attr(spec1, "ei_n")
    j = match(pred, attr(spec1, "ei_x"))
    k = match(out, attr(spec1, "ei_y"))

    fit_aug = ei_riesz_impl(x, cbind(z_riesz, a), total, rep(1, nrow(spec1)), rr1$penalty)

    delta_fit = attr(ss, "eta")[, j + (k - 1L) * ncol(x)] -
        attr(ss, "b_loc")[, j + (k - 1L) * ncol(x)]
    delta_rr = fit_aug$alpha[, j] - rr1$weights[, j]
    actual = cor(delta_fit, delta_rr)^2

    expect_lt(abs(actual - target), tol)
})

test_that("bias bound holds for specific predictor/outcome combo", {
    set.seed(1968)
    c_out = 0.15; c_pred = 0.10; conf = 0.5
    ss = ei_semisynthetic(m1, rr1, spec1, predictor = "vap_white", outcome = "pres_ind_wal",
                          est = est1, b_cov = bc1,
                          c_outcome = c_out, c_predictor = c_pred, confounding = conf)

    est_ss = ei_est(m1, rr1, ss)
    true_vals = attr(ss, "est_true")$true[1]
    actual_bias = abs(est_ss$estimate[1] - true_vals)

    sens = ei_sens(est1, c_outcome = c_out, c_predictor = c_pred, confounding = conf)
    expect_true(actual_bias <= sens$bias_bound[1])
})

# --------------------------------------------------------------------------
# Contrast handling
# --------------------------------------------------------------------------

test_that("c_predictor calibration works with a predictor contrast", {
    set.seed(1968)
    target = 0.12
    q = c(1, -1, 0)  # vap_white - vap_black
    contr = list(predictor = q)
    pred_nm = check_contrast(contr, attr(spec4, "ei_x"), attr(spec4, "ei_y"))$x_nm[1]

    ss = ei_semisynthetic(m4, rr4, spec4, predictor = pred_nm, outcome = "pres_dem_hum",
                          est = ei_est(m4, rr4, spec4, contrast = contr),
                          b_cov = bc4,
                          c_outcome = 0.1, c_predictor = target, contrast = contr)
    a = ss$.confounder

    z_riesz = make_z_riesz(spec4, rr4)
    x = as.matrix(spec4[attr(spec4, "ei_x")]); x = x / rowSums(x)
    total = attr(spec4, "ei_n")

    nu2_short = mean((rr4$weights %*% q)^2)
    z_aug = cbind(z_riesz, a)
    fit_aug = ei_riesz_impl(x, z_aug, total, rep(1, nrow(spec4)), rr4$penalty)
    nu2_full  = mean((fit_aug$alpha %*% q)^2)
    actual = 1 - nu2_short / nu2_full
    expect_lt(abs(actual - target), tol)
})

test_that("c_outcome calibration works with an outcome contrast", {
    set.seed(1968)
    target = 0.15
    cv = c(1, -1, 0, 0)  # pres_dem_hum - pres_rep_nix
    contr = list(outcome = cv)
    out_nm = check_contrast(contr, attr(spec4, "ei_x"), attr(spec4, "ei_y"))$y_nm[1]

    ss = ei_semisynthetic(m4, rr4, spec4, predictor = "vap_white", outcome = out_nm,
                          est = ei_est(m4, rr4, spec4, contrast = contr),
                          b_cov = bc4,
                          c_outcome = target, c_predictor = 0.1, contrast = contr)
    a = ss$.confounder

    y_nm = attr(spec4, "ei_y")
    y_contrast = (as.matrix(ss[y_nm]) - m4$fitted) %*% cv
    r2 = cor(as.vector(y_contrast), a)^2
    expect_lt(abs(r2 - target), tol)
})

# --------------------------------------------------------------------------
# predictor/outcome argument validation
# --------------------------------------------------------------------------

test_that("invalid predictor name is rejected", {
    expect_error(
        ei_semisynthetic(m1, rr1, spec1, predictor = "nonexistent",
                         est = est1, b_cov = bc1),
        "predictor"
    )
})

test_that("invalid outcome name is rejected", {
    expect_error(
        ei_semisynthetic(m1, rr1, spec1, outcome = "nonexistent",
                         est = est1, b_cov = bc1),
        "outcome"
    )
})

# --------------------------------------------------------------------------
# Attribute correctness
# --------------------------------------------------------------------------

test_that("b attribute has correct dimensions", {
    set.seed(1968)
    ss = ei_semisynthetic(m4, rr4, spec4, est = est4, b_cov = bc4,
                          c_outcome = 0.15, c_predictor = 0.1)
    b = attr(ss, "b")
    n_x = length(attr(spec4, "ei_x"))
    n_y = length(attr(spec4, "ei_y"))
    expect_equal(dim(b), c(nrow(spec4), n_x * n_y))
    expect_true(all(is.finite(b)))
})

test_that("eta attribute has correct dimensions", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.2, c_predictor = 0.1)
    eta = attr(ss, "eta")
    n_x = length(attr(spec1, "ei_x"))
    n_y = length(attr(spec1, "ei_y"))
    expect_equal(dim(eta), c(nrow(spec1), n_x * n_y))
})

test_that("stored c_outcome, c_predictor, confounding match inputs", {
    set.seed(1968)
    ss = ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                          c_outcome = 0.3, c_predictor = 0.15, confounding = 0.8)
    expect_equal(attr(ss, "c_outcome"), 0.3)
    expect_equal(attr(ss, "c_predictor"), 0.15)
    expect_equal(attr(ss, "confounding"), 0.8)
})

# --------------------------------------------------------------------------
# Input validation 
# --------------------------------------------------------------------------

test_that("non-ei_spec input is rejected", {
    expect_error(ei_semisynthetic(spec = data.frame(x = 1)), "ei_spec")
})

test_that("out-of-range c_outcome/c_predictor/confounding are rejected", {
    expect_error(
        ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                         c_outcome = -0.1, c_predictor = 0.1),
        "between 0 and 1"
    )
    expect_error(
        ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                         c_outcome = 0.1, c_predictor = 1.5),
        "between 0 and 1"
    )
    expect_error(
        ei_semisynthetic(m1, rr1, spec1, est = est1, b_cov = bc1,
                         c_outcome = 0.1, c_predictor = 0.1, confounding = 2),
        "between 0 and 1"
    )
})

test_that("wrong regr type is rejected", {
    expect_error(
        ei_semisynthetic(list(x = 1), rr1, spec1, est = est1, b_cov = bc1),
        "ei_ridge"
    )
})
