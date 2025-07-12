test_that("ei_spec() produces correct specification", {
    data(elec_1968)
    expect_error(ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth), "required")

    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_oth, pres_total)
    expect_no_error(validate_ei_spec(spec))
    expect_equal(attr(spec, "ei_x"), c("white", "black", "other")) # tests str_strip_prefix() too
    expect_equal(attr(spec, "ei_y"), c("dem_hum", "rep_nix", "ind_wal", "abs", "oth"))
    expect_equal(attr(spec, "ei_z"), character(0))
    expect_equal(attr(spec, "ei_n"), elec_1968$pres_total)
})


test_that("EI formulas are parsed correctly", {
    res = ei_forms(y ~ x)
    res2 = ei_forms(y ~ x | 1)
    expect_equal(res, list(outcome=expr(y), predictors=expr(x), covariates=expr(1)))
    expect_equal(res, res2)

    res2 = ei_forms(y1 + y2 ~ x1 + x2 | w1 * w2)
    expect_equal(res2, list(outcome=expr(y1 + y2), predictors=expr(x1 + x2),
                            covariates=expr(w1 * w2)))
})


test_that("Bounds are inferred correctly", {
    expect_error(ei_bounds(c(1, 0)), "less")
    expect_error(ei_bounds(c(1, 1)), "strictly less")
    expect_warning(ei_bounds(c(0, 1), -3:3), "outside")

    outcomes = data.frame(y = seq(0, 1, 0.1))
    expect_equal(ei_bounds(c(0, 1), outcomes), c(0, 1))
    expect_equal(ei_bounds(NULL, outcomes), c(0, 1))
    expect_equal(ei_bounds(NULL, outcomes + 1), c(0, Inf))
    expect_equal(ei_bounds(NULL, outcomes - 1), c(-Inf, Inf))
})
