test_that("local estimates satisfy constraints", {
    spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
                   total = pres_total,
                   covariates = c(pop_urban, farm, educ_elem, educ_coll, inc_00_03k))

    m = ei_ridge(spec)

    suppressWarnings(
        ests <- ei_est_local(m, spec, bounds=c(0, 1), sum_one = TRUE)
    )
    expect_true(all(ests$estimate > -1e-12))
    expect_true(all(ests$estimate < 1 + 1e-12))
    ea = as.array(ests)
    expect_identical(dim(ea), c(nrow(spec), 3L, 4L))
    expect_lt(mean(abs(rowSums(ea, dims = 2) - 1)), 5e-16)
})
