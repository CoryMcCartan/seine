test_that("ei_proportions() runs correctly", {
    data(elec_1968)
    # Make a data frame with counts
    d_unnorm = with(elec_1968, data.frame(
        vap = vap,
        vap_white = vap * vap_white,
        vap_black = vap * vap_black,
        vap_other = vap * vap_other
    ))

    d1 = ei_proportions(d_unnorm, vap_white:vap_black, .total=vap)
    d2 = ei_proportions(d_unnorm, vap_white:vap_other, .total="vap_copy")
    d3 = ei_proportions(d_unnorm, white=vap_white, black=vap_black,
                        .total=c(total=vap), .other="other")

    expect_equal(d1$vap_white + d1$vap_black + d1$.other, rep(1, nrow(d_unnorm)))
    expect_equal(d1$vap_white, d2$vap_white)
    expect_equal(d1$vap_black, d2$vap_black)
    expect_contains(names(d2), "vap_copy")
    expect_equal(d1$vap, d2$vap_copy)
    expect_setequal(names(d3), c("total", "white", "black", "vap_other", "other"))

    expect_error(ei_proportions(d_unnorm), "at least one")
    expect_error(ei_proportions(d_unnorm, vap_white, .total=c(vap, vap_white)), "single")
    d_unnorm$vap_white[1] = -100
    expect_error(ei_proportions(d_unnorm, vap_white, .total=vap), "negative")
    d_unnorm$vap_white[1] = 1e5
    expect_error(ei_proportions(d_unnorm, vap_white, .total=vap), "more than 1")
})
