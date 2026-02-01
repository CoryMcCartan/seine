devtools::load_all(".")
library(tidyverse)

# try on simulated data
n = 1000
n_perm = 100

run_test <- function(n = 1000, n_x = 3, p = 4, n_perm = 1000) {
    spec0 = ei_synthetic(n, p, n_x = n_x, b_cov = 0.0004 * (1 + diag(n_x)))
    spec = ei_spec(spec0, starts_with("x"), starts_with("y"), total = attr(spec0, "ei_n"))
    # attr(spec, "b") |> hist(breaks=50)

    m = ei_ridge(spec)

    XZ = bases::b_tpsob(spec[, c(attr(spec, "ei_x"), attr(spec, "ei_z"))], p = 200)
    # ZZ = bases::b_tpsob(spec[, c(attr(spec, "ei_z"))], p = 200)

    udv = svd(XZ)

    fit0 = ridge_auto(udv, resid(m), rep(1, n), vcov = FALSE)
    r2 = diag(as.matrix(cor(fit0$fitted, resid(m))^2))

    perm = replicate(
        n_perm,
        {
            yy = resid(m)[sample(n), , drop = FALSE]
            fit = ridge_svd(udv, yy, rep(1, n), penalty = fit0$penalty, vcov = FALSE)
            diag(as.matrix(cor(fit$fitted, yy)^2))
        },
        simplify = FALSE
    ) |>
        do.call(cbind, args = _)

    tibble::tibble_row(
        r2 = r2,
        p = (rowSums(perm >= r2) + 1) / (n_perm + 1)
    )
}

res = map(1:200, ~ run_test(n = 2000, n_perm = 50), .progress = TRUE) |>
    bind_rows()

hist(res$p)
mean(res$p <= 0.05)


# try on wallace
data(elec_1968)
elec_1968 = elec_1968 |>
    mutate(
        vap_nonwhite = 1 - vap_white,
        z = bases::b_bart(pop_urban, pop_rural, educ_elem, educ_hsch, educ_coll, farm,
            inc_00_03k, inc_03_08k, inc_08_25k, inc_25_99k)
    )

spec = ei_spec(elec_1968, c(vap_white, vap_black, vap_other), c(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs),
               total = pres_total, covariates = c(state, z))

m = ei_ridge(spec, bounds=F, sum_one=F)

Z = cbind(elec_1968$z, model.matrix(~ state - 1, elec_1968))
XZ = with(elec_1968, bases::b_tpsob(vap_white, vap_black, vap_other, pop_urban, pop_rural, educ_elem, educ_hsch, educ_coll, farm,
            inc_00_03k, inc_03_08k, inc_08_25k, inc_25_99k, p=200)) |>
  cbind(model.matrix(~ state - 1, elec_1968))

udv = svd(XZ)
n = nrow(Z)
fit0 = ridge_auto(udv, resid(m), rep(1, n), vcov = FALSE)
r2 = diag(as.matrix(cor(fit0$fitted, resid(m))^2))

perm = replicate(1000, {
    yy = resid(m)[sample(n), , drop=FALSE]
    fit = ridge_svd(udv, yy, rep(1, n), penalty = fit0$penalty, vcov = FALSE)
    diag(as.matrix(cor(fit$fitted, yy)^2))
})

r2
(rowSums(perm >= r2) + 1) / (1000 + 1)
