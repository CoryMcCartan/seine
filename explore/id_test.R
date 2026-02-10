devtools::load_all(".")
library(tidyverse)

# try on simulated data
n = 1000
n_perm = 100

run_test <- function(n = 1000, n_x = 3, p = 2, r2 = 0.5, iter = 1000) {
    spec0 = ei_synthetic(
        n,
        p,
        n_x = n_x,
        z = c(1, rep(0, p - 1)),
        r2_xz = r2,
        r2_bz = r2,
        b_cov = 0.0004 * (1 + diag(n_x))
    )
    spec = ei_spec(
        spec0,
        predictors = starts_with("x"),
        outcome = starts_with("y"),
        total = attr(spec0, "ei_n"),
        covariates = starts_with("z"),
        # covariates = "z1",
        preproc = function(z) {
            if (ncol(z) == 0) {
                matrix(nrow=nrow(z), ncol=0)
            } else {
                bases::b_tpsob(z, p = 50)
            }
        }
    )

    ei_test_car(spec, iter = iter, use_chisq = F)
}

res = map(1:400, ~ run_test(n = 500, iter = 100), .progress = TRUE) |>
    bind_rows()

hist(res$df, breaks=50)
hist(res$p.value, breaks=50)
hist(res$p.value[res$df > 0], breaks=50)
mean(res$p.value <= 0.05)
hist(res$W, breaks=50)
hist(pchisq(res$W, res$df, lower.tail=F), breaks=50)
plot(pchisq(res$W, res$df, lower.tail=F), res$p.value)


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
