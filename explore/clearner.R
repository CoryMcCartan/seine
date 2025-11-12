devtools::load_all(".")
library(tidyverse)

data(elec_1968)
elec_1968$vap_nonwhite = 1 - elec_1968$vap_white

spec = ei_spec(elec_1968, c(vap_white, vap_nonwhite), c(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs),
               total = pres_total, covariates = c(state, pop_urban, farm))
d_pred = tibble(
    predictor = attr(spec, "ei_x"),
    wt = map_dbl(attr(spec, "ei_x"), ~ weighted.mean(spec[[.x]], weights(spec)))
)
d_out = tibble(
    outcome = attr(spec, "ei_y"),
    wt = map_dbl(attr(spec, "ei_y"), ~ weighted.mean(spec[[.x]], weights(spec)))
)
n = nrow(spec)

m = ei_ridge(spec)
rr = ei_riesz(spec, bounds=0:1, penalty=m$penalty)
mc = ei_ridge(spec, bounds=0:1, riesz=rr)

colMeans(weights(rr))
crossprod(resid(m), weights(rr)) / n
crossprod(resid(mc), weights(rr)) / n

plot(fitted(m)[, 1], fitted(mc)[, 1])
plot(fitted(m)[, 2], fitted(mc)[, 2])
plot(fitted(m)[, 3], fitted(mc)[, 3])
plot(fitted(m)[, 4], fitted(mc)[, 4])

ei_est(m, data=spec)
ei_est(m, rr, data=spec)
est = ei_est(mc, rr, data = spec) |>
    dplyr::mutate(estimate = round(estimate, 4)) |>
    print()

rbind(
    true = d_out$wt,
    est = colSums(d_pred$wt * as.matrix(est))
)
rowSums(as.matrix(est))
sum(d_pred$wt * as.matrix(est))
