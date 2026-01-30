devtools::load_all(".")
library(tidyverse)

data(elec_1968)
elec_1968 = elec_1968 |>
    mutate(vap_nonwhite = 1 - vap_white, pres_abs = pmax(1e-6, pres_abs)) |>
    ei_proportions(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs, clamp = 1e-12) |>
    select(-.total)

spec = ei_spec(elec_1968, c(vap_white, vap_black, vap_other), c(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs),
               total = pres_total, covariates = c(state, pop_urban, pop_rural, educ_elem:educ_coll, farm, inc_00_03k:inc_25_99k))

m = ei_ridge(spec, bounds=0:1, sum_one=F)
m = ei_ridge(spec, bounds=F)
rr = ei_riesz(spec, penalty = m$penalty)

# mean(c(y - rowSums(eta * x)) * weights(rr)[, 2])
# wx = x * weights(spec) / rep(colMeans(x * weights(spec)), each=n)

ei_est(m, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)
ei_est(m, rr, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)

# eif = eta_proj * wx
# est = colMeans(eif)
# vcov = crossprod(shift_cols(eif, est)) / (n - 1)^2
# cbind(estimate=est, std.error=sqrt(diag(vcov)))

ei_est_local(m, spec, conf_level = 0.95, bounds=c(0, 1), sum_one = F) |>
# ei_est_local(m, spec, conf_level = 0.95, bounds=F, sum_one = F) |>
    (\(x) { print(attr(x, "proj_misses")); x })() |>
    # dplyr::filter(estimate < -1e-6 | estimate > 1)
    # print()
    summarize(err = sum(estimate) - 1, .by = c(.row, predictor)) |>
    # arrange(-err)
    pull() |>
    hist()

k = 1
n = nrow(spec)

eta = vapply(rl$preds, function(p) p[, k], numeric(n))

bounds = 0:1
zeros = rep(0, 2)
Amat = cbind(rl$x[i, ], diag(2), -diag(2))
eta_diff = 0 * eta
for (i in seq_len(n)) {
    eps = y[i, ] - rl$yhat[i, ]
    Amat[, 1] = rl$x[i, ]
    b0 = c(eps, bounds[1] - eta[i, ], -bounds[2] + eta[i, ])
    eta_diff[i, ] = quadprog::solve.QP(r_cov[[1]], zeros, Amat, b0, meq=1)$solution
}
eta_proj = eta + eta_diff

range(eta[, 1])
range(eta_proj[, 2])
plot(eta[, 1], eta_proj[, 1]); abline(a=0, b=1, col="red")
plot(eta[, 2], eta_proj[, 2]); abline(a=0, b=1, col="red")

colSums(eta * rl$x * weights(spec)) / colSums(rl$x * weights(spec))
colSums(eta_proj * rl$x * weights(spec)) / colSums(rl$x * weights(spec))



idx = which.max(rl$preds$vap_white < 0)
y = spec$pres_dem_hum[idx]
x = rl$x[idx, ]
b0 = sapply(rl$preds, \(x) x[idx, ])
all.equal(rl$yhat[idx], sum(x * b))

fn_dist <- function(b, x, y) {
    eps = y - sum(x * b)
    1e3*(eps^2 + sum(b[b < 0]^2) + sum((b[b > 1] - 1)^2)) + sum((b - b0)^2)
}
gr_dist <- function(b, x, y) {
    eps = y - sum(x * b)
    1e3*(-2*eps*x + 2*b*(b < 0) + 2*(b-1)*(b > 1)) + 2*(b - b0)
}
(c(fn_dist(b0 + c(1, 0)*1e-6, x, y), fn_dist(b0 + c(0, 1)*1e-6, x, y)) - fn_dist(b0, x, y)) / 1e-6
gr_dist(b0, x, y)

b = optim(b0, fn_dist, gr_dist, x=x, y=y, method="L-BFGS-B")$par |> print()
fn_dist(b, x, y)

crossing(
    b0 = seq(-0.5, 1.5, 0.01),
    b1 = seq(-0.5, 1.5, 0.01),
) |>
    mutate(
        eps = y - b0*x[1] - b1*x[2],
        pen = eps^2 + if_else(b0 < 0, b0^2, 0) + if_else(b1 < 0, b1^2, 0) +
            if_else(b0 > 1, (b0-1)^2, 0) + if_else(b1 > 1, (b1-1)^2, 0)
    ) |>
ggplot(aes(b0, b1, z=sqrt(pen))) +
    geom_contour()
    geom_raster()



