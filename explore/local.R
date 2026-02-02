devtools::load_all(".")
library(tidyverse)

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
rr = ei_riesz(spec, penalty = m$penalty)

left_join(ei_est(m, rr, spec), ei_bounds(spec, bounds=0:1, global=TRUE)) |>
    mutate(
        within = min <= estimate & estimate <= max,
        se_max = (max - min) / sqrt(12),
        se_within = std.error <= se_max
    )

(colMeans(est_check_regr(m, spec, nrow(spec), NULL, 4, T)$vcov) |>
    matrix(3, 3) * m$sigma2[1]) |>
    diag() |>
    sqrt()

ei_est(m, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)
ei_est(m, rr, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)


b_cov = ei_local_cov(m, spec)
cov_x = cov2cor(c(20, 1, 2) %o% c(20, 1, 2) + diag(3)*1e-2)
cov_y = cov(m$y) + diag(4)*1e-4
b_cov2 = cov_y %x% cov_x

ei_est_local(m, spec, b_cov = b_cov2, bounds=c(0, 1), sum_one = TRUE) |>
    # summarize(est = weighted.mean(estimate, wt), .by = c(predictor, outcome)) |>
    # arrange(predictor)
    ggplot(aes(estimate)) +
    facet_grid(outcome ~ predictor, scales="free") +
    geom_histogram(bins=20)


ei_est_local(m, spec, b_cov = b_cov2, bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = T) |>
    dplyr::filter(.row == 592) |>
    # dplyr::filter(.row == 1) |>
    ggplot(aes(estimate, paste0(outcome, " | ", predictor))) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size=0.25)

ei_est_local(m, spec, b_cov = b_cov2, contrast = list(outcome=c(1, -1, 0, 0)), bounds = c(0, 1), sum_one = TRUE, conf_level = 0.95, regr_var = T) |>
    dplyr::filter(.row == 592) |>
        print() |>
    ggplot(aes(estimate, predictor)) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high))


H = diag(n_y) %x% local_basis(x[1, ])
R = chol(b_cov)

eig = eigen(oblique_proj(H, R))
ev = ifelse(eig$values > 1e-8, 1/eig$values, 0)
eig$vectors %*% diag(ev) %*% t(eig$vectors)

oblique_proj(H, R)

sqrt(diag(oblique_proj(H, R)))

Sigma = matrix(c(4, 1, 1, 1), 2, 2)
n = 1e4
z = matrix(rnorm(2*n), nrow=n) %*% chol(Sigma)
plot(z, cex=0.25, pch=16)

alphas = seq(0.05, 1, 0.05)^2
plot(alphas, sapply(alphas, \(alpha) {
#   1 - mean(rowSums((z %*% solve(Sigma)) * z) <= (2 - 1)/alpha)
    1 - mean(t(abs(z)) <= 2/3 * sqrt(diag(Sigma)) / sqrt(alpha))
})); abline(0, 1, col='red')

sqrt(diag(Sigma))



