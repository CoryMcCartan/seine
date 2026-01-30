devtools::load_all(".")
library(tidyverse)

data(elec_1968)
elec_1968 = elec_1968 |>
    mutate(
        vap_nonwhite = 1 - vap_white,
        pres_abs = pmax(1e-6, pres_abs),
        z = bases::b_bart(pop_urban, pop_rural, educ_elem, educ_hsch, educ_coll, farm,
            inc_00_03k, inc_03_08k, inc_08_25k, inc_25_99k)
    ) |>
    ei_proportions(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs, clamp = 1e-12) |>
    select(-.total)

spec = ei_spec(elec_1968, c(vap_white, vap_black, vap_other), c(pres_dem_hum, pres_rep_nix, pres_ind_wal, pres_abs),
               total = pres_total, covariates = c(state, z))

m = ei_ridge(spec)#, bounds=0:1, sum_one=T)
rr = ei_riesz(spec, penalty = m$penalty)

ei_est(m, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)
ei_est(m, rr, data = spec) |>
    summarize(err = sum(estimate) - 1, .by = predictor)


r_cov = ei_resid_cov(m, spec)
cov_x = cov2cor(c(20, 1, 1) %o% c(20, 1, 1) + diag(3)*1e-2)
# cov_y = cov(m$y) + diag(4)*1e-4
cov_y = diag(m$sigma2)
r_cov2 = cov_y %x% cov_x

ei_est_local(m, spec, r_cov = r_cov2, bounds=c(0, 1), sum_one = T) |>
    # summarize(est = weighted.mean(estimate, wt), .by = c(predictor, outcome)) |>
    # arrange(predictor)
    ggplot(aes(estimate)) +
    facet_grid(outcome ~ predictor, scales="free") +
    geom_histogram(bins=20)

H = diag(n_y) %x% local_basis(x[1, ])
