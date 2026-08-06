devtools::load_all(".")
library(dplyr)

data(elec_1968)
elec_1968$log_total = log(elec_1968$pres_total)

spec = ei_spec(
    elec_1968,
    vap_white:vap_other,
    pres_dem_hum:pres_abs,
    total = pres_total,
    covariates = c(state, pop_urban, farm, inc_00_03k, inc_25_99k, educ_elem, educ_coll, log_total),
    preproc = ~ bases::b_tpsob(model.matrix(~ 0 + ., .x), p = 96)
)

m0 = ei_ridge(spec)
m = ei_ridge(spec, bounds = 0:1, sum_one = TRUE)#, weights = weights(spec))
rr = ei_riesz(spec, penalty = m$penalty)#, weights = weights(spec))
mt0 = ei_target(m0, rr, spec, bounds=0:1)
mt = ei_target(m, rr, spec, bounds=0:1)

# Three estimators
est_plugin0 = ei_est(regr = m0, data = spec)
est_dml0    = ei_est(regr = m0, riesz = rr, data = spec)
est_plugin = ei_est(regr = m, data = spec)
est_dml    = ei_est(regr = m, riesz = rr, data = spec)
est_tmle0   = ei_est(regr = mt0, riesz = rr, data = spec) 
est_tmle   = ei_est(regr = mt, riesz = rr, data = spec) 
ei_est_local(mt, spec, b_cov = 0.9, bounds = 0:1) |> 
    filter(predictor == "vap_black", outcome == "pres_ind_wal") |> 
    pull(estimate) |> hist(breaks=64)

# TMLE estimates respect the simplex: in [0,1] and sum to 1
mat_tmle = as.matrix(est_tmle)
cat("\nTMLE in [0, 1]:", all(mat_tmle >= -1e-3 & mat_tmle <= 1 + 1e-3), "\n")
cat("TMLE row sums: ", rowSums(mat_tmle), "\n")

# Compare standard errors
cat("\n=== Std errors (plug-in vs TMLE) ===\n")
print(data.frame(
    predictor = est_plugin$predictor,
    outcome   = est_plugin$outcome,
    se_plugin = round(est_plugin$std.error, 4),
    se_tmle   = round(est_tmle$std.error, 4)
))

# other identities
px = colMeans(spec[attr(spec, "ei_x")] * weights(spec))
py = colMeans(spec[attr(spec, "ei_y")] * weights(spec))
rbind(
    plugin0 = px %*% as.matrix(est_plugin0) - py,
    dml0 = px %*% as.matrix(est_dml0) - py,
    plugin = px %*% as.matrix(est_plugin) - py,
    dml = px %*% as.matrix(est_dml) - py,
    tmle0 = px %*% as.matrix(est_tmle0) - py,
    tmle = px %*% as.matrix(est_tmle) - py
) |> 
    zapsmall()
