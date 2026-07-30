devtools::load_all(".")

data(elec_1968)
elec_1968$log_total = sqrt(elec_1968$pres_total)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm, pres_total))

m0 = ei_ridge(spec)#, weights = weights(spec))
m = ei_ridge(spec, bounds = 0:1, sum_one = TRUE)#, weights = weights(spec))
rr = ei_riesz(spec, penalty = m$penalty)#, weights = weights(spec))

# Three estimators
est_plugin0 = ei_est(regr = m0, data = spec)
est_dml0    = ei_est(regr = m0, riesz = rr, data = spec)
est_plugin = ei_est(regr = m, data = spec)
est_dml    = ei_est(regr = m, riesz = rr, data = spec)
est_tmle   = ei_est(regr = m, riesz = rr, data = spec, bounds=0:1)  # TMLE
# ei_est_local(m0, spec, b_cov = 0.9)

cat("=== Plug-in ===\n");  print(as.matrix(est_plugin0))
cat("\n=== DML (one-step) ===\n"); print(as.matrix(est_dml0))
cat("=== Plug-in ===\n");  print(as.matrix(est_plugin))
cat("\n=== DML (one-step) ===\n"); print(as.matrix(est_dml))
cat("\n=== TMLE ===\n");   print(as.matrix(est_tmle))

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
    tmle = px %*% as.matrix(est_tmle) - py
) |> 
    zapsmall()
