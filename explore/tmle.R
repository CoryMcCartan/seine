devtools::load_all(".")
library(rlang)
# library(tidyverse)

data(elec_1968)

spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
               total = pres_total, covariates = c(state, pop_urban, farm))

m = ei_ridge(spec, bounds=0:1, sum_one=TRUE)
rr = ei_riesz(spec, penalty = m$penalty)

y = est_check_outcome(m, spec, quo(NULL))
n = nrow(y)
n_y = ncol(y)

subset = eval_tidy(quo(NULL), spec)
subset = as.numeric(check_subset(subset, n))
total = attr(spec, "ei_n")

w = check_make_weights(!!quo(total), spec, n)
w = subset * w / mean(subset * w)

# build predictions and RR
rm = est_check_riesz(rr, spec, w, n, m) # riesz matrix
rl = est_check_regr(m, spec, n, colnames(rm), n_y, vcov = FALSE) # regr list
n_x = ncol(rm)

offsets = log(pmin(pmax(rl$yhat, 1e-8), 1 - 1e-8))
offsets = offsets - offsets[, 1]

eps = matrix(0.01 * rnorm(n_x * n_y), n_x, n_y)

nllik = function(eps) {
    eps = matrix(eps, n_x, n_y)
    linpred = offsets + rm %*% eps
    max_lp = linpred[cbind(seq_len(n), max.col(linpred))]
    linpred = linpred - max_lp
    -sum(y * (linpred - log(rowSums(exp(linpred)))))
}

(nllik(eps + c(1e-6, rep(0, 11))) - nllik(eps)) / 1e-6

e0 = coef(nnet::multinom(y ~ offset(offsets) + rm))
e1 = optim(matrix(rnorm(n_x*n_y, 0, 0.001), n_x, n_y), nllik, method = "L-BFGS-B")$par

# linpred = offsets + rm %*% e1
# max_lp = linpred[cbind(seq_len(n), max.col(linpred))]
# fluct = exp(linpred - max_lp)
fluct = exp(rm %*% e1)

# plot(rl$yhat, fluct)

x_nm = names(rl$preds)
y_nm = colnames(y)
eif0 = matrix(nrow=n, ncol=n_y * n_x) # x varies faster than y
eif1 = eif0
for (i in seq_len(n_x)) {
    plugin = rl$preds[[x_nm[i]]] * rl$x[, i] * w / mean(rl$x[, i] * w)
    wtd = rm[, i] * (y - rl$yhat)
    eif0[, i + (seq_len(n_y) - 1)*n_x] = plugin + wtd

    pred = rl$preds[[x_nm[i]]] * fluct
    pred = pred / rowSums(pred)
    plugin = pred * rl$x[, i] * w / mean(rl$x[, i] * w)
    pred = rl$yhat * fluct
    pred = pred / rowSums(pred)
    wtd = rm[, i] * (y - pred) 
    eif1[, i + (seq_len(n_y) - 1)*n_x] = plugin + wtd
}
matrix(colMeans(eif0), n_x, n_y, dimnames=list(x_nm, y_nm))
matrix(colMeans(eif1), n_x, n_y, dimnames=list(x_nm, y_nm))

# may need to get coefficients for RR