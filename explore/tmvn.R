ei_tmvn(pres_ind_wal ~ vap_black, elec_1968)

d = elec_1968
d$vap_black[d$pres_ind_wal < 0.2] = 1
ei_tmvn(pres_ind_wal ~ vap_black, d)

pars_to_L = function(pars) {
    matrix(c(pars$sigma[1], pars$sigma[2]*pars$rho,
                0, pars$sigma[2]*sqrt(1 - pars$rho^2)), 2, 2)
}

# save(pars, x, y, file="tests/testthat/llik.rda")
load("tests/testthat/llik.rda")
n = length(y)
L = pars_to_L(pars)
eps_loc = y - rowSums(x * eta)
var_loc = diag(tcrossprod(x %*% L))
wt = rep(1, n)

i = 1
log_Z = R_ep_moments(eta[i, , drop=FALSE], L, numeric(0), 0, 1e-7)[[1]]

bad_i = integer()
for (i in 1:n) {
    val = R_llik(pars$eta, L, y[i], t(x[i, ]), wt[i], 1e-6)
    if (!is.finite(val)) {
        bad_i = c(bad_i, i)
    }
}

bad_i
y[bad_i]
x[bad_i, ]

i = bad_i[1]
j = if (x[i, 1] > 0) 1 else 2
L_eta_proj = r_proj_mvn(pars$eta, L, x[i, ], eps_loc[i])
lb = (y[i] - x[i, 3 - j]) / x[i, j]
ub = y[i] / x[i, j]
lsd = if (j == 1) L_eta_proj[j, 1] else sqrt(sum(L_eta_proj[, 1]^2))

diff(pnorm(c(lb, ub), L_eta_proj[j, 2], lsd))
R_utn_logZ(L_eta_proj[j, 2], lsd, max(0, lb), min(1, ub))
