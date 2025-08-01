# EI model fitting function
ei_tmvn_impl <- function(x, y, z, weights, bounds, penalty) {
    if (ncol(x) == 2 && ncol(y) == 1) {
        ei_tmvn_2x2_impl(x, c(y), z, weights, bounds, draws_local=0L)
    } else {
        cli_abort("{.fn ei_tmvn} not implemented yet beyond 2x2 case or for covariates.",
                  call=parent.frame())
    }
}

tmvn_pr = list(
    rho_shp = 4,
    sigma_shp = 2,
    sigma_loc = 0.5,
    eta_sd = 2
)

make_nlpost = function(x, y, weights) {
    function(upars) {
        pars = constr_pars(upars)

        out = -(
            R_llik(pars$eta, pars_to_L(pars), y, x, weights, 1e-8) +
                # LKJ(4) prior on rho
                dbeta((pars$rho + 1)/2, tmvn_pr$rho_shp, tmvn_pr$rho_shp, log=TRUE) +
                # prior on sigmas
                sum(dgamma(pars$sigma, tmvn_pr$sigma_shp,
                    tmvn_pr$sigma_shp / tmvn_pr$sigma_loc, log=TRUE)) +
                # prior on etas
                sum(dnorm(pars$eta, 0.5, tmvn_pr$eta_sd, log=TRUE))
        )
        out
    }
}
make_nlpost_gr = function(x, y) {
    function(upars) {
        pars = constr_pars(upars)
        jac_c = constr_jac(upars)

        out = -(
            R_llik(pars$eta, pars_to_L(pars), y, x, weights, 1e-8) +
                # LKJ(4) prior on rho
                ld_dbeta((pars$rho + 1) / 2, tmvn_pr$rho_shp, tmvn_pr$rho_shp) +
                # prior on sigmas
                sum(ld_dgamma(pars$sigma, tmvn_pr$sigma_shp,
                    tmvn_pr$sigma_shp / tmvn_pr$sigma_loc)) +
                # prior on etas
                sum(ld_dnorm(pars$eta, 0.5, tmvn_pr$eta_sd))
        )
        out
    }
}

ei_tmvn_2x2_impl = function(x, y, z, weights, bounds, penalty, draws_local=1000L) {
    if (ncol(z) > 0) {
        cli_abort("Covariates not yet supported.")
    }

    n = length(y)
    upars = unconstr_pars(list(
        eta = c(0.5, 0.5),
        sigma = c(0.1, 0.1),
        rho = 0.5,
        z = rep(0, ncol(z))
    ))

    nlpost = make_nlpost(x, y, weights)

    res = optim(
        upars,
        nlpost,
        method = "L-BFGS-B",
        hessian = TRUE,
        lower = c(-10, -10, -Inf, -Inf, -4), # bounds help avoid pathologies
        upper = c(10, 10, 5, 5, 4),
        control = list(factr = 1e8)
    )
    if (res$convergence > 0) {
        cli_abort(c("Optimization did not converge.", "i"=res$message))
    }
    pars = constr_pars(res$par)

    n_sim = 1000L
    chol_hess = chol(res$hessian)
    zz = matrix(rnorm(n_sim)*5, nrow=5)
    pars_dr = t(backsolve(chol_hess, zz) + res$par)
    # TODO maybe importance sampling, though looks bad
    # lpost_dr = apply(pars_dr, 1, function(x) tryCatch(-nlpost(x), error=\(e) -Inf))
    # lq_dr = sum(log(diag(chol_hess))) - 0.5 * colSums(zz^2)
    pars_dr[, 3:4] = exp(pars_dr[, 3:4])
    pars_dr[, 5] = tanh(pars_dr[, 5])
    vcov = cov(pars_dr)

    if (draws_local > 0) {
        L = pars_to_L(pars)
        b = R_draw_local(draws_local, pars$eta, L, y, x, warmup, 1e-6)
        b = array(b, c(draws_local, n, ncol(x)))
    }

    out = list(
        est = pars,
        est_draws = pars_dr,
        se = sqrt(diag(vcov)),
        vcov = vcov
    )
    if (draws_local > 0) {
        out$b = b
        out$b_global = colSums(colMeans(b) * x * weights) / colSums(x * weights)
    }
    out
}

# parameter transformations
constr_pars = function(upars) {
    idx_z = 5L + seq_len(length(upars) - 5)
    list(
        eta = upars[1:2],
        sigma = exp(upars[3:4]),
        rho = tanh(upars[5]),
        beta = upars[idx_z]
    )
}
constr_jac = function(pars) {
    nz = length(pars$z)
    jac = matrix(0, nrow=6 + nz, ncol=5 + nz)
    idx_z = 5L + seq_len(nz)
    rho_c = 1 - pars$rho^2

    jac[1:2, 1:2] = diag(2) # eta
    jac[3, 3] = pars$sigma[1] # row 3 is L[1, 1] = sigma[1]
    # row 4 is L[2, 1] = sigma[2] * rho
    jac[4, 4] = pars$sigma[2] * pars$rho
    jac[4, 5] = pars$sigma[2] * (rho_c)
    # row 5 is L[2, 2] = sigma[2] * sqrt(1 - rho^2)
    jac[5, 4] = pars$sigma[2] * sqrt(rho_c)
    jac[5, 5] = -pars$sigma[2] * pars$rho * sqrt(rho_c)
    jac[6, 5] = rho_c # rho
    jac[idx_z, idx_z] = diag(nz)

    jac
}

unconstr_pars = function(pars) {
    c(pars$eta, log(pars$sigma), atanh(pars$rho), pars$z)
}
pars_to_L = function(pars) {
    matrix(c(pars$sigma[1], pars$sigma[2]*pars$rho,
                0, pars$sigma[2]*sqrt(1 - pars$rho^2)), 2, 2)
}

# derivatives
ld_dbeta <- function(x, shape1, shape2) {
    (shape1 - 1) * log(x) + (shape2 - 1) * log1p(-x)
}
ld_dgamma <- function(x, shape, rate) {
    (shape - 1) * log(x) - x*rate
}
ld_dnorm <- function(x, mean, sd) {
    -(x - mean) / sd
}