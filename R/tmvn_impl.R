# EI model fitting function
ei_tmvn_impl <- function(x, y, z, weights, bounds, penalty) {
    penalty = max(1e-12, penalty)

    if (ncol(x) == 2 && ncol(y) == 1) {
        z = cbind(rep(1, length(y)), z)
        ei_tmvn_2x2_impl(x, c(y), z, weights, bounds, penalty, draws_local=1000L)
    } else {
        cli_abort("{.fn ei_tmvn} not implemented yet beyond 2x2 case or for covariates.",
                  call=parent.frame())
    }
}

tmvn_pr = list(
    rho_shp = 4,
    sigma_shp = 2,
    sigma_loc = 0.5,
    int_sd = 2
)

make_nlpost = function(x, y, z, weights, beta_sd) {
    function(upars) {
        pars = constr_pars(upars)
        eta = tcrossprod(z, pars$beta)

        out = -(
            R_llik(eta, pars_to_L(pars), y, x, weights, p=ncol(pars$beta), tol=1e-8) +
                # LKJ(4) prior on rho
                dbeta((pars$rho + 1)/2, tmvn_pr$rho_shp, tmvn_pr$rho_shp, log=TRUE) +
                # prior on sigmas
                sum(dgamma(pars$sigma, tmvn_pr$sigma_shp,
                    tmvn_pr$sigma_shp / tmvn_pr$sigma_loc, log=TRUE)) +
                # prior on linear coefs
                sum(dnorm(pars$beta[1, ], 0.5, tmvn_pr$int_sd, log=TRUE)) +
                sum(dnorm(pars$beta[-1, ], 0, beta_sd, log=TRUE))
        )
        # if (!is.finite(out)) browser()
        out
    }
}
make_nlpost_gr = function(x, y, z, weights, beta_sd) {
    cli_abort("Not implemented")
    function(upars) {
        pars = constr_pars(upars)
        jac_c = constr_jac(upars)

        out = -(
            R_llik(eta, pars_to_L(pars), y, x, weights, p=ncol(pars$beta), tol=1e-8) +
                # LKJ(4) prior on rho
                ld_dbeta((pars$rho + 1) / 2, tmvn_pr$rho_shp, tmvn_pr$rho_shp) +
                # prior on sigmas
                sum(ld_dgamma(pars$sigma, tmvn_pr$sigma_shp,
                    tmvn_pr$sigma_shp / tmvn_pr$sigma_loc)) +
                # prior on linear coefs
                sum(ld_dnorm(pars$beta[1, ], 0.5, tmvn_pr$int_sd)) +
                sum(ld_dnorm(pars$beta[-1, ], 0, beta_sd))
        )
        out
    }
}

ei_tmvn_2x2_impl = function(x, y, z, weights, bounds, penalty, draws_local=1000L) {
    n = length(y)
    p = ncol(z)

    upars = unconstr_pars(list(
        sigma = c(0.1, 0.1),
        rho = 0.5,
        beta = cbind(c(0.5, 0.5), matrix(0, nrow=2, ncol=p - 1))
    ))

    nlpost = make_nlpost(x, y, z, weights, sqrt(1 / (2 * penalty)))

    res = optim(
        upars,
        nlpost,
        method = "L-BFGS-B",
        hessian = TRUE,
        lower = c(-Inf, -Inf, -4, rep(-10, p)), # bounds help avoid pathologies
        upper = c(5, 5, 4, rep(10, p)),
        control = list(factr = 1e8)
    )
    if (res$convergence > 0) {
        cli_abort(c("Optimization did not converge.", "i"=res$message))
    }
    pars = constr_pars(res$par)

    n_sim = 1000L
    chol_hess = chol(res$hessian)
    par_unorm = matrix(rnorm(n_sim*ncol(chol_hess)), ncol=n_sim)
    pars_dr = t(backsolve(chol_hess, par_unorm) + res$par)
    # TODO maybe importance sampling, though looks bad
    # lpost_dr = apply(pars_dr, 1, function(x) tryCatch(-nlpost(x), error=\(e) -Inf))
    # lq_dr = sum(log(diag(chol_hess))) - 0.5 * colSums(par_unorm^2)
    pars_dr[, 1:2] = exp(pars_dr[, 1:2])
    pars_dr[, 3] = tanh(pars_dr[, 3])
    vcov = cov(pars_dr)

    if (draws_local > 0) {
        eta = tcrossprod(z, pars$beta)
        L = pars_to_L(pars)
        b = R_draw_local(draws_local, eta, L, y, x, warmup, 1e-6)
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
    list(
        sigma = exp(upars[1:2]),
        rho = tanh(upars[3]),
        beta = matrix(upars[-1:-3], nrow=2)
    )
}
constr_jac = function(pars) {
    n_beta = 2*ncol(pars$beta)
    jac = matrix(0, nrow=4 + n_beta, ncol=3 + n_beta)
    rho_c = 1 - pars$rho^2

    jac[1, 1] = pars$sigma[1] # row 1 is L[1, 1] = sigma[1]
    # row 2 is L[2, 1] = sigma[2] * rho
    jac[2, 2] = pars$sigma[2] * pars$rho
    jac[2, 3] = pars$sigma[2] * (rho_c)
    # row 3 is L[2, 2] = sigma[2] * sqrt(1 - rho^2)
    jac[3, 2] = pars$sigma[2] * sqrt(rho_c)
    jac[3, 3] = -pars$sigma[2] * pars$rho * sqrt(rho_c)
    jac[4, 3] = rho_c # rho
    jac[-1:-4, -1:-3] = diag(n_beta)

    jac
}

unconstr_pars = function(pars) {
    c(log(pars$sigma), atanh(pars$rho), c(pars$beta))
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