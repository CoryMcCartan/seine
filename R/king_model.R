# EI model fitting function
ei_model_impl <- function(x, y, z, weights, bounds) {
    # TODO REMOVE
    # x <<- x; y <<- y; z <<- z; assign("weights", weights, envir=rlang::global_env());
    # assign("bounds", bounds, envir=rlang::global_env())

    if (ncol(x) == 2 && ncol(y) == 1 && ncol(z) == 0) {
        ei_model_2x2_impl(x, c(y), rep(1, nrow(x)), weights, bounds, draw_local=TRUE)
    } else {
        cli_abort("{.fn ei_model} not implemented yet beyond 2x2 case or for covariates.",
                  call=parent.frame())
    }
}


ei_model_2x2_impl = function(x, y, z, weights, bounds, draw_local=TRUE) {
    if (!identical(bounds, c(0, 1))) {
        cli_abort("For now, bounds must be [0, 1].", call=parent.frame())
    }

    if (any(x[, 1] == 0)) {
        cli_warn("Homogeneous units found; this special case is not yet implemented.");
    }
    if (any(y == 0 | y == 1)) {
        cli_abort("Unanimous units found; this special case is not yet implemented.");
    }

    n = length(y)
    wt = sum(weights)
    lb = pmax(0, (y - x[, 2])/x[, 1])
    ub = pmin(1, y/x[, 1])
    # init
    unconstrain = function(pars) c(pars$eta, log(pars$sigma), atanh(pars$rho))
    constrain = function(upars) {
        list(eta = upars[1:2], sigma = exp(upars[3:4]), rho = tanh(upars[5]))
    }
    upars = unconstrain(list(eta = c(0.5, 0.5), sigma = c(0.1, 0.1), rho = 0.5))
    pars_to_L = function(pars) {
        matrix(c(pars$sigma[1], pars$sigma[2]*pars$rho,
                 0, pars$sigma[2]*sqrt(1 - pars$rho^2)), 2, 2)
    }

    current = upars
    nlpost = function(upars) {
        current <<- upars
        pars = constrain(upars)

        -(
            R_llik_intonly(pars$eta, pars_to_L(pars), y, x, weights, 1e-6) +
                dbeta((pars$rho + 1)/2, 4, 4, log=TRUE) + # LKJ(4) prior on rho
                sum(dgamma(pars$sigma, 2, 2/0.5, log=TRUE)) + # prior on sigmas
                sum(dnorm(pars$eta, 0.5, 2, log=TRUE)) # prior on etas
        )
    }

    res = optim(upars, nlpost, method="L-BFGS-B", hessian=TRUE)
    if (res$convergence > 0) {
        cli_abort(c("Optimization did not converge.", "i"=res$message))
    }
    n_sim = 1000L
    pars = constrain(res$par)
    chol_hess = chol(res$hessian)
    zz = matrix(rnorm(n_sim)*5, nrow=5)
    pars_dr = t(backsolve(chol_hess, zz) + res$par)
    # TODO maybe importance sampling, though looks bad
    # lpost_dr = apply(pars_dr, 1, function(x) tryCatch(-nlpost(x), error=\(e) -Inf))
    # lq_dr = sum(log(diag(chol_hess))) - 0.5 * colSums(zz^2)
    pars_dr[, 3:4] = exp(pars_dr[, 3:4])
    pars_dr[, 5] = tanh(pars_dr[, 5])
    vcov = cov(pars_dr)

    if (draw_local)  {
        warmup = 5L
        draws = 1000L
        L = pars_to_L(pars)
        b = R_draw_local(draws, pars$eta, L, y, x, warmup, 1e-6)
        b = array(b, c(draws, n, ncol(x)))
    }

    out = list(
        est = pars,
        est_draws = pars_dr,
        se = sqrt(diag(vcov)),
        vcov = vcov
    )
    if (draw_local) {
        out$b = b
        out$b_global = colSums(colMeans(b) * x * weights) / colSums(x * weights)
    }
    out
}
