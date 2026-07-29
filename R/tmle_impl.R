tmle_target = function(y, rl, rm, rm_cf, bounds, sum_one) {
    fluctuation = tmle_fluctuation(bounds, sum_one)
    rl$yhat = fluctuation$project(rl$yhat)
    rl$preds = lapply(rl$preds, fluctuation$project)

    eps = fluctuation$fit_eps(y, rl$yhat, rm)
    rl$yhat = fluctuation$update(rl$yhat, rm, eps)

    score = crossprod(rm, y - rl$yhat) / nrow(y)
    if (max(abs(score)) > 1e-3) {
        cli_warn("Targeting step may not have converged properly.")
    }

    for (group in names(rl$preds)) {
        rl$preds[[group]] = fluctuation$update(rl$preds[[group]], rm_cf[[group]], eps)
    }
    rl
}

# Construct the counterfactual 'clever covariates' for a fitted Riesz object.
tmle_riesz_cf = function(riesz, data) {
    data = hardhat::forge(data, riesz$blueprint)$predictors
    xcols = riesz$blueprint$ei_x
    idx_x = match(xcols, colnames(data))
    z = as.matrix(data[, -idx_x, drop = FALSE])
    p = ncol(z)
    n_x = length(xcols)

    z = shift_cols(z, riesz$z_shift)
    z = scale_cols(z, riesz$z_scale)
    if (anyNA(z)) {
        cli_abort("Missing values found in covariates.", call=parent.frame())
    }
    z = cbind(riesz$int_scale, z)
    w = riesz$blueprint$ei_wgt / mean(riesz$blueprint$ei_wgt)
    w = w %o% (riesz$int_scale / colMeans(riesz$weights))

    rm_cf = list()
    for (group in seq_along(xcols)) {
        use = c(group, n_x + p*(group - 1) + seq_len(p))
        rm_cf[[xcols[group]]] = w * (z %*% riesz$coef[use, ])
    }
    rm_cf
}



# Select a fluctuation submodel which preserves the requested outcome support.
tmle_fluctuation = function(bounds, sum_one) {
    if (isTRUE(sum_one)) {
        if (!isTRUE(all.equal(bounds, c(0, 1)))) {
            cli_abort(
                "TMLE with {.arg sum_one = TRUE} currently supports only {.code bounds = c(0, 1)}.",
                call = parent.frame()
            )
        }
        list(
            project = tmle_project_simplex,
            fit_eps = function(y, q, H) tmle_fit_eps(y, q, H, tmle_loss_multinomial),
            update = tmle_update_softmax
        )
    } else if (all(is.finite(bounds))) {
        scale = bounds[2] - bounds[1]
        list(
            project = function(q) tmle_project_interval(q, bounds[1], bounds[2]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(
                    (y - bounds[1]) / scale,
                    (q - bounds[1]) / scale,
                    H,
                    tmle_loss_binomial
                )
            },
            update = function(q, H, eps) tmle_update_logistic(q, H, eps, bounds)
        )
    } else if (is.finite(bounds[1])) {
        list(
            project = function(q) tmle_project_lower(q, bounds[1]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(y - bounds[1], q - bounds[1], H, tmle_loss_poisson)
            },
            update = function(q, H, eps) tmle_update_log(q, H, eps, bounds[1])
        )
    } else if (is.finite(bounds[2])) {
        list(
            project = function(q) tmle_project_upper(q, bounds[2]),
            fit_eps = function(y, q, H) {
                tmle_fit_eps(bounds[2] - y, bounds[2] - q, H, tmle_loss_poisson)
            },
            update = function(q, H, eps) tmle_update_log(q, H, eps, bounds[2], upper = TRUE)
        )
    } else{
        cli_abort("TMLE requires at least one finite outcome bound.", call = parent.frame())
    }
}


# Fit a J by K fluctuation parameter matrix after column-normalizing H.
tmle_fit_eps = function(y, q, H, loss) {
    n_x = ncol(H)
    n_y = ncol(y)
    H_scale = pmax(sqrt(colMeans(H^2)), 1e-8)
    H_norm = scale_cols(H, 1 / H_scale)

    nllik = function(eps_vec) {
        eps = matrix(eps_vec, n_x, n_y)
        loss(y, q, H_norm, eps)
    }
    opt = optim(
        par = rep(0, n_x * n_y),
        fn = nllik,
        method = "L-BFGS-B",
        control = list(maxit = 1000, factr = 1e3)
    )
    matrix(opt$par, n_x, n_y) / H_scale
}

tmle_loss_multinomial = function(y, q, H, eps) {
    linpred = log(q) + H %*% eps
    max_lp = linpred[cbind(seq_len(nrow(linpred)), max.col(linpred))]
    linpred = linpred - max_lp
    sum(log(rowSums(exp(linpred))) - rowSums(y * linpred))
}

tmle_loss_binomial = function(y, q, H, eps) {
    linpred = qlogis(q) + H %*% eps
    sum(pmax(linpred, 0) + log1p(exp(-abs(linpred))) - y * linpred)
}

tmle_loss_poisson = function(y, q, H, eps) {
    linpred = log(q) + H %*% eps
    sum(exp(linpred) - y * linpred)
}

tmle_update_softmax = function(q, H, eps) {
    linpred = log(q) + H %*% eps
    max_lp = linpred[cbind(seq_len(nrow(linpred)), max.col(linpred))]
    preds = exp(linpred - max_lp)
    preds / rowSums(preds)
}

tmle_update_logistic = function(q, H, eps, bounds) {
    scale = bounds[2] - bounds[1]
    bounds[1] + scale * plogis(qlogis((q - bounds[1]) / scale) + H %*% eps)
}

tmle_update_log = function(q, H, eps, bound, upper = FALSE) {
    q = if (upper) bound - q else q - bound
    q = q * exp(H %*% eps)
    if (upper) bound - q else bound + q
}

# Duchi et al. (2008) fast simplex projection onto its delta-interior.
tmle_project_simplex = function(x, delta = 1e-8) {
    U = t(apply(x - delta, 1, function(z) sort.int(z, decreasing = TRUE, method = "quick")))
    z = 1 - ncol(x) * delta

    L = diag(ncol(x))
    L[upper.tri(L)] = 1
    cssv = U %*% L
    idxs = rep(seq_len(ncol(x)), each = nrow(x))
    rho = rowSums(U - (cssv - z) / idxs > 0)

    cssv_rho = cssv[cbind(seq_len(nrow(x)), rho)]
    theta = (cssv_rho - z) / rho
    x = x - theta
    x[x < 0] = 0
    x + delta
}

tmle_project_interval = function(q, lower, upper, delta = 1e-8) {
    delta = min(delta, (upper - lower) / 4)
    pmin(pmax(q, lower + delta), upper - delta)
}

tmle_project_lower = function(q, lower, delta = 1e-8) {
    pmax(q, lower + delta)
}

tmle_project_upper = function(q, upper, delta = 1e-8) {
    pmin(q, upper - delta)
}