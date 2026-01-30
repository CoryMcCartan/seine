# Ridge regression implementations ---------------------------------------------

# Does ridge regression with the normal equations
ridge_naive <- function(X, y, weights, penalty=0, vcov=TRUE) {
    Lambda = diag(rep(penalty, ncol(X)))
    coef = solve(crossprod(X, weights * X) + Lambda, crossprod(X, weights * y))
    fitted = X %*% coef
    vcov_u = if (vcov) {
        tcrossprod(solve(crossprod(X, weights * X) + Lambda, t(X)))
    } else {
        NULL
    }
    # Neyman-orthogonal estimate of residual variance
    sigma2 = colMeans((y - fitted)^2 * weights)

    list(
        coef = coef,
        vcov_u = vcov_u,
        fitted = fitted,
        sigma2 = sigma2,
        penalty = penalty
    )
}

# Does ridge regression given SVD of design matrix
# `sqrt_w` (square root of unit weights) must be incorporated into udv,
# i.e., `udv = svd(xz * sqrt(weights))`
ridge_svd <- function(udv, y, sqrt_w, penalty=0, vcov=TRUE) {
    d_pen_c = udv$d / (udv$d^2 + penalty)
    d_uy = d_pen_c * crossprod(udv$u, sqrt_w * y)
    fitted = (udv$u / sqrt_w) %*% (udv$d * d_uy)
    vcov_u = if (vcov) {
        tcrossprod(scale_cols(udv$v, d_pen_c^2), udv$v)
    } else {
        NULL
    }
    # Neyman-orthogonal estimate of residual variance
    sigma2 = colMeans(((y - fitted) * sqrt_w)^2)

    list(
        coef = udv$v %*% d_uy,
        vcov_u = vcov_u,
        fitted = fitted,
        sigma2 = sigma2,
        penalty = penalty
    )
}

# Does ridge regression with LOOCV-minimized penalty given SVD of design matrix
# Any weights must be incorporated into udv
ridge_auto <- function(udv, y, sqrt_w, vcov=TRUE) {
    uy = crossprod(udv$u, y * sqrt_w)
    uow = udv$u / sqrt_w

    loo_mse = function(lpen) {
        d_pen_f = udv$d^2 / (udv$d^2 + 10^lpen)
        hat1m = 1 - rowSums(scale_cols(udv$u, d_pen_f) * udv$u)
        resid = y - uow %*% (d_pen_f * uy)
        mean((sqrt_w * resid / hat1m)^2)
    }

    penalty = 10^(optimize(loo_mse, c(-8, 8), tol=0.01)$minimum)
    d_pen_c = udv$d / (udv$d^2 + penalty)
    fitted = uow %*% (d_pen_c * udv$d * uy)
    vcov_u = if (vcov) {
        tcrossprod(scale_cols(udv$v, d_pen_c^2), udv$v)
    } else {
        NULL
    }
    # Neyman-orthogonal estimate of residual variance
    sigma2 = colMeans(((y - fitted) * sqrt_w)^2)

    list(
        coef = udv$v %*% (d_pen_c * uy),
        vcov_u = vcov_u,
        fitted = fitted,
        sigma2 = sigma2,
        penalty = penalty
    )
}

# Ridge regression with bounds on the fitted values _for each group_
#   and option to enforce sum-to-one on outcomes _within_ each group
# Sets ups as QP minimizing -d_vec %*% beta + (1/2) * t(beta) %*% Dmat %*% beta
# see quadprog::solve.QP documentation
ridge_bounds <- function(xz, z, y, weights, bounds, sum_one=FALSE, penalty=0) {
    n = nrow(xz)
    p = ncol(z)
    n_x = ncol(xz) %/% (1L + p)
    n_y = ncol(y)

    if (n_y == 1 && isTRUE(sum_one)) {
        cli_abort("{.fn ridge_bounds} cannot enforce sum-to-one constraint when there is a single outcome variable.")
    }

    dvecs = as.matrix(crossprod(xz, weights * y))
    Dmat = crossprod(xz, weights * xz) + diag(ncol(xz)) * penalty
    R = backsolve(chol(Dmat), diag(nrow(Dmat)))

    Amat = matrix(0, nrow = p + 1, ncol = n * n_x)
    Aind = matrix(0, nrow = p + 2, ncol = n * n_x)
    int_scale = sum(xz[1, 1:n_x])
    Amat[1, 1:(n*n_x)] = int_scale
    Aind[1, ] = p + 1
    for (j in seq_len(n_x)) {
        idx = (j - 1) * n + seq_len(n)
        Aind[-1, idx] = c(j, n_x + p*(j - 1) + seq_len(p))
        Amat[-1, idx] = t(z)
    }

    enforce = is.finite(bounds)
    if (all(enforce)) {
        Amat = cbind(Amat, -Amat)
        Aind = cbind(Aind, Aind)
        bvec = rep(bounds * c(1, -1), each = n*n_x)
    } else if (enforce[1]) {
        bvec = rep(bounds[1], n*n_x)
    } else if (enforce[2]) {
        Amat = -Amat
        bvec = rep(-bounds[2], n*n_x)
    } else {
        cli_abort("{.fn ridge_bounds} requires at least one finite bound.")
    }

    fit_err = \(e) {
        cli_abort(c(
            "Constrained ridge regression failed with inconsistent constraints.",
            ">" = "Try setting {.arg sum_one=FALSE} or relaxing the bounds."
        ), call = NULL)
    }
    if (isFALSE(sum_one)) {
        coefs = matrix(nrow = nrow(dvecs), ncol = ncol(dvecs))
        for (i in seq_len(n_y)) {
            fit = tryCatch(
                quadprog::solve.QP.compact(R, dvecs[, i], Amat, Aind, bvec, factorized = TRUE),
                error = fit_err
            )
            coefs[, i] = fit$solution
        }
    } else {
        R_y = diag(n_y) %x% R
        # first set is sum-to-one; then copy in
        Aind_y = matrix(0, nrow = (p + 1) * n_y + 1, ncol = n * n_x + ncol(Aind) * n_y)
        Amat_y = matrix(0, nrow = (p + 1) * n_y, ncol = n * n_x + ncol(Amat) * n_y)
        idx_st1 = seq_len(n * n_x)
        Amat_y[, idx_st1] = matrix(1, n_y, n_x) %x% rbind(int_scale, t(z))
        Aind_y[1, idx_st1] = n_y * (p + 1)
        for (i in seq_len(n_y)) {
            # sum-to-1
            idx_row = 1 + (i - 1) * (p + 1) + seq_len(p + 1)
            Aind_y[idx_row, idx_st1] = Aind[-1, idx_st1] + (ncol(xz) * (i - 1))

            # copy in constraints
            idx = n * n_x + (i - 1) * ncol(Aind) + seq_len(ncol(Aind))
            Amat_y[seq_len(p + 1), idx] = Amat
            Aind_y[1 + seq_len(p + 1), idx] = Aind[-1, ] + (ncol(xz) * (i - 1))
            Aind_y[1, idx] = Aind[1, ]
        }
        bvec_y = c(rep(1, n * n_x), rep(1, n_y) %x% bvec)

        do_fit = function(eq) {
            quadprog::solve.QP.compact(R_y, c(dvecs), Amat_y, Aind_y, bvec_y, meq = eq, factorized = TRUE)
        }

        # relax to inequality constraint if sum-to-one fails
        eq_constr = n * n_x
        repeat {
            fit = tryCatch(do_fit(eq_constr), error = \(e) NULL)
            if (!is.null(fit)) break
            if (eq_constr > 0) {
                eq_constr = max(eq_constr - n, 0) # reduce by one group
            } else {
                fit_err()
                break
            }
        }
        if (eq_constr < n * n_x) {
            cli_warn(
                "Relaxing {n * n_x - eq_constr} sum-to-one constraint{?s} to inequality to achieve feasible solution.",
                call = NULL
            )
        }
        coefs = matrix(fit$solution, nrow = nrow(dvecs), ncol = ncol(dvecs))
    }

    fitted = xz %*% coefs
    sigma2 = colMeans((y - fitted)^2 * weights)

    list(
        coef = coefs,
        vcov_u = NULL,
        fitted = fitted,
        sigma2 = sigma2,
        penalty = penalty
    )
}

# Calculate diagonal of hat matrix
ridge_hat_svd <- function(udv, penalty=0) {
    d_pen_f = c(udv$d^2 / (udv$d^2 + penalty))
    rowSums(scale_cols(udv$u, d_pen_f) * udv$u)
}
ridge_hat_naive <- function(
        X,
        weights,
        XXinv = solve(crossprod(X, weights * X) + diag(rep(penalty, ncol(X)))),
        penalty=0
) {
    rowSums(weights * X * (X %*% XXinv))
}

# Riesz regression implementations ---------------------------------------------

# Closed-form Riesz regression via the normal equations
riesz_naive <- function(xz, p, total, weights, group=1, penalty=0) {
    n_x = ncol(xz) %/% (1L + p)
    use = c(group, n_x + p*(group-1) + seq_len(p))
    Dz = crossprod(xz[, use, drop=FALSE], total) / mean(xz[, group] * total)
    Lambda = diag(rep(penalty, ncol(xz)))

    XXinv = solve(crossprod(xz, weights*xz) + Lambda)
    xzAinv = xz %*% XXinv[, use, drop=FALSE]
    alpha = c(xzAinv %*% Dz)

    h1m = 1 - ridge_hat_naive(xz, weights, XXinv, penalty)
    xzi = xz[, use, drop=FALSE] * total / mean(xz[, group] * total)
    loo = c((alpha - rowSums(xzAinv * xzi)) / h1m)

    list(alpha = alpha, loo = loo, nu2 = NA)
}

# Closed-form Riesz regression via SVD
# Any weights must be incorporated into udv
riesz_svd <- function(xz, udv, p, total, weights, sqrt_w, group=1, penalty=0) {
    n_x = ncol(xz) %/% (1L + p)
    use = c(group, n_x + p*(group-1) + seq_len(p))
    Dz = colSums(xz[, use, drop=FALSE] * total) / mean(xz[, group] * total)
    d_pen = c(udv$d / (udv$d^2 + penalty))

    xzAinv = (udv$u / sqrt_w) %*% (d_pen * t(udv$v[use, , drop=FALSE]))
    alpha = xzAinv %*% Dz

    h1m = 1 - ridge_hat_svd(udv, penalty)
    xzi = xz[, use, drop=FALSE] * total / mean(xz[, group] * total)
    loo = rowSums(-xzAinv * shift_cols(xzi, Dz)) / h1m

    # Neyman-orthogonal estimate of criterion fn
    nu2 = sum(crossprod(Dz, udv$v[use, , drop=FALSE])^2 *
                  (2/(udv$d^2 + penalty) - d_pen^2)) / nrow(xz)

    list(alpha = alpha, loo = loo, nu2 = nu2)
}
