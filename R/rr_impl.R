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
# Sets ups as QP minimizing -d_vec %*% beta + (1/2) * t(beta) %*% Dmat %*% beta
# see quadprog::solve.QP documentation
ridge_bounds <- function(xz, z, y, weights, bounds, riesz=NULL, penalty=0) {
    n = nrow(xz)
    p = ncol(z)
    n_x = ncol(xz) %/% (1L + p)

    dvecs = as.matrix(crossprod(xz, weights * y))
    Dmat = crossprod(xz, weights * xz) + diag(ncol(xz)) * penalty
    R = backsolve(chol(Dmat), diag(nrow(Dmat)))

    # fitted values
    Amat = matrix(0, nrow = nrow(Dmat), ncol = n * n_x)
    int_scale = sum(xz[1, 1:n_x])
    for (j in seq_len(n_x)) {
        use = n_x + p*(j-1) + seq_len(p)
        idx = (j - 1) * n + seq_len(n)
        Amat[j, idx] = int_scale
        Amat[use, idx] = t(z)
    }

    enforce = is.finite(bounds)
    if (all(enforce)) {
        Amat = cbind(Amat, -Amat)
        bvec = rep(bounds * c(1, -1), each = n*n_x)
    } else if (enforce[1]) {
        bvec = rep(bounds[1], n*n_x)
    } else if (enforce[2]) {
        Amat = -Amat
        bvec = rep(-bounds[2], n*n_x)
    } else {
        cli_abort("{.fn ridge_bounds} requires at least one finite bound.")
    }

    coefs = matrix(nrow = nrow(dvecs), ncol = ncol(dvecs))
    for (i in seq_len(ncol(dvecs))) {
        Amat_i = Amat
        bvec_i = bvec
        meq = 0
        if (!is.null(riesz)) {
            xz_avg = t(crossprod(riesz, xz) / nrow(xz))
            Amat_i = cbind(xz_avg, Amat)
            bvec_i = c(colMeans(riesz * y[, i]), bvec)
            meq = ncol(xz_avg)
        }

        fit = quadprog::solve.QP(R, dvecs[, i], Amat_i, bvec_i, meq=meq, factorized=TRUE)
        coefs[, i] = fit$solution
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

# Riesz regression with bounds on the values for each observation
# Sets ups as QP minimizing -d_vec %*% beta + (1/2) * t(beta) %*% Dmat %*% beta
# see quadprog::solve.QP documentation
riesz_bounds <- function(xz, z, total, weights, bounds, y, group=1, penalty=0) {
    n = nrow(xz)
    p = ncol(z)
    n_x = ncol(xz) %/% (1L + p)
    n_y = ncol(y)
    use = c(group, n_x + p*(group-1) + seq_len(p))

    dvec = numeric(ncol(xz))
    dvec[use] = crossprod(xz[, use, drop=FALSE], total) / mean(xz[, group] * total)
    Dmat = crossprod(xz, weights * xz) + diag(ncol(xz)) * penalty
    R = backsolve(chol(Dmat), diag(nrow(Dmat)))

    int_scale = sum(xz[1, 1:n_x])
    enforce = is.finite(bounds)
    avg = crossprod(xz, weights) / n
    avg_y = crossprod(xz, weights * y) / n
    if (all(enforce)) {
        Amat = cbind(avg, avg_y, avg_y)
        bvec = c(1 / int_scale, rep(bounds[1], n_y), -rep(bounds[2], n_y))
    } else if (enforce[1]) {
        Amat = cbind(avg, avg_y)
        bvec = c(1 / int_scale, rep(bounds[1], n_y))
    } else if (enforce[2]) {
        Amat = cbind(avg, avg_y)
        bvec = c(1 / int_scale, -rep(bounds[2], n_y))
    } else {
        cli_abort("{.fn ridge_bounds} requires at least one finite bound.")
    }

    fit = quadprog::solve.QP(R, dvec, Amat, bvec, factorized=TRUE, meq=1)

    list(alpha = c(xz %*% fit$solution), loo = NA, nu2 = -2 * fit$value)
}