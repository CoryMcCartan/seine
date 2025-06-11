# Ridge regression implementations ---------------------------------------------

# Does ridge regression with the normal equations
ridge_naive <- function(X, y, weights, penalty=0) {
    Lambda = diag(rep(penalty, ncol(X)))
    coef = solve(crossprod(X, weights * X) + Lambda, crossprod(X, weights * y))
    list(coef = coef, fitted = X %*% coef, penalty = penalty)
}

# Does ridge regression given SVD of design matrix
# `sqrt_w` (square root of unit weights) must be incorporated into udv,
# i.e., `udv = svd(xz * sqrt(weights))`
ridge_svd <- function(udv, y, sqrt_w, penalty=0) {
    d_pen_c = udv$d / (udv$d^2 + penalty)
    d_uy = d_pen_c * crossprod(udv$u, sqrt_w * y)
    list(
        coef = udv$v %*% d_uy,
        fitted = (udv$u / sqrt_w) %*% (udv$d * d_uy),
        penalty = penalty
    )
}

# Does ridge regression with LOOCV-minimized penalty given SVD of design matrix
# Any weights must be incorporated into udv
ridge_auto <- function(udv, y, sqrt_w) {
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
    list(
        coef = udv$v %*% (d_pen_c * uy),
        fitted = uow %*% (d_pen_c * udv$d * uy),
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
    Dz = crossprod(xz[, use], total) / mean(xz[, group] * total)
    Lambda = diag(rep(penalty, ncol(xz)))

    XXinv = solve(crossprod(xz, weights*xz) + Lambda)
    xzAinv = xz %*% XXinv[, use]
    alpha = c(xzAinv %*% Dz)
    h1m = 1 - ridge_hat_naive(xz, weights, XXinv, penalty)
    loo = c((alpha - rowSums(xzAinv * xz[, use] * weights)) / h1m)
    list(alpha = alpha, loo = loo)
}

# Closed-form Riesz regression via SVD
# Any weights must be incorporated into udv
riesz_svd <- function(xz, udv, p, total, weights, sqrt_w, group=1, penalty=0) {
    n_x = ncol(xz) %/% (1L + p)
    use = c(group, n_x + p*(group-1) + seq_len(p))
    Dz = colSums(xz[, use] * total) / mean(xz[, group] * total)
    d_pen = c(udv$d / (udv$d^2 + penalty))

    xzAinv = (udv$u / sqrt_w) %*% (d_pen * t(udv$v[use, ]))
    alpha = xzAinv %*% Dz
    h1m = 1 - ridge_hat_svd(udv, penalty)
    loo = rowSums(-xzAinv * shift_cols(xz[, use] * weights, Dz)) / h1m
    list(alpha = alpha, loo = loo)
}
