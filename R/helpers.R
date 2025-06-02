# These scaling functions have been benchmarked on R 3.4 and found to be the
# fastest way to accomplish this task.

# Subtract from each column a provided factor
shift_cols = function(X, shift) {
    X - rep(shift, each=nrow(X))
}
# Scale each column by a provided factor
scale_cols = function(X, scale) {
    X * rep(scale, each=nrow(X))
}

# Rowwise Kronecker product, with intercept added and weight scaling
# `int_scale` = 0 implies no intercept
# `int_scale` > 0 creates an intercept by additionally prepending `x` itself
#     The scale number is used to scale `x`, to adjust how these terms are
#     affected by shrinkage
row_kronecker = function(x, z, int_scale=1) {
    n_x = ncol(x)
    n_z = ncol(z)

    offset = if (int_scale > 0) n_x else 0
    out = matrix(nrow=nrow(x), ncol=offset + n_x * n_z)
    if (int_scale > 0) {
        out[, 1:n_x] = x * int_scale
    }
    for (j in seq_len(n_x)) {
        idx = offset + (j - 1)*n_z + seq_len(n_z)
        out[, idx] = x[, j] * z
    }

    colnames(out) = c(
        rep_len(colnames(x), offset),
        paste0(rep(colnames(x), each=n_z), ":", rep(colnames(z), n_x), recycle0=TRUE)
    )
    out
}

# Cholesky decomposition
#
# Takes transpose and pivots as needed to account for degenerate covariance matrices.
# Multiply on right by iid Gaussians to simulate.
chol_pivot = function(Sigma) {
    L = chol.default(Sigma, pivot=TRUE)
    t(L[1:attr(L, "rank"), order(attr(L, "pivot")), drop=FALSE])
}



# These functions are currently only used in `explore/`

# MVN density, Cholesky parametrization
# Assumes `x` is a matrix with 1 row per obs.
# Assumes `L` is a lower triangular matrix
# CONSTANTS DROPPED
dmvnorm_chol = function(x, mu, L, log=FALSE) {
    rss = colSums(forwardsolve(L, t(x) - mu)^2)
    out = -sum(log(diag(L))) - 0.5 * rss
    if (!log) exp(out) else out
}
# gradient of log density wrt mu
dmvnorm_chol_lgr = function(x, mu, L) {
    t(chol2inv(t(L)) %*% (t(x) - mu))
}
dmvnorm_chol_rk1 = function(x, mu, L, log=FALSE) {
    rss = colSums(forwardsolve(L, t(x) - mu)^2)
    out = -log(sum(L^2)) - 0.5 * rss
    if (!log) exp(out) else out
}
