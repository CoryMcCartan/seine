warmup = 10
N = 1e4
keep = warmup + 1:N
ALPHA = 1e-6

test_that("Sampling TMVN far from bounds works", {
    Sigma = 0.01 + 0.01*diag(2)
    L = chol_pivot(Sigma)
    z = R_ess_tmvn(warmup + N, c(0.5, 0.5), 0.1*L, init=c(0.1, 0.1))[keep, ]
    expect_gt(ks.test(z[, 1], z[, 2])$p.value, ALPHA)
    expect_gt(shapiro.test(z[1:5000, 1])$p.value, ALPHA*0.1)
})

test_that("Sampling TMVN far outside bounds is approximately uniform", {
    z = R_ess_tmvn(warmup + N, c(0.5, 0.5), 100*diag(2), init=c(0.1, 0.1))[keep, ]
    expect_gt(ks.test(z[, 1], z[, 2])$p.value, ALPHA)
    expect_gt(ks.test(z[, 1], "punif")$p.value, ALPHA)
    expect_true(all(z >= 0))
    expect_true(all(z <= 1))
})

test_that("Symmetric TMVN samples have expected moments", {
    z = R_ess_tmvn(warmup + N, c(0.05, 0.5), 0.1*diag(2), init=c(0.1, 0.1))[keep, ]
    expect_gt(t.test(z[, 2], mu=0.5)$p.value, ALPHA)
})


test_that("Degenerate TMVN samples correctly", {
    z = R_ess_tmvn(warmup + N, c(0.0, 0.5), 0.02*diag(2)-0.01, init=c(0.1, 0.1))[keep, ]
    expect_true(all(z[, 1] >= 0))
    expect_true(all(z[, 2] <= 0.5))
    expect_equal(qr(scale(z))$rank, 1)
})

test_that("Random seed synchronizes", {
    set.seed(16801)
    x1 = R_ess_tmvn(2, c(0.5, 0.5), diag(2), c(0.5, 0.5))
    set.seed(16801)
    x2 = R_ess_tmvn(2, c(0.5, 0.5), diag(2), c(0.5, 0.5))

    set.seed(16801)
    R_sync_rng()
    x3 = R_ess_tmvn(2, c(0.5, 0.5), diag(2), c(0.5, 0.5))
    set.seed(16801)
    R_sync_rng()
    x4 = R_ess_tmvn(2, c(0.5, 0.5), diag(2), c(0.5, 0.5))

    expect_false(all(x1 == x2))
    expect_equal(x3, x4)
})

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

# infer mean for known Sigma
# return z-score
# Note: normal approx doesn't always hold well for truncated likelihood esp when OOB
tmvn_infer = function(mu, Sigma=matrix(0.01*c(2, -1, -1, 2), 2), n=1e4) {
    L = chol_pivot(Sigma)
    X = R_ess_tmvn(n + 10, mu, L, init=c(0.1, 0.1))[10 + 1:n, ]

    nllik = function(theta) {
        -sum(dmvnorm_chol(X, theta, L, log=TRUE)) + n*R_ep_moments(theta, L, numeric(), 0, tol=1e-5)[[1]]
    }
    nllik_gr = function(theta) {
        -colSums(dmvnorm_chol_lgr(X, theta, L)) + n*R_ep_moments(theta, L, numeric(), 0, tol=1e-5)[[4]]
    }

    res = optim(c(0.5, 0.5), nllik, nllik_gr, method="L-BFGS-B", hessian=TRUE)

    (res$par - mu) / sqrt(diag(solve(res$hessian)))
}

test_that("Inference for TMVN is correct", {
    expect_true(all( abs(tmvn_infer(c(0.5, 0.5))) < 5 ))
    expect_true(all( abs(tmvn_infer(c(0, 0))) < 5 ))
    expect_true(all( abs(tmvn_infer(c(-0.5, 0.5))) < 8 ))
    expect_true(all( abs(tmvn_infer(c(-1, -1))) < 10 ))
})
