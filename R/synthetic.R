#' Generate synthetic ecological data
#'
#' Samples data from a following truncated Normal ecological model. The data can
#' be generated completely at random, or can be generated conditional on
#' provided predictors `x` and/or covariates `z`.
#'
#' This function samples data from the following truncated Normal ecological
#' model: \deqn{
#'      \begin{pmatrix}x_i\\ z_i\end{pmatrix} \stackrel{\text{iid}}{\sim}
#'          \mathcal{N}_{[0,1]^{n_x} \times \mathbb{R}^p}\left(
#'          \begin{pmatrix}\mu_x\\ 0\end{pmatrix},
#'          \begin{pmatrix}\Sigma_x & \Gamma \\ \Gamma & T\end{pmatrix}\right)
#' } \deqn{
#'      \eta = z_i^\top \Lambda + \mathtt{b_{loc}}
#' } \deqn{
#'      b_i \stackrel{\text{iid}}{\sim} \mathcal{N}_{[0, 1]^{n_x}}(\eta, \mathtt{B_{cov}})
#' } \deqn{
#'      y_i = b_i^\top x_i,
#' } where \eqn{\mu_x} and \eqn{\Sigma_x} are the mean and covariance of the
#' Normal approximation to a Dirichlet distribution with parameters supplied by
#' the `x` argument below, and \eqn{\Gamma}, \eqn{T}, and \eqn{\Gamma} are
#' matrices sampled to have certain properties, as described below.
#' The subscripts on \eqn{\mathcal{N}} indicate truncation; i.e., both the
#' predictors `x` and the unit-level parameters `b` are truncated to the
#' *n_x*-dimensional hypercube.
#'
#' The matrix \eqn{T} is a symmetric Toeplitz matrix with diagonals provided by
#' the `z` argument. Generally, a decreasing set of nonnegative values will be
#' sufficient for a positive definite \eqn{T}.
#'
#' The matrices \eqn{\Gamma} and \eqn{\Lambda} are initially filled with
#' independent samples from a standard Normal distribution. \eqn{\Gamma} is then
#' projected so that its rows sum to zero, preserving the sum-to-1 requirement
#' on `x`, and so that its columns are scaled to produce the correct \eqn{R^2}
#' value matching `r2_xz`. The matrix \eqn{\Lambda} is likewise scaled to
#' produce the correct \eqn{R^2} value matching `r2_bz`. Due to the truncation
#' in the sampling of `x` and `b`, the in-sample \eqn{R^2} values will generally
#' be slightly smaller than the provided arguments.
#'
#' Aspects of the model can be replaced with data provided to the function.
#' If `x` or `z` is provided as a matrix or data frame, then the other value is
#' sampled from its marginal distribution. If both are provided, then the first
#' row of the model is skipped.
#'
#'
#' @param n The number of rows (geographies) to generate. Defaults to the number
#'   of rows in `x` or `z`, if they are a matrix or data frame.
#' @param p The number of covariates. Defaults to the number of columns in `z`,
#'   if it is a matrix or data frame, or the length of `z`, if it is a vector
#'   of singular values.
#' @param n_x The number of predictor variables. Defaults to the number of
#'   columns in `x`, if it is a matrix or data frame, or the length of `x`, if
#'   it is a vector of mean parameters for the softmax-transformed Normal
#'   distribution.
#' @param x Either a matrix or data frame containing the predictor percentages
#'   in each row, or a vector containing Dirichlet parameters to use in sampling
#'   predictor percentages.
#' @param z A matrix or data frame containing geography-level covariates, or a
#'   vector of values to form a Toeplitz covariance matrix for the random
#'   covariates.
#' @param r2_xz The approximate \eqn{R^2} of the covariates `z` and predictors
#'   `x`. See the model specification for details. If either `r2_xz` or `r2_bz`
#'   are zero, then there is no confounding, and an unadjusted Goodman
#'   regression will estimate the global parameters correctly.
#' @param r2_bz The approximate \eqn{R^2} of the covariates `z` and unit-level
#'   parameters `b`. See the model specification for details. If either `r2_xz`
#'   or `r2_bz` are zero, then there is no confounding, and an unadjusted
#'   Goodman regression will estimate the global parameters correctly.
#' @param b_loc The center of the distribution of geography-level parameters.
#'   Defaults to a linearly spaced sequence across groups from 0.5 to 0.9.
#'   Because of the truncation, this will not exactly be the mean of the
#'   geography-level parameters.
#' @param b_cov The residual covariance matrix for geography-level parameters.
#'   Defaults to `0.02 * (1 + diag(n_x))`.
#'
#' @returns An `ei_spec` object with additional attributes:
#'   - `b_loc` and `b_cov`
#'   - `Lambda` with the coefficients of `z`
#'   - `eta`, the linear predictor for `b`
#'   -  `est_true`, the mean of the geography-level parameters, formatted
#'   similarly to the output from  [ei_est()]
#'   - `r2_xz_act` and `r2_bz_act`, containing the actual (sample) \eqn{R^2}
#'   values for `x` and `z`, and `b` and `z`, respectively.
#'
#' @examples
#' ei_synthetic(n = 10)
#'
#' ei_synthetic(n = 10, p = 2, n_x = 3)
#'
#' # Manual hyperparameters: x2 dominant and z1, z2 very correlated
#' ei_synthetic(n = 10, x = c(1, 95, 4), z = c(10, 9.999))
#'
#' # Condition on provided x but not z
#' data(elec_1968)
#' ei_synthetic(
#'     x = cbind(elec_1968$pop_white, 1 - elec_1968$pop_white),
#'     p = 5,
#'     b_loc = c(0.3, 0.9),
#'     b_cov = matrix(c(0.02, 0.016, 0.016, 0.2), nrow=2)
#' )
#'
#' @export
ei_synthetic = function(n, p = 0, n_x = 2, x = n_x:1, z = 0.25 * exp(-(seq_len(p) - 1)/2),
                        r2_xz = 0.5, r2_bz = 0.5, b_loc = NULL, b_cov = NULL) {
    if (is.matrix(x) || is.data.frame(x)) {
        need_x = FALSE
        x = as.matrix(x)
        n = nrow(x)
        n_x = ncol(x)
    } else {
        need_x = TRUE
        n_x = length(x)
    }
    if (is.matrix(z) || is.data.frame(z)) {
        need_z = FALSE
        z = as.matrix(x)
        n = nrow(z)
        p = ncol(z)
    } else {
        if (is.null(z)) {
            z = exp(-(seq_len(p) - 1)/2)
        }
        need_z = TRUE
        p = length(z)
    }
    if (missing(n)) {
        cli_abort("If {.arg x} and {.arg z} are not a matrix or data frame,
                  you must specify {.arg n}.")
    }

    dd = p*need_z + n_x*need_x
    iz = seq_len(p*need_z)
    ix = p*need_z + seq_len(n_x*need_x)
    mu = numeric(dd)
    L = matrix(0, nrow = dd, ncol = dd)

    if (need_x) {
        if (n_x < 2) {
            cli_abort("{.arg n_x} must be at least 2.")
        }
        # normal approx .to Dirichlet
        sx = sum(x)
        mu[ix] = x / sx
        L[ix, ix] = chol(((1 + 1e-15)*diag(mu[ix]) - tcrossprod(mu[ix])) / (sx + 1))
        vxb = diag(crossprod(L[ix, ix]))
    }
    if (need_z && p > 0) {
        mu[iz] = 0.5
        L[iz, iz] = chol(toeplitz(z))
    }
    if (need_z && need_x && p > 0) {
        # random coefficients, projected to constraints:
        # - rows should sum to 0
        # - variance injected by Z should be prop. to variance of X
        L[iz, ix] = rnorm(p * n_x * need_x * need_z)
        L[iz, ix] = scale_cols(L[iz, ix, drop=FALSE], sqrt(vxb))
        res = optim(L[iz, ix], function(x) {
            Lxz = matrix(x, nrow=p, ncol=n_x)
            2*sum(abs(rowSums(Lxz))) + sum(abs(colSums(Lxz^2) - vxb)) + 1e-9*sum((Lxz - L[iz, ix])^2)
        }, method="BFGS", control=list(abstol=1e-4))
        L[iz, ix] = res$par

        # make R^2 work
        if (length(r2_xz) != 1 || r2_xz < 0 || r2_xz > 1) {
            cli_abort("{.arg r2_xz} must be a single value between 0 and 1.")
        }
        L[ix, ix] = scale_cols(L[ix, ix], sqrt(1 - r2_xz))
        L[iz, ix] = scale_cols(L[iz, ix, drop=FALSE], sqrt(r2_xz))
    }

    warmup = 10L
    thin = 5L
    L[, iz] = 1e-1 * L[, iz]
    R_sync_rng()
    xz = R_ess_tmvn(warmup + thin*n, mu, t(L), init=mu)
    if (p > 0) {
        if (need_z) {
            z = 1e1 * (xz[seq(warmup + 1, nrow(xz), by=thin), iz, drop=FALSE] - 0.5)
        }
        z = shift_cols(z, colMeans(z))
    }
    if (need_x) {
        x = xz[seq(warmup + 1, nrow(xz), by=thin), ix, drop=FALSE]
        x = x / rowSums(x)
    }

    if (is.null(b_loc)) {
        b_loc = seq(0.5, 0.9, length.out=n_x)
    }
    if (is.null(b_cov)) {
        b_cov = 0.02 * (1 + diag(n_x))
    }

    if (p > 0) {
        Lambda = matrix(rnorm(p * n_x), nrow=p, ncol=n_x)
        # make R^2 work
        if (length(r2_bz) != 1 || r2_bz < 0 || r2_bz > 1) {
            cli_abort("{.arg r2_xz} must be a single value between 0 and 1.")
        }
        r2_sc = sqrt(diag(b_cov) * r2_bz / (1 - r2_bz))
        Lambda = scale_cols(Lambda, r2_sc / sqrt(colMeans((z %*% Lambda)^2)))
        eta = shift_cols(z %*% Lambda, -b_loc)
    } else {
        Lambda = NULL
        eta = rep(1, n) %o% b_loc
    }

    b = matrix(nrow=n, ncol=n_x)
    L = t(chol(b_cov))
    for (i in seq_len(n)) {
        # b[i, ] = R_ess_tmvn(warmup, eta[i, ], L, init=eta[i, ])[warmup, ]
        b[i, ] = R_ess_tmvn(warmup, eta[i, ], L, init=rep(0.5, n_x))[warmup, ]
        eta[i, ] = R_ep_moments(eta[i, ], L, numeric(0), 0, 1e-4)[[2]]
    }

    colnames(x) = paste0("x", seq_len(n_x))
    if (p > 0) colnames(z) = paste0("z", seq_len(p))
    d = as.data.frame(cbind(y=rowSums(b * x), x, z))
    true = new_tibble(list(
        predictor = colnames(x),
        outcome = rep("y", n_x),
        true = colSums(eta * x) / colSums(x)
    ))

    if (p > 0) {
        r2_xz_act = unname(sapply(summary(lm(x ~ z)), \(s) s$r.squared))
        r2_bz_act = unname(sapply(summary(lm(b ~ z)), \(s) s$r.squared))
    } else {
        r2_xz_act = rep(0, 3)
        r2_bz_act = rep(0, 3)
    }

    new_tibble(
        d,
        ei_x = colnames(x),
        ei_y = "y",
        ei_z = colnames(z),
        ei_n = rep(1, n),
        b = b,
        b_loc = b_loc,
        b_cov = b_cov,
        Lambda = Lambda,
        eta = eta,
        est_true = true,
        r2_xz_act = r2_xz_act,
        r2_bz_act = r2_bz_act,
        class = "ei_spec"
    )
}
