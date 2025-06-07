#' Generate Synthetic Ecological Data
#'
#' Samples data from the following truncated Normal ecological model:
#' \deqn{
#'      x_i \stackrel{\text{iid}}{\sim} \mathrm{Dirichlet}(\alpha) \\
#'      \eta = z_i^\top \mathtt{coefs} + \mathtt{b_{loc}} \\
#'      b_i \stackrel{\text{iid}}{\sim} \mathcal{N}_{[0, 1]^d}(\eta, \mathtt{B_{cov}}) \\
#'      y_i = b_i^\top x_i,
#' }
#' where \eqn{\alpha} is a vector of hyperparameters that may be optionally
#' specified as the `x` argument to the function, and \eqn{\mathcal{N}_{[0,
#' 1]^d}} indicates a multivariate Normal distribution truncated to the
#' *d*-dimensional hypercube.
#' Aspects of the model can be replaced with data provided to the function.
#'
#'
#' @param n The number of rows (geographies) to generate. Defaults to the number
#'   of rows in `x`, if it is a matrix or data frame.
#' @param p The number of covariates. Defaults to the number of columns in `z`,
#'   if it is a matrix or data frame, or the length of `z`, if it is a vector
#'   of singular values.
#' @param n_x The number of predictor variables. Defaults to the number of
#'   columns in `x`, if it is a matrix or data frame, or the length of `x`, if
#'   it is a vector of Dirichlet parameters.
#' @param x Either a matrix or data frame containing the predictor percentages
#'   in each row, or a vector containing Dirichlet parameters to use in sampling
#'   predictor percentages.
#' @param z A matrix or data frame containing geography-level covariates, or a
#'   vector of singular values to use in generating the random covariates.
#'   Singular values default to `exp(-(seq_len(p) - 1)/2)`. They are applied
#'   to scale a random orthogonal matrix.
#' @param b_loc The center of the distribution of geography-level parameters.
#'   Defaults to a linearly spaced sequence across groups from 0.5 to 0.9.
#' @param coefs Either A `p`-by-`n_x` matrix of coefficients for `z`, or
#'   The product with `z` is added to `ctr` to form the mean of the individual
#'   geography-level parameters. Defaults to a matrix of standard Normals.
#' @param b_cov The covariance matrix for geography-level parameters.
#'   Defaults to `0.02 * (1 + diag(n_x))`.
#'
#' @returns An `ei_spec` object with additional attributes `b_loc`, `coefs`, and
#'   `b_cov`, `eta` (the linear predictor for `b`), and `est_true` (the mean of
#'   the geography-level parameters, formatted similarly to the output from
#'   [ei_est()]).
#'
#' @examples
#' ei_synthetic(n = 10)
#'
#' ei_synthetic(n = 10, p = 2, n_x = 3)
#'
#' data(elec_1968)
#' ei_synthetic(
#'     x = cbind(elec_1968$pop_white, 1 - elec_1968$pop_white),
#'     b_loc = c(0.3, 0.9),
#'     b_cov = matrix(c(0.02, 0.016, 0.016, 0.2), nrow=2)
#' )
#'
#' @export
ei_synthetic = function(n, p = 0, n_x = 2, x = n_x:1, z = NULL,
                        b_loc = NULL, coefs = NULL, b_cov = NULL) {
    if (!is.matrix(x) && !is.data.frame(x)) {
        if (missing(n))
            cli_abort("If {.arg x} is not a matrix or data frame, you must specify {.arg n}.")

        x = matrix(rgamma(n*length(x), shape=rep(x, each=n)),
                   nrow=n, ncol=length(x))
        x = x / rowSums(x)
    } else {
        x = as.matrix(x)
    }
    n = nrow(x)
    n_x = ncol(x)

    if (!is.matrix(z) && !is.data.frame(z)) {
        if (is.null(z)) {
            z = exp(-(seq_len(p) - 1)/2)
        }
        z0 = matrix(rnorm(n * p), nrow=n, ncol=p)
        if (p > 0) {
            udv = svd(qr.Q(qr(z0)))
            udv$d = z
            z = tcrossprod(udv$u, scale_cols(udv$v, udv$d))
        } else {
            z = z0
        }
    } else {
        z = as.matrix(z)
    }
    p = ncol(z)

    if (is.null(b_loc)) {
        b_loc = seq(0.5, 0.9, length.out=n_x)
    }
    if (is.null(coefs)) {
        coefs = matrix(rnorm(p * n_x), nrow=p, ncol=n_x)
    } else {
        if (!all(dim(coefs) == c(p, n_x))) {
            cli_abort("If {.arg coefs} is provided, it must be a {p}-by-{n_x} matrix.")
        }
    }
    if (is.null(b_cov)) {
        b_cov = 0.02 * (1 + diag(n_x))
    }

    eta = shift_cols(z %*% coefs, -b_loc)

    warmup = 10L
    b = matrix(nrow=n, ncol=n_x)
    L = t(chol(b_cov))
    for (i in seq_len(n)) {
        b[i, ] = R_ess_tmvn(warmup, eta[i, ], L, init=eta[i, ])[warmup, ]
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

    new_tibble(
        d,
        ei_x = colnames(x),
        ei_y = "y",
        ei_z = colnames(z),
        ei_wgt = rep(1, n),
        b = b,
        b_loc = b_loc,
        coefs = coefs,
        b_cov = b_cov,
        eta = eta,
        est_true = true,
        class = "ei_spec"
    )
}
