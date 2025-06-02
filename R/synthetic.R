#' Generate Synthetic Ecological Data
#'
#' Samples data from the following truncated Normal ecological model:
#' \deqn{
#'      x_i \stackrel{\text{iid}}{\sim} \mathrm{Dirichlet}(\alpha) \\
#'      b_i \stackrel{\text{iid}}{\sim} \mathcal{N}_{[0, 1]^d}(\eta, \Sigma) \\
#'      y_i = b_i^\top x_i,
#' }
#' where \eqn{\alpha} is a vector of hyperparameters that may be optionally
#' specified as the `x` argument to the function, and \eqn{\mathcal{N}_{[0,
#' 1]^d}} indicates a multivariate Normal distribution truncated to the
#' *d*-dimensional hypercube.
#'
#'
#' @param x Either a matrix or data frame containing the predictor percentages
#'   in each row, or a vector containing Dirichlet parameters to use in sampling
#'   predictor percentages.
#' @param eta The center of the distribution of geography-level parameters.
#'   Defaults to a sequence across groups.
#' @param Sigma The covariance matrix for geography-level parameters.
#' @param n The number of rows (geographies) to generate. Defaults to the number
#'   of rows in `x`, if the latter is a matrix or data frame.
#'
#' @returns A list with elements `x` (predictors), `y` (outcome), `b`
#'   (geography-level parameters), and `eta` and `Sigma` (global parameters).
#'
#' @examples
#' ei_synthetic(n=100)
#'
#' data(elec_1968)
#' ei_synthetic(cbind(elec_1968$pop_white, 1-elec_1968$pop_white),
#'      eta=c(0.3, 0.9), Sigma=matrix(c(0.02, 0.016, 0.016, 0.2), nrow=2))
#'
#' @export
ei_synthetic = function(x=c(3, 1), eta=seq(0.5, 0.9, length.out=ncol(x) %||% length(x)),
                        Sigma=0.02*(1 + diag(ncol(x) %||% length(x))), n=nrow(x)) {
    if (!is.matrix(x) && !is.data.frame(x)) {
        if (missing(n))
            cli_abort("If {.arg x} is not a matrix or data frame, you must specify {.arg n}.")

        x = matrix(rgamma(n*length(x), shape=rep(x, each=n)),
                   nrow=n, ncol=length(x))
        x = x / rowSums(x)
    }
    nx = ncol(x)

    warmup = 10L
    thin = 4L
    b = R_ess_tmvn(warmup + thin*n, eta, chol_pivot(Sigma), init=eta)
    b = b[warmup + seq(1, thin*n, by=4), ]

    list(
        x = x,
        y = rowSums(b * x),
        b = b,
        eta = eta,
        Sigma = Sigma
    )
}
