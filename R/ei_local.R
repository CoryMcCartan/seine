#' Produce local ecological estimates
#'
#' Projects predictions from a fitted regression model onto the accounting
#' constraint using a provided residual covariance matrix. This ensures that
#' each set of local estimates satisfies the accounting identity. Local
#' estimates may be truncated to variable bounds.
#'
#' Local estimates are produced independently for each outcome variable.
#' Truncation to bounds, if used, will in general lead to estimates that do
#' not satisfy the accounting identity.
#'
#' @param regr A fitted regression model, from [ei_ridge()], or another kind
#'    of regression model wrapped with [ei_wrap_model()].
#' @param data The data frame, matrix, or [ei_spec] object that was used to fit
#'   the regression.
#' @param r_cov A covariance matrix of the residuals to use in projecting the
#'   local estimates onto the accounting constraint, or a list of matrices, one
#'   for each outcome variable. Defaults to the identity matrix scaled by the
#'   residual variance of `regr`, corresponding to orthogonal projection. Set
#'   `r_cov=1` to use a degenerate covariance matrix corresponding to a (local)
#'   neighborhood model. When there are multiple outcome variables and `r_cov` is
#'   a matrix, it will be applied identically to each outcome.
#' @param bounds A vector `c(min, max)` of bounds for the outcome, to which the
#'   local estimates will be truncated. In general, truncation will lead to
#'   violations of the accounting identity. If `bounds = NULL`, they will be
#'   inferred from the outcome variable: if it is contained within \eqn{[0, 1]},
#'   for instance, then the bounds will be `c(0, 1)`. Setting `bounds = FALSE`
#'   forces unbounded estimates.
#' @param conf_level A numeric specifying the level for confidence intervals.
#'   If `FALSE` (the default), no confidence intervals are calculated.
#'   For `regr` arguments from [ei_wrap_model()], confidence intervals will not
#'   incorporate uncertainty in the prediction itself, just the residual. This
#'   will trigger a warning periodically.
#' @param unimodal If `TRUE`, assume a unimodal residual distribution. Improves
#'   width of confidence intervals by a factor of 4/9.
#'
#' @returns A data frame with estimates. The `.row` column in the output
#'   corresponds to the observation index in the input. It has class
#'   `ei_est_local`, supporting several methods.
#'
#' @examples \dontrun{
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_dem_hum:pres_abs,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#'
#' m = ei_ridge(spec)
#'
#' ei_est_local(m, spec, conf_level = 0.95)
#' suppressWarnings(ei_est_local(m, spec, bounds=c(0.01, 0.2)))
#' }
#' @export
ei_est_local = function(regr, data, r_cov=NULL, bounds=NULL, conf_level=FALSE, unimodal=TRUE) {
    y = est_check_outcome(regr, data, NULL)
    n = nrow(y)
    n_y = ncol(y)

    cli_warn("Local confidence intervals do not yet incorporate prediction uncertainty.",
             .frequency="regularly", .frequency_id="ei_est_local_temp")

    rl = est_check_regr(regr, data, n, NULL, n_y, sd = TRUE)
    rl <<- rl
    n_x = length(rl$preds)
    if (inherits(regr, "ei_wrapped") && !isFALSE(conf_level)) {
        cli_warn("Local confidence intervals with wrapped model objects
                  do not incorporate prediction uncertainty.",
                 .frequency="regularly", .frequency_id="ei_est_local")
    }

    bounds = ei_bounds(bounds, y)

    # Process r_cov; TODO: heteroskedastic model
    if (is.null(r_cov)) {
        r_cov = lapply(regr$sigma2, function(s2) s2 * diag(n_x))
    } else if (length(r_cov) == 1 && r_cov == 1) {
        r_cov = lapply(regr$sigma2, function(s2) s2 * (1 + diag(n_x) * 1e-8))
    }
    if (!is.list(r_cov)) {
        r_cov = lapply(seq_len(n_y), function(i) r_cov)
    }
    for (r in r_cov) {
        if (!is.matrix(r) || nrow(r) != ncol(r) || nrow(r) != n_x) {
            cli_abort("Invalid {.arg r_cov} found.")
        }
    }
    r_cov = lapply(r_cov, chol)

    ests = list()
    for (k in seq_len(n_y)) {
        eta = vapply(rl$preds, function(p) p[, k], numeric(n))
        eta <<- eta
        eta_proj = local_proj(rl$x, eta, y[, k] - rl$yhat[, k], r_cov[[k]], bounds)
        eta_proj <<- eta_proj

        ests[[k]] = tibble::new_tibble(list(
            .row = rep(seq_len(n), n_x),
            predictor = rep(colnames(rl$x), each=n),
            outcome = rep(colnames(y)[k], n * n_x),
            estimate = c(eta_proj),
            std.error = NA #sqrt(c(proj[[2]]))
        ), class="ei_est_local")
    }

    ests = do.call(rbind, ests)

    if (!isFALSE(conf_level)) {
        fac = if (isTRUE(unimodal)) 4/9 else 1
        chebyshev = sqrt(fac / (1 - conf_level))
        ests$conf.low = ests$estimate - chebyshev * ests$std.error
        ests$conf.high = ests$estimate + chebyshev * ests$std.error

        ests$conf.low[ests$conf.low < bounds[1]] = bounds[1]
        ests$conf.high[ests$conf.high < bounds[1]] = bounds[1]
        ests$conf.low[ests$conf.low > bounds[2]] = bounds[2]
        ests$conf.high[ests$conf.high > bounds[2]] = bounds[2]
    }
    ests$std.error = NULL

    ests
}

#' @describeIn ei_est_local Format estimates an array with dimensions
#'   `<rows>*<predictors>*<outcomes>`. Does not work if the object has been sorted.
#' @param x An object of class `ei_est_local`
#' @param ... Additional arguments (ignored)
#' @export
as.array.ei_est_local = function(x, ...) {
    nm_x = unique(x$predictor)
    nm_y = unique(x$outcome)
    n_x = length(nm_x)
    n_y = length(nm_y)
    n = nrow(x) / n_x / n_y
    array(x$estimate, dim=c(n, n_x, n_y), dimnames=list(NULL, nm_x, nm_y))
}

# Solve QP to project estimates onto tomography plane and into bounds
# Not the fastest possible implementation (pure C++ would be better), but fast enough
local_proj = function(x, eta, eps, r_cov, bounds) {
    n = nrow(eta)
    n_x = ncol(x)
    eta_diff = matrix(nrow = n, ncol = n_x)

    zeros = rep(0, n_x)
    Amat = cbind(zeros)
    b0 = cbind(eps)
    if (!is.infinite(bounds[1])) {
        Amat = cbind(Amat, diag(n_x))
        b0 = cbind(b0, bounds[1] - eta)
    }
    if (!is.infinite(bounds[2])) {
        Amat = cbind(Amat, -diag(n_x))
        b0 = cbind(b0, -bounds[2] + eta)
    }

    for (i in seq_len(n)) {
        Amat[, 1] = x[i, ]
        eta_diff[i, ] = tryCatch({
            quadprog::solve.QP(
                Dmat = r_cov,
                dvec = zeros,
                Amat = Amat,
                bvec = b0[i, ],
                meq = 1,
                factorized = TRUE
            )$solution
        }, error = \(e) eps[i])
    }

    eta + eta_diff
}

local_basis = function(x) {
    qr.Q(qr(cbind(x, diag(length(x)))))[,  drop=FALSE]
}

