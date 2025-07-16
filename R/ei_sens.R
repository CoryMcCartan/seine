#' Conduct a sensitivity analysis for estimated ecological quantities
#'
#' Relates confounding of an omitted variable with predictor or outcome to
#' bias in ecological estimates, using the nonparametric sensitivity analysis
#' of Chernozhukov et al. (2022).
#'
#' The parameter `c_predictor` equals \eqn{1 - R^2_{\alpha\sim\alpha_s}}, where
#' \eqn{\alpha} is the true Riesz representer and \eqn{\alpha_s} is the Riesz
#' representer with the observed covariates. The RR can be equivalently
#' expressed as \deqn{
#'   \alpha = \partial_x \log f(X\mid Z, U),
#' } where \eqn{U} is the unobserved confounder and \eqn{f} is the conditional
#' density. The corresponding `c_predictor` is then \deqn{
#'   1 - R^2_{\alpha\sim\alpha_s} = 1 - \
#'   \frac{\mathbb{E}[(\partial_x \log f(X\mid Z))^2]}{
#'   \mathbb{E}[(\partial_x \log f(X\mid Z, U))^2]}.
#' } When \eqn{X\mid Z,U} and \eqn{X\mid Z} are homoscedastic Gaussian, this
#' simplifies to \deqn{
#'   1 - R^2_{\alpha\sim\alpha_s} =
#'   1 - \frac{\mathbb{E}[X - \mathbb{E}[X\mid Z, U]]^2}{
#'   \mathbb{E}[X - \mathbb{E}[X\mid Z]]^2}
#'   = R^2_{X\sim U\mid Z}.
#' }
#'
#' The bounds here are plug-in estimates and do not incorporate sampling
#' uncertainty. As such, they may fail to cover the true value in finite
#' samples, even under large enough sensitivity parameters; see Section 5 of
#' Chernozhukov et al. (2022).
#'
#' @param est A set of estimates from [ei_est()] using both regression and Riesz
#'   representer.
#' @param c_outcome The (nonparametric) partial \eqn{R^2} of the omitted
#'   variables with the outcome variables. Must be between 0 and 1.
#'   Can be a vector, in which case all combinations of values with
#'   `c_predictor` are used.
#' @param c_predictor How much variation latent variables create in the Riesz
#'   representer, i.e. \eqn{1-R^2} of the true Riesz representer on the
#'   estimated one without the omitted variable. Must be between 0 and 1.
#'   Can be a vector, in which case all combinations of values with `c_outcome`
#'   are used.
#' @param bias_bound If provided, overrides `c_predictor` and finds values of
#'   `c_predictor` that correspond to (the absolute value of) the provided
#'   amount of bias.
#' @param confounding The confounding parameter (\eqn{\rho}), which must be
#'   between 0 and 1 (the adversarial worst-case).
#' @param expand_ci If `TRUE` and confidence intervals are present in `est`,
#'   expand the width of the intervals in each direction by the calculated bias
#'   bound.
#'
#' @returns A data frame of the same format as `est`, but with additional
#'   columns: `c_outcome` and `c_predictor`, matching all combinations of those
#'   arguments, and `bias_bound`, containing the bound on the amount of bias.
#'   The data frame has additional class `ei_sens`, which supports a
#'   [plot.ei_sens()] method.
#'
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2022).
#' *Long story short: Omitted variable bias in causal machine learning*
#' (No. w30302). National Bureau of Economic Research.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#' est = ei_est(m, rr, spec)
#'
#' ei_sens(est, c_outcome=0.2)
#'
#' # How much variation would the regression residual need to explain of
#' # Riesz representer to cause bias of 0.4?
#' ei_sens(est, c_outcome=1, bias_bound=0.4)
#'
#' # Update confidence intervals and extract as matrix
#' est = ei_est(m, rr, spec, conf_level=0.95)
#' sens = ei_sens(est, c_outcome=0.5, c_predictor=0.2)
#' as.matrix(sens, "conf.high")
#'
#' @export
ei_sens <- function(est, c_outcome=seq(0, 1, 0.01)^2, c_predictor=seq(0, 1, 0.01)^2,
                    bias_bound=NULL, confounding=1, expand_ci=TRUE) {
    is_01 <- function(value, arg) {
        if (!is.numeric(value) || all(value < 0) || all(value > 1)) {
            cli_abort("{.arg {arg}} must be between 0 and 1.",
                      call = parent.frame())
        }
    }

    is_01(c_outcome, "c_outcome")
    is_01(c_predictor, "c_predictor")
    is_01(confounding, "confounding")
    if (length(confounding) != 1) {
        cli_abort("{.arg confounding} must be a single value between 0 and 1.",
                  call = parent.frame())
    }

    idx = as.matrix(as.data.frame(est[, c("predictor", "outcome")]))
    sens_s = attr(est, "sens_s")
    bounds_inf = attr(est, "bounds_inf")
    if (is.null(bias_bound)) {
        cc = expand.grid(c_outcome=c_outcome, c_predictor=c_predictor)
        est = merge(est, cc)
        est$bias_bound = sens_s[idx] * confounding *
            sqrt(est$c_outcome * est$c_predictor / (1 - est$c_predictor))
    } else {
        if (!is.numeric(bias_bound)) {
            cli_abort("If provided, {.arg bias_bound} must be a numeric vector.",
                      call = parent.frame())
        }

        cc = expand.grid(c_outcome=c_outcome, c_predictor=1, bias_bound=abs(bias_bound))
        est = merge(est, cc)

        cp = est$bias_bound^2 / (sens_s[idx]^2 * confounding^2 * est$c_outcome)
        est$c_predictor = cp / (1 + cp)
    }

    if (isTRUE(expand_ci) && all(c("conf.low", "conf.high") %in% names(est))) {
        est$conf.low = est$conf.low - est$bias_bound
        est$conf.high = est$conf.high + est$bias_bound
    }

    tibble::new_tibble(est, bounds_inf=bounds_inf, sens_s=sens_s,
                       class=c("ei_sens", "ei_est"))
}

#' Robustness values for ecological inference
#'
#' TODO fill in...
#'
#' @param bias_bound <[`data-masking`][rlang::args_data_masking]> A bias bound:
#'   an amount of bias which is considered substantial. Evaluated in the context
#'   of `est`, so that one can to refer to `std.error` and `estimate` as needed.
#' @inheritParams ei_sens
#'
#' @returns A data frame of the same format as `est`, but with a new `rv` column
#'   containing the robustness values.
#'
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2022).
#' *Long story short: Omitted variable bias in causal machine learning*
#' (No. w30302). National Bureau of Economic Research.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#' est = ei_est(m, rr, spec)
#'
#' ei_sens_rv(est, 0.1) # how much confounding for bias of 0.1
#' ei_sens_rv(est, 2 * std.error) # how much confounding for bias of 2 SE
#'
#' # How much confounding to equalize all estimates (no polarization)
#' y_avg = weighted.mean(elec_1968$pres_ind_wal, elec_1968$pres_total)
#' ei_sens_rv(est, estimate - y_avg)
#'
#' # Extract as matrix
#' as.matrix(ei_sens_rv(est, 0.2), "rv")
#'
#' @export
ei_sens_rv <- function(est, bias_bound, confounding=1) {
    bias_bound = eval_tidy(enquo(bias_bound), est)
    if (!is.numeric(bias_bound)) {
        cli_abort("{.arg bias_bound} must be numeric.", call = parent.frame())
    }

    idx = as.matrix(as.data.frame(est[, c("predictor", "outcome")]))
    a = bias_bound^2 / attr(est, "sens_s")[idx]^2 / confounding^2
    est$rv = 0.5 * (-a + sqrt(a^2 + 4*a))

    est
}

#' Bias contour plot for ecological inference estimates
#'
#' TODO fill in...
#'
#' @param x An [ei_sens] object
#' @param y An outcome variable, as a character vector. Defaults to first.
#' @param predictor A predictor variable to plot, as a character vector. Defaults to first.
#' @param bounds A vector `c(min, max)` of bounds for the outcome, which will
#'   affect the contours which are plotted. In general, truncation will lead to
#'   violations of the accounting identity. If `bounds = NULL` (the default),
#'   they will be inferred from the outcome variable: if it is contained within
#'   \eqn{[0, 1]}, for instance, then the bounds will be `c(0, 1)`. Setting
#'   `bounds = FALSE` forces unbounded estimates.
#' @param plot_se A vector of multiples of the standard error to plot as contours.
#' @param ... Additional arguments passed on to [contour()]
#' @param lwd A scaling factor for the contour line widths
#'
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2022).
#' *Long story short: Omitted variable bias in causal machine learning*
#' (No. w30302). National Bureau of Economic Research.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
#'                total = pres_total, covariates = c(state, pop_urban, farm))
#' m = ei_ridge(spec)
#' rr = ei_riesz(spec, penalty = m$penalty)
#' est = ei_est(m, rr, spec)
#' sens = ei_sens(est)
#'
#' plot(sens)
#'
#' @export
plot.ei_sens <- function(x, y=NULL, predictor=NULL, bounds=NULL, plot_se=1:3, ..., lwd=1) {
    if (is.null(y)) y = x$outcome[1]
    if (is.null(predictor)) predictor = x$predictor[1]

    if (!is.character(y) && length(y) != 1) {
        cli_abort("{.arg y} must be a character vector with a single outcome name.",
                  call = parent.frame())
    }
    if (!y %in% x$outcome) {
        cli_abort("{.arg y} must be one of the outcomes in {.arg x}.",
                  call = parent.frame())
    }
    if (!is.character(predictor) && length(predictor) != 1) {
        cli_abort("{.arg predictor} must be a character vector with a single predictor name.",
                  call = parent.frame())
    }
    preds = unique(x$predictor)
    if (!predictor %in% preds) {
        cli_abort("{.arg predictor} must be one of the predictors in {.arg x}.",
                  call = parent.frame())
    }

    x = x[x$outcome == y & x$predictor == predictor, ]
    cx = unique(x$c_outcome)
    cy = unique(x$c_predictor)

    if (!all(diff(cx) > 0) || !all(diff(cy) > 0)) {
        cli_abort("{.arg x} must be sorted by {.var c_outcome} and {.var c_predictor}.",
                  call = parent.frame())
    }
    cz = matrix(x$bias_bound, nrow=length(cx), byrow = TRUE)

    if (is.null(bounds)) {
        bounds = ei_bounds(attr(x, "bounds_inf"), NULL)
    }
    bounds[is.infinite(bounds)] = range(x$estimate)[is.infinite(bounds)]

    n_om = 3 # orders of magnitude
    breaks = diff(bounds) * (10^-seq_len(n_om) %x% c(2:4, 6:9))
    oldmar = graphics::par()$mar
    graphics::par(mar = c(4.2, 5.2, 3, 1.1))
    graphics::contour(
        cx, cy, cz, levels=breaks, drawlabels=FALSE, col="#bbb", lwd=lwd,
        xlab = bquote(1 - {R^2}[alpha ~ "~" ~ alpha[s]]),
        ylab = bquote({R^2}[.(y) ~ "~ confounder |" ~
                            .(paste(preds, collapse = " + ")) ~ " + covariates" ]),
        main = paste0("Sensitivity bounds for E[", y, " | ", predictor, "]"),
        xaxs="i", yaxs="i", cex.lab=1.5
    )
    graphics::grid()

    breaks = c(10^-seq_len(n_om) %x% c(1, 5), 1)
    labels = as.character(breaks)
    special = c(abs(bounds - x$estimate[1]), x$std.error[1] * plot_se)
    dists = apply(abs(outer(special, breaks, `/`) - 1), 2, min)
    labels[dists < 0.05] = ""
    graphics::contour(
        cx, cy, cz, levels=breaks, labels=labels,
        lwd = lwd * c(rep(c(1.6, 1.0), n_om), 1.6),
        labcex=0.8, col = "#444", add=TRUE, method="edge"
    )
    graphics::contour(
        cx, cy, cz, lwd=2*lwd, labcex=1.0, col="#a42",
        levels = abs(bounds - x$estimate[1]),
        labels = paste("Estimate =", bounds),
        add=TRUE, method="edge"
    )
    if (length(plot_se) > 0) {
        graphics::contour(
            cx, cy, cz, lwd=2*lwd, lty="dashed", labcex=1.0, col="#46b",
            levels = x$std.error[1] * plot_se,
            labels = paste("\u00b1", plot_se, "SE"),
            add=TRUE, method="edge"
        )
    }
    graphics::par(mar = oldmar)
}


