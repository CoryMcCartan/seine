#' Conduct a sensitivity analysis for estimated ecological quantities
#'
#' Relates confounding of an omitted variable with predictor or outcome to
#' bias in ecological estimates, using the nonparametric sensitivity analysis
#' of Chernozhukov et al. (2024).
#'
#' The parameter `c_predictor` equals \eqn{1 - R^2_{\alpha\sim\alpha_s}}, where
#' \eqn{\alpha} is the true Riesz representer and \eqn{\alpha_s} is the Riesz
#' representer with the observed covariates. The RR can be equivalently
#' expressed as \deqn{
#'   \alpha = \partial_{\bar x_j} \log f(\bar x_j\mid z, u),
#' } where \eqn{U} is the unobserved confounder and \eqn{f} is the conditional
#' density. The corresponding `c_predictor` is then \deqn{
#'   1 - R^2_{\alpha\sim\alpha_s} = 1 - \
#'   \frac{\mathbb{E}[(\partial_{\bar x_j} \log f(\bar x_j\mid z))^2]}{
#'   \mathbb{E}[(\partial_{\bar x_j} \log f(\bar x_j\mid z, u))^2]}.
#' }
#'
#' The bounds here are plug-in estimates and do not incorporate sampling
#' uncertainty. As such, they may fail to cover the true value in finite
#' samples, even under large enough sensitivity parameters; see Section 5 of
#' Chernozhukov et al. (2024).
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
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2024).
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
#' # Works for contrasts as well
#' est = ei_est(m, rr, spec, contrast = list(predictor=c(1, -1, 0)))
#' ei_sens(est, c_outcome=0.5, c_predictor=0.5)
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
#' The robustness value is the minimum bound for both `c_outcome` and
#' `c_predictor` in [ei_sens()] such that the bias bound is a certain value.
#' For example, if the provided bias bound is 0.5, meaning a bias of magnitude
#' 0.5 would be considered problematic, then the robustness value is the minimum
#' amount of confounding of outcome and predictor (more specifically, the Riesz
#' representer) that would lead to bias of magnitude 0.5.
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
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2024).
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
#' # Works for contrasts as well
#' est = ei_est(m, rr, spec, contrast = list(predictor=c(1, -1, 0)))
#' ei_sens_rv(est, estimate) # how much to eliminate disparity
#'
#' @export
ei_sens_rv <- function(est, bias_bound, confounding=1) {
    if (missing(bias_bound)) {
        cli_abort("{.arg bias_bound} is required.")
    }
    bias_bound = eval_tidy(enquo(bias_bound), est)
    if (!is.numeric(bias_bound)) {
        cli_abort("{.arg bias_bound} must be numeric.")
    }

    idx = as.matrix(as.data.frame(est[, c("predictor", "outcome")]))
    a = bias_bound^2 / attr(est, "sens_s")[idx]^2 / confounding^2
    est$rv = 0.5 * (-a + sqrt(a^2 + 4*a))

    est
}

#' Bias contour plot for ecological inference estimates
#'
#' Displays bias bound as a function of `c_outcome` and `c_predictor` in
#' [ei_sens()] on a contour plot. Bounds on the outcome, and standard errors of
#' the point estimate, can be overlaid as contours on the plot to aid in
#' interpretation. Benchmarked values of `c_outcome` and `c_predictor` based on
#' the observed covariates can also be overlaid.
#'
#' @param x An [ei_sens] object
#' @param y An outcome variable, as a character vector. Defaults to first.
#' @param predictor A predictor variable to plot, as a character vector. Defaults to first.
#' @param bounds A vector `c(min, max)` of bounds for the outcome, which will
#'   affect the contours which are plotted. If `bounds = NULL` (the default),
#'   they will be inferred from the outcome variable: if it is contained within
#'   \eqn{[0, 1]}, for instance, then the bounds will be `c(0, 1)`. Setting
#'   `bounds = FALSE` forces unbounded estimates.
#' @param bench A data frame of benchmark values, from [ei_bench()], to plot.
#' @param plot_se A vector of multiples of the standard error to plot as contours.
#' @param contour_exp Powers of 10 for which to plot contours of the bias bound.
#' @param ... Additional arguments passed on to [contour()]
#' @param lwd Scaling factor for the contour line widths
#' @param cex Scaling factor for the benchmark points and labels, if provided
#' @param pch The point type (see [points()]) for the benchmark values, if
#'   provided
#'
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2024).
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
#' plot(sens, bench = ei_bench(spec), plot_se=NULL)
#' @export
plot.ei_sens <- function(x, y=NULL, predictor=NULL, bounds=NULL, bench=NULL,
                         plot_se=1:3, contour_exp=-2:-1, ..., lwd=1, pch=8, cex=1) {
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

    if (is.null(bounds)) {
        bounds = ei_bounds(attr(x, "bounds_inf"), NULL)
    }
    bounds[is.infinite(bounds)] = range(x$estimate)[is.infinite(bounds)]

    x = x[x$outcome == y & x$predictor == predictor, ]
    cx = unique(x$c_outcome)
    cy = unique(x$c_predictor)

    if (!all(diff(cx) > 0) || !all(diff(cy) > 0)) {
        cli_abort("{.arg x} must be sorted by {.var c_outcome} and {.var c_predictor}.",
                  call = parent.frame())
    }
    cz = matrix(x$bias_bound, nrow=length(cx), byrow = TRUE)

    breaks = 10^contour_exp %x% c(2:4, 6:9)
    oldmar = graphics::par()$mar
    graphics::par(mar = c(4.2, 5.2, 3, 1.1))
    graphics::contour(
        cx, cy, cz, levels=breaks, drawlabels=FALSE, col="#d0d0d0", lwd=lwd,
        xlab = bquote(1 - {R^2}[alpha ~ "~" ~ alpha[s]]),
        ylab = bquote({R^2}[.(y) ~ "~ confounder | predictors, covariates"]),
        # ylab = bquote({R^2}[.(y) ~ "~ confounder |" ~
        #                .(paste(preds, collapse = ", ")) ~ ", covariates" ]),
        main = paste0("Sensitivity bounds for E[", y, " | ", predictor, "]"),
        xaxs="i", yaxs="i", cex.lab=1.5
    )
    graphics::grid(col="#dfdfdf")

    breaks = c(10^contour_exp %x% c(1, 5), 1)
    labels = as.character(breaks)
    special = c(abs(bounds - x$estimate[1]), x$std.error[1] * plot_se)
    dists = apply(abs(outer(special, breaks, `/`) - 1), 2, min)
    labels[dists < 0.05] = ""
    graphics::contour(
        cx, cy, cz, levels=breaks, labels=labels,
        lwd = lwd * c(rep(c(1.4, 1.0), length(contour_exp)), 1.4),
        labcex=0.8, col = "#666", add=TRUE, method="edge"
    )
    graphics::contour(
        cx, cy, cz, lwd=2*lwd, labcex=0.85*cex, col="#a42",
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
    if (!missing(bench)) {
        if (!inherits(bench, "ei_bench")) {
            cli_abort("{.arg bench} must be an {.cls ei_bench} object.",
                      call = parent.frame())
        }

        bench = bench[bench$outcome == y & bench$predictor == predictor, ]
        graphics::points(bench$c_predictor, bench$c_outcome,
                         col="#a42", pch=pch, cex=1.5*cex)
        tpos = ifelse(bench$c_predictor < bench$c_outcome, 4, 3)
        graphics::text(bench$c_predictor, bench$c_outcome, bench$covariate,
                       pos=tpos, cex=0.85*cex, font=2)
    }
    graphics::par(mar = oldmar)
}

#' Benchmark sensitivity parameters from observed covariates
#'
#' Produces a table of benchmark values for `c_outcome` and `c_predictor` in
#' [ei_sens()] for each covariate, following the methodology of Chernozhukov
#' et al. (2024).
#'
#' @param spec An [ei_spec] object.
#' @param preproc An optional function which takes in a `ei_spec` object (`spec`
#'   with one covariate removed) and returns a modified object that includes
#'   modified object. Useful to apply any preprocessing, such as a basis
#'   transformation, as part of the benchmarking process.
#' @param subset Passed on to [ei_est()].
#'
#' @references
#' Chernozhukov, V., Cinelli, C., Newey, W., Sharma, A., & Syrgkanis, V. (2024).
#' *Long story short: Omitted variable bias in causal machine learning*
#' (No. w30302). National Bureau of Economic Research.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(elec_1968, vap_white:vap_other, pres_ind_wal,
#'                total = pres_total, covariates = c(educ_elem, pop_urban, farm))
#'
#' ei_bench(spec)
#'
#' # preprocess to add all 2-way interactions
#' ei_bench(spec, preproc = function(s) {
#'     z_cols = match(attr(s, "ei_z"), names(s))
#'     s_out = s[-z_cols]
#'     z_new = model.matrix(~ .^2 - 1, data = s[z_cols])
#'     s_out = cbind(s_out, z_new)
#'     ei_spec(s_out, vap_white:vap_other, pres_ind_wal,
#'             total = attr(s, "ei_n"), covariates = colnames(z_new))
#' })
#' @export
ei_bench <- function(spec, preproc = NULL, subset = NULL) {
    validate_ei_spec(spec)

    if (!missing(preproc)) {
        if (!is.function(preproc)) {
            cli_abort("{.arg preproc} must be a function.")
        }
    } else {
        preproc = function(x) x
    }

    n_x = length(attr(spec, "ei_x"))
    n_y = length(attr(spec, "ei_y"))
    subs = eval_tidy(enquo(subset), spec)
    var_resid = function(regr) {
        apply(regr$y - regr$fitted, 2, var)
    }

    spec_proc = preproc(spec)
    regr0 = ei_ridge(spec_proc, vcov=FALSE)
    riesz0 = ei_riesz(spec_proc, penalty=regr0$penalty)
    est0 = ei_est(regr0, riesz0, spec_proc, subset = subs)
    vy = apply(regr0$y, 2, var)
    var_resid0 = var_resid(regr0)
    r2_out0 = 1 - var_resid0 / vy

    covs = attr(spec, "ei_z")
    benches = lapply(covs, function(cv) {
        spec_loo = reconstruct_ei_spec(spec[setdiff(names(spec), cv)], spec)
        spec_loo = preproc(spec_loo)

        regr_loo = ei_ridge(spec_loo, vcov=FALSE)
        riesz_loo = ei_riesz(spec_loo, penalty=regr_loo$penalty)
        est_loo = ei_est(regr_loo, riesz_loo, spec_loo, subset = subs)
        var_resid_loo = var_resid(regr_loo)
        r2_out_loo = 1 - var_resid_loo / vy
        r2_riesz = riesz_loo$nu2 / riesz0$nu2

        c_outcome = pmin((r2_out0 - r2_out_loo) / r2_out0, 1)
        c_predictor = pmin((1 - r2_riesz) / r2_riesz, 1)
        est_chg = est_loo$estimate - est0$estimate
        sd_diff = sqrt(pmax(var_resid_loo - var_resid0, 0))
        nu_diff = sqrt(pmax(riesz0$nu2 - riesz_loo$nu2, 0))
        confounding = est_chg / rep(sd_diff, each=n_x) / rep(nu_diff, n_y)
        confounding = pmax(pmin(confounding, 1), -1)

        est_loo$covariate = cv
        est_loo = est_loo[c("covariate", "predictor", "outcome")]
        est_loo$c_outcome = rep(c_outcome, each=n_x)
        est_loo$c_predictor = rep(c_predictor, n_y)
        est_loo$confounding = confounding
        est_loo$est_chg = est_chg
        est_loo
    })

    tibble::new_tibble(do.call(rbind, benches), class="ei_bench")
}

