#' Generate semisynthetic ecological data from a fitted model
#'
#' Takes an [ei_spec] object and generates a new dataset with a synthetic
#' hidden confounder whose statistical properties are calibrated to the
#' sensitivity parameters `c_outcome` and `c_predictor`. The real predictors
#' and observed covariates are preserved unchanged; the outcome variable
#' is replaced with a simulated version driven by the confounded model.
#'
#' The function constructs a multivariate latent confounder \eqn{A} for each unit
#' such that:
#' - \eqn{A} is orthogonal (in the linear regression sense) to the observed
#'   covariates \eqn{Z}.
#' - The partial \eqn{R^2} of \eqn{A} with the Riesz representer weights
#'   (controlling for \eqn{Z}) is approximately `c_predictor`.
#' - The partial \eqn{R^2} of \eqn{A} with the simulated outcome (controlling
#'   for \eqn{Z}) is approximately `c_outcome`.
#'
#' Unit-level parameters \eqn{b_i} are drawn from a Normal distribution
#' centered at the original regression's fitted values plus a (calibrated)
#' linear effect of the confounder.
#' New outcomes are then computed as \eqn{y_i = \sum_j b_{ij} x_{ij}}.
#'
#' The `confounding` parameter is the \eqn{R^2} of the change in the regression
#' predictions when including the confounder with the change in the Riesz 
#' representer weights when including the confounder. Setting ' `confounding =
#' 1` (the default) produces the maximum adversarial bias consistent ' with the
#' specified `c_outcome` and `c_predictor`.
#'
#' @param regr A fitted regression model from [ei_ridge()].
#' @param spec An [ei_spec] object.
#' @param riesz A fitted Riesz representer from [ei_riesz()].
#' @inheritParams ei_est_local
#' @param predictor The predictor to target for calibration.  Must match an
#'   entry in the `predictor` column of `est` (which may be a contrast name).
#'   Defaults to the first predictor.
#' @param outcome The outcome to target for calibration.  Must match an entry
#'   in the `outcome` column of `est` (which may be a contrast name).  Defaults
#'   to the first outcome.
#' @param est An optional set of DML estimates from [ei_est()].
#'   If `NULL`, it will be computed automatically using `regr`, `riesz`, and
#'   `contrast`.
#' @param c_outcome The approximate partial \eqn{R^2} of the hidden confounder
#'   with the target outcome, controlling for the observed covariates. Must be
#'   between 0 and 1.
#' @param c_predictor The approximate partial \eqn{R^2} of the hidden confounder
#'   with the Riesz representer weights for the target predictor (controlling for
#'   the observed covariates), i.e., \eqn{1 - R^2_{\alpha \sim \alpha_s}}. Must
#'   be between 0 and 1.
#' @param confounding The correlation parameter between the confounder's effect
#'   on the predictor and the outcome, between 0 and 1. Scales `Lambda` as a
#'   whole.
#' @param include If `TRUE`, the returned `ei_spec` object will include the
#'   confounder as a covariate; if `FALSE`, the confounder will be a column in the
#'   `ei_spec` but not marked as a confounder.
#'
#' @returns An `ei_spec` object with:
#'   - The original predictor columns (unchanged).
#'   - New simulated outcome columns (replacing the original outcomes).
#'   - The original covariate columns (unchanged).
#'   - A `.confounder` column containing the hidden confounder values 
#'
#'   Additional attributes:
#'   - `b`: true unit-level parameters (\eqn{n \times n_x n_y} matrix, ordered
#'     as Y1|X1, Y1|X2, ..., Y2|X1, Y2|X2, ...).
#'   - `b_loc`: center of the b distribution (\eqn{n_x \times n_y} matrix).
#'   - `b_cov`: residual covariance of b (\eqn{n_x n_y \times n_x n_y} matrix).
#'   - `a`: the hidden confounder values (as a matrix with \eqn{n} rows).
#'   - `eta`: New linear predictor for b, including effect of `a`.
#'     (\eqn{n \times n_x n_y} matrix).
#'   - `c_outcome`, `c_predictor`, `confounding`: the specified parameters.
#'
#' @examples
#' data(elec_1968)
#'
#' spec = ei_spec(
#'     elec_1968, 
#'     predictors = vap_white:vap_other, 
#'     outcome = pres_dem_hum:pres_abs,
#'     total = pres_total,
#'     covariates = c(pop_urban, farm, educ_hsch, inc_25_99k),
#'     preproc = function(z) poly(model.matrix(~ 0 + ., z), degree = 2)
#' )
#' regr = ei_ridge(spec)
#' riesz = ei_riesz(spec, penalty = regr$penalty)
#'
#' spec_synth = ei_semisynthetic(
#'     regr, riesz, spec, 
#'     predictor = "vap_white", outcome = "pres_dem_hum", 
#'     c_outcome = 0.5, c_predictor = 0.5
#' )
#' 
#' # compute estimates on original and synthetic data
#' regr_synth = ei_ridge(spec_synth)
#' riesz_synth = ei_riesz(spec_synth, penalty = regr_synth$penalty)
#' est = ei_est(regr_synth, riesz_synth, spec_synth, conf_level = FALSE)
#' 
#' est_hum_white = subset(est, predictor == "vap_white" & outcome == "pres_dem_hum")
#' true_hum_white = subset(
#'     attr(spec_synth, "est_true"), 
#'     predictor == "vap_white" & outcome == "pres_dem_hum"
#' )$true
#' transform(
#'     ei_sens(est_hum_white, c_outcome = 0.5, c_predictor = 0.5),
#'     bias_act = abs(estimate - true_hum_white)
#' )
#' @export
ei_semisynthetic = function(
    regr,
    riesz,
    spec,
    predictor = NULL,
    outcome = NULL,
    c_outcome = 0.1,
    c_predictor = 0.1,
    confounding = 1,
    subset = NULL,
    contrast = NULL,
    include = FALSE,
    est = NULL,
    b_cov = NULL
) {
    validate_ei_spec(spec)
    if (!inherits(regr, c("ei_ridge", "ei_model"))) {
        cli_abort("{.arg regr} must be an {.cls ei_ridge} object.")
    }

    # Input validation
    for (arg in c("c_outcome", "c_predictor", "confounding")) {
        val = get(arg)
        if (!is.numeric(val) || length(val) != 1 || val < 0 || val > 1) {
            cli_abort("{.arg {arg}} must be a single number between 0 and 1.")
        }
    }

    # Extract basic dimensions
    y = est_check_outcome(regr, spec, NULL)
    n = nrow(y)
    n_y = ncol(y)
    rm = est_check_riesz(riesz, spec, NULL, n, regr)  # n x n_x, normalized
    rl = est_check_regr(regr, spec, n, colnames(rm), n_y, vcov = FALSE)

    n_x = length(rl$preds)
    x = rl$x  
    x_nm = colnames(x)
    y_nm = colnames(y)

    # Process subset: used to weight the sensitivity-parameter calibration
    subset_lgl = check_subset(eval_tidy(enquo(subset), spec), n)
    subset_cal = which(subset_lgl)

    total = attr(spec, "ei_n")
    weights = riesz$blueprint$ei_wgt
    if (is.null(weights)) weights = rep(1, n)
    w = weights / mean(weights)

    # Determine valid predictor/outcome names for this contrast configuration
    contr_info = check_contrast(contrast, x_nm, y_nm)
    valid_pred = unique(contr_info$x_nm)
    valid_out  = unique(contr_info$y_nm)

    if (is.null(predictor)) {
        predictor = valid_pred[1]
    } else if (!predictor %in% valid_pred) {
        cli_abort("{.arg predictor} must be one of {.val {valid_pred}}.")
    }
    if (is.null(outcome)) {
        outcome = valid_out[1]
    } else if (!outcome %in% valid_out) {
        cli_abort("{.arg outcome} must be one of {.val {valid_out}}.")
    }
    j_pred = match(predictor, valid_pred)
    k_out  = match(outcome, valid_out)

    # Predictor and outcome contrast vectors (lengths n_x, n_y)
    q_pred = if (!is.null(contrast$predictor)) {
        as.vector(contrast$predictor)
    } else {
        replace(rep(0, n_x), j_pred, 1L)
    }
    c_out = if (!is.null(contrast$outcome)) {
        as.vector(contrast$outcome)
    } else {
        replace(rep(0, n_y), k_out, 1L)
    }

    # Fit objects if not provided
    if (is.null(est)) {
        est = ei_est(regr, riesz, spec, contrast = contrast)
    } else if (!inherits(est, "ei_est")) {
        cli_abort("{.arg est} must be an {.cls ei_est} object.")
    }
    if (is.null(b_cov)) {
        b_cov = ei_local_cov(regr, spec)
    } else if (!is.matrix(b_cov)) {
        cli_abort("{.arg b_cov} must be a matrix.")
    }

    # --- Step 1: Residualize Riesz weights against Z with ridge ---
    # Use riesz$weights (unnormalized) for consistency with nu2 computations.
    alpha = riesz$weights  # n x n_x
    z_proc = as.matrix(attr(spec, "ei_z_proc"))
    p = ncol(z_proc)
    z_riesz = if (p > 0) {
        scale_cols(shift_cols(z_proc, riesz$z_shift), riesz$z_scale)
    } else {
        matrix(nrow = n, ncol = 0)
    }

    if (p > 0) {
        svd_z = svd(z_riesz * sqrt(w))
        d2 = svd_z$d^2
        ridge_d = d2 / (d2 + riesz$penalty)
        Utz = crossprod(svd_z$u, alpha * sqrt(w))
        hat_alpha = (svd_z$u %*% (ridge_d * Utz)) / sqrt(w)
        alpha_perp = alpha - hat_alpha
    } else {
        alpha_perp = alpha
    }

    # --- Step 2: Confounder direction for target predictor ---
    target_alpha = as.vector(alpha_perp %*% q_pred)  # n-vector
    norm_ta = sqrt(sum(target_alpha^2))
    u_a = if (norm_ta > 1e-12) target_alpha / norm_ta else target_alpha
    sd_ua = sd(u_a)
    u_a_v1 = if (sd_ua > 1e-14) u_a / sd_ua else u_a  # variance-1 signal direction

    # Orthogonal noise direction with variance 1 (used for dilution below)
    r = rnorm(n)
    r = r - c(crossprod(r, u_a_v1) / crossprod(u_a_v1)) * u_a_v1  # Gram-Schmidt
    sd_r = sd(r)
    e_perp = if (sd_r > 1e-14) r / sd_r else r

    # Diluted confounder: a(dil) = (u_a_v1 + dil * e_perp) normalised to variance 1.
    # dil=0 → pure signal (maximum c_predictor); dil→∞ → pure noise (zero c_predictor).
    # var(a(dil)) = 1 for all dil by construction.
    make_a_dil = function(dil) {
        v = u_a_v1 + dil * e_perp
        v / sd(v)
    }

    # --- Step 3: Binary search on dilution parameter to calibrate c_predictor ---
    # nu2 for target predictor before augmentation.
    # When a subset is used, refit Riesz on the subset so that nu2_short and
    # nu2_full (inside measure_c_pred) share the same Dz normalisation.
    if (all(subset_lgl)) {
        nu2_short_target = mean((alpha %*% q_pred)^2)
    } else {
        fit_short_sub = ei_riesz_impl(x[subset_cal, , drop=FALSE],
                                      z_riesz[subset_cal, , drop=FALSE],
                                      total[subset_cal], weights[subset_cal],
                                      riesz$penalty)
        nu2_short_target = mean((fit_short_sub$alpha %*% q_pred)^2)
    }

    measure_c_pred = function(dil) {
        a_dil = make_a_dil(dil)
        z_aug = cbind(z_riesz[subset_cal, , drop=FALSE], a_dil[subset_cal])
        fit_aug = ei_riesz_impl(x[subset_cal, , drop=FALSE], z_aug,
                                total[subset_cal], weights[subset_cal], riesz$penalty)
        nu2_full_target = mean((fit_aug$alpha %*% q_pred)^2)
        if (nu2_full_target <= 0) return(0)
        1 - nu2_short_target / nu2_full_target
    }

    if (c_predictor < 1e-10) {
        a = e_perp  # no signal direction; variance-1 pure-noise confounder
    } else {
        # c_pred(dil) decreases from c_pred_max at dil=0 toward ~0 as dil→∞
        c_pred_max = measure_c_pred(0)
        if (c_pred_max <= c_predictor) {
            dil_star = 0  # target exceeds max achievable at variance 1; use full signal
        } else {
            f = function(dil) measure_c_pred(dil) - c_predictor  # f(0) > 0, f(large) < 0
            dil_hi = 1
            max_hi = 1e12
            while (f(dil_hi) > 0 && dil_hi < max_hi) {
                dil_hi = dil_hi * 2
            }
            if (f(dil_hi) > 0) {
                dil_star = dil_hi
            } else {
                res = tryCatch(
                    uniroot(f, lower = 0, upper = dil_hi, tol = 1e-6, maxiter = 100),
                    error = function(e) NULL
                )
                dil_star = if (!is.null(res)) res$root else dil_hi
            }
        }
        a = make_a_dil(dil_star)  # variance-1 calibrated confounder
    }

    # --- Step 4: Build target outcome shift and calibrate b_noise ---
    # Draw the noise component of b first, then construct the target contrast
    # shift so that:
    #   1. cor(y_new[c_out] - yhat[c_out], a)^2 = c_outcome exactly, and
    #   2. cor(delta_fit, delta_riesz)^2 = confounding for the target contrast,
    # where delta_fit is the fitted-value change induced by A and delta_riesz is
    # the change in the target Riesz representer when augmenting Z with A.
    L_cov = chol(b_cov + diag(1e-10, nrow(b_cov)))
    b_noise = matrix(rnorm(n * n_x * n_y), nrow = n) %*% L_cov

    # eps_c: the noise component of (y[c_out] - yhat[c_out]) for target contrast.
    # eps_c[i] = sum_k c_out[k] * sum_j b_noise[i, j+(k-1)*n_x] * x[i,j]
    eps_c = numeric(n)
    for (k in seq_len(n_y)) {
        idx = seq_len(n_x) + (k - 1L) * n_x
        eps_c = eps_c + c_out[k] * rowSums(b_noise[, idx, drop = FALSE] * x)
    }

    c_out_ss = sum(c_out^2)
    sigma2_target = as.numeric(crossprod(c_out^2, regr$sigma2))

    # center_scale_cal: centres/scales using subset statistics, returns full-n vector.
    # When subset_cal == seq_len(n) this is identical to the naive centre_scale.
    center_scale_cal = function(v) {
        v_sub = v[subset_cal]
        v = v - mean(v_sub)
        s = sd(v_sub)
        if (s > 1e-14) v / s else v
    }
    cov_cal = function(u, v) cov(u[subset_cal], v[subset_cal])
    var_cal = function(v) var(v[subset_cal])
    sd_cal  = function(v) sd(v[subset_cal])
    cor_cal = function(u, v) cor(u[subset_cal], v[subset_cal])

    fit_aug = ei_riesz_impl(x, cbind(z_riesz, a), total, weights, riesz$penalty)
    delta_riesz = as.vector((fit_aug$alpha - alpha) %*% q_pred)

    if (confounding < 1e-10 || c_predictor < 1e-10 || c_outcome < 1e-10) {
        delta_fit = numeric(n)
    } else {
        a_std = center_scale_cal(a)
        rr_std = center_scale_cal(delta_riesz)
        rr_sd = sd_cal(rr_std)
        if (rr_sd < 1e-14) {
            cli_abort("Unable to calibrate {.arg confounding}: Riesz change is degenerate.")
        }

        r_ar = cor_cal(a_std, rr_std)
        a_perp_raw = a_std - r_ar * rr_std
        a_perp_sd = sd_cal(a_perp_raw)
        if (a_perp_sd > 1e-14) {
            a_perp = a_perp_raw / a_perp_sd
            if (cor_cal(a_perp, a_std) < 0) a_perp = -a_perp
        } else {
            a_perp = rr_std
        }

        # Build a unit-variance fitted-value direction whose squared correlation
        # with delta_riesz equals confounding, while remaining sufficiently
        # aligned with a to achieve the requested c_outcome after scaling.
        delta_fit_dir =
            sign(r_ar) * sqrt(confounding) * rr_std +
            sqrt(1 - confounding) * a_perp
        delta_fit_dir = center_scale_cal(delta_fit_dir)

        r_da = cor_cal(delta_fit_dir, a_std)
        if (r_da^2 + 1e-10 < c_outcome) {
            # Some extreme confounding requests are infeasible for the realised
            # A/RR geometry. Fall back to the fully outcome-aligned direction so
            # c_outcome remains calibrated exactly.
            delta_fit_dir = a_std
            r_da = 1
        }

        delta_fit = delta_fit_dir * sqrt(c_outcome * sigma2_target) / r_da
    }

    # Orthogonalize b_noise w.r.t. both a and delta_fit so that the target
    # contrast residual has the desired correlation with a and total variance.
    a_basis = center_scale_cal(a)
    fit_basis = center_scale_cal(delta_fit)
    fit_basis = fit_basis - cov_cal(fit_basis, a_basis) * a_basis

    proj = numeric(n)
    if (sd_cal(a_basis) > 1e-14) {
        proj = proj + cov_cal(eps_c, a_basis) * a_basis
    }
    if (sd_cal(fit_basis) > 1e-14) {
        proj = proj + cov_cal(eps_c, fit_basis) / var_cal(fit_basis) * fit_basis
    }
    if (any(abs(proj) > 0)) {
        for (k in seq_len(n_y)) {
            if (abs(c_out[k]) > 1e-10) {
                idx = seq_len(n_x) + (k - 1L) * n_x
                b_noise[, idx] = b_noise[, idx] -
                    (c_out[k] / c_out_ss) * matrix(proj, nrow = n, ncol = n_x)
            }
        }
    }

    # Recompute eps_c after orthogonalization.
    eps_c = numeric(n)
    for (k in seq_len(n_y)) {
        idx = seq_len(n_x) + (k - 1L) * n_x
        eps_c = eps_c + c_out[k] * rowSums(b_noise[, idx, drop = FALSE] * x)
    }
    var_eps_perp = var_cal(eps_c)

    var_eps_target = sigma2_target - var_cal(delta_fit)
    if (var_eps_target < -1e-10) {
        cli_abort(
            "Unable to jointly calibrate {.arg c_outcome} and {.arg confounding}: implied noise variance is negative."
        )
    }
    if (var_eps_perp > 1e-30 && var_eps_target > 0) {
        scale_eps = sqrt(var_eps_target / var_eps_perp)
        b_noise = b_noise * scale_eps
    }

    # --- Step 5: Build eta, apply signal, compute y_new ---
    # eta: n x (n_x * n_y), ordered as Y1|X1, Y1|X2, ..., Y2|X1, ...
    eta = matrix(nrow = n, ncol = n_x * n_y)
    for (k in seq_len(n_x)) {
        idx = k + seq(0, by = n_x, length.out = n_y)
        eta[, idx] = rl$preds[[k]]
    }

    eta_new = eta
    for (k in seq_len(n_y)) {
        idx = seq_len(n_x) + (k - 1L) * n_x
        eta_new[, idx] = eta[, idx] + (delta_fit * c_out[k] / c_out_ss)
    }
    colnames(eta_new) = c(outer(x_nm, y_nm, paste, sep = ":"))

    b = eta_new + b_noise
    colnames(b) = colnames(eta_new)

    y_new = matrix(nrow = n, ncol = n_y, dimnames = list(NULL, y_nm))
    for (k in seq_len(n_y)) {
        idx = seq_len(n_x) + (k - 1L) * n_x
        y_new[, k] = rowSums(b[, idx, drop = FALSE] * x)
    }

    est_true_vals = numeric(n_x * n_y)
    for (j in seq_len(n_x)) {
        for (k in seq_len(n_y)) {
            col_idx = j + (k - 1L) * n_x
            w_jk = x[, j] * total
            # Use eta_new (the expected b), not the drawn b. The EI estimator
            # targets E_b[beta_j | x, z] = sum x_j * eta_new_j * n / sum(n).
            # Using the drawn b would add b_cov noise that inflates "actual_bias."
            est_true_vals[col_idx] = sum(eta_new[, col_idx] * w_jk) / sum(w_jk)
        }
    }
    est_true = tibble::new_tibble(list(
        predictor = rep(x_nm, n_y),
        outcome = rep(y_nm, each = n_x),
        true = est_true_vals
    ))

    # --- Step 6: Package as ei_spec ---
    new_data = as.data.frame(spec)
    new_data[y_nm] = y_new
    new_data[".confounder"] = a
    ei_z = c(attr(spec, "ei_z"), if (isTRUE(include)) ".confounder")

    new_tibble(
        new_data,
        ei_x = attr(spec, "ei_x"),
        ei_y = y_nm,
        ei_z = ei_z,
        ei_n = total,
        ei_preproc = attr(spec, "ei_preproc"),
        ei_z_proc = attr(spec, "ei_z_proc"),
        b = b,
        b_loc = eta,
        b_cov = b_cov,
        eta = eta_new,
        est_true = est_true,
        c_outcome = c_outcome,
        c_predictor = c_predictor,
        confounding = confounding,
        class = "ei_spec"
    )
}