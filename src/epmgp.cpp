// (c) 2025 Cory McCARTAN

#include <cmath>
#include "armadillo.hpp"

#include "epmgp.h"

const double M_SQRT_2PI = 2.5066282746310002;

using namespace arma;

// helper functions
double log_erf_diff(double a, double b);

double utn_logZ(double mu, double sigma, double l, double u) {
    double alpha = (l - mu) / (M_SQRT2 * sigma);
    double beta = (u - mu) / (M_SQRT2 * sigma);
    return -M_LN2 + log_erf_diff(alpha, beta);
}

/*
 * Calculate moments for truncated univariate Normal given mean and standard deviation
 */
Moments utn_moments(double mu, double sigma2) {
    // fixing these so extensions of this code can make more general
    const double l = 0.0;
    const double u = 1.0;

    double sigma = std::sqrt(sigma2);
    double alpha = (l - mu) / (M_SQRT2 * sigma);
    double beta = (u - mu) / (M_SQRT2 * sigma);
    double log_m0 = -M_LN2 + log_erf_diff(alpha, beta);

    if (std::isnan(log_m0)) {
        return {0.0, mu > u ? u : l, 0.0};
    }

    double alpha_updf = std::exp(-alpha*alpha);
    double beta_updf = std::exp(-beta*beta);
    double m0_s2 = M_SQRT_2PI * std::exp(log_m0);
    double m1 = mu + sigma * (alpha_updf - beta_updf) / m0_s2;
    double m2 = mu*mu + sigma2 + sigma * (
        (l + mu)*alpha_updf - (u + mu)*beta_updf) / m0_s2 - m1*m1;

    return {log_m0, m1, m2};
}

// Implementation notes (see <https://arxiv.org/pdf/1111.6832>)
//
// The paper is wrong in claiming the Appendix B procedure reduces computational
// complexity from O(n^3) to O(n^2); it replaces a single O(n^3) update with n
// O(n^2) updates. The original update seems numerically stable enough and is
// implemented here.
//
// The derivatives were calculated by hand and checked numerically; the paper's
// version of these calculations ignores that many subexpressions are computed
// as part of the EP.
EPResult ep_moments(const vec& mu, const mat& L, const vec& addl_c, double addl_shift,
                    bool gr, double tol) {
    bool has_addl = addl_c.n_elem > 0;
    int n = mu.n_elem;
    int m = n + has_addl;
    // constants
    mat K_inv = trimatl(L).i().t() * trimatl(L).i();
    vec Kmu = K_inv * mu;
    mat CCt = addl_c * addl_c.t();
    // message parameters
    vec tau = zeros(m);
    vec nu = zeros(m);
    // cavity parameters
    vec cav_loc = zeros(m);
    vec cav_var = zeros(m);
    // moment estimates
    vec m1 = mu;
    vec prev_m1 = mu;
    mat m2_inv(n, n);
    mat m2 = L * L.t();
    vec m2_diag(m);
    vec m1_resc(m);
    double log_Z;

    for (int iter = 0; iter < 32; ++iter) {
        // form cavity locals
        m2_diag.head(n) = diagvec(m2);
        if (has_addl) m2_diag(n) = as_scalar(addl_c.t() * m2 * addl_c);
        cav_var = 1.0 / (1.0 / m2_diag - tau);
        m1_resc.head(n) = m1;
        if (has_addl) m1_resc(n) = dot(addl_c, m1);
        cav_loc = cav_var % (m1_resc / m2_diag - nu);

        log_Z = 0.0;
        for (int i = 0; i < m; ++i) {
            // Moment computation
            if (has_addl && i == n) cav_loc(i) += addl_shift;
            auto [log_Z_hat, mu_hat, sigma2_hat] = utn_moments(cav_loc(i), cav_var(i));
            if (has_addl && i == n) {
                cav_loc(i) -= addl_shift;
                mu_hat -= addl_shift;
            }

            log_Z += log_Z_hat;
            // match moments
            tau(i) = 1.0/sigma2_hat - 1.0/cav_var(i);
            nu(i) = mu_hat/sigma2_hat - cav_loc(i)/cav_var(i);
        }

        m2_inv = K_inv + diagmat(tau.head(n));
        if (has_addl) m2_inv += CCt * tau(n);
        m2 = inv_sympd(m2_inv);
        m1 = Kmu + nu.head(n);
        if (has_addl) m1 += nu(n) * addl_c;
        m1 = m2 * m1;

        if (vecnorm(m1 - prev_m1) < tol) break;
        prev_m1 = m1;
    }

    // eq. 59; we += for sum of Z_i
    // the first det terms is outside the -0.5(...) because
    // it is of a cholesky factor and so needs to be doubled
    log_Z += -sum(log(L.diag())) + 0.5 * (
        log_det_sympd(m2)
        -as_scalar(mu.t() * K_inv * mu) + as_scalar(m1.t() * m2_inv * m1) +
        sum(log1p(tau % cav_var)) +
        sum((square(cav_loc) % tau - 2 * cav_loc % nu - square(nu) % cav_var) /
            (1 + tau % cav_var))
    );

    // derivatives
    if (gr) {
        vec dlogZ_mu = K_inv * (m1 - mu);
        mat dlogZ_L = (dlogZ_mu * dlogZ_mu.t() + (K_inv * m2 * K_inv - K_inv)) * L;
        // mat dlogZ_L = trimatl((dlogZ_mu * dlogZ_mu.t() + (K_inv * m2 * K_inv - K_inv)) * trimatl(L));

        return {log_Z, m1, m2, dlogZ_mu, dlogZ_L};
    } else {
        return {log_Z, m1, m2, m1, m2};
    }
}

// default handler
EPResult ep_moments(const vec& mu, const mat& L, bool gr, double tol) {
    const vec addl_c;
    return ep_moments(mu, L, addl_c, 0.0, gr, tol);
}


/*
 * Numerically stable log(erf(b) - erf(a))
 */
double log_erf_diff(double a, double b) {
    // easy case
    if (std::signbit(a) != std::signbit(b)) {
        return std::log(std::erfc(a) - std::erfc(b));
    }
    // erfc stable for positive values
    if (a < 0 && b < 0) {
        double tmp = -a;
        a = -b;
        b = tmp;
    }

    double erfc_a = std::erfc(a);
    return std::log(erfc_a) + std::log1p(- std::erfc(b) / erfc_a);
}



#include <testthat.h>
context("EPMGP Tests") {
    test_that("Log erf difference is correctly calculated") {
        const double tol = 1.4901161e-8;
        double stdlib = std::log(std::erf(0.5) - std::erf(-0.5));
        double stable = log_erf_diff(-0.5, 0.5);
        expect_true(stdlib == stable);

        stdlib = std::log(std::erf(0.5) - std::erf(0.2));
        stable = log_erf_diff(0.2, 0.5);
        expect_true(abs(stdlib - stable) < tol);

        stdlib = std::log(std::erf(-6) - std::erf(-7));
        stable = log_erf_diff(-7, -6);
        expect_true(std::isinf(stdlib));
        expect_true(abs(stable - -38.3775592290446) < 1e-5);
    }

    test_that("Truncated Normal moments are correctly calculated") {
        const double tol = 1.4901161e-8;
        auto [log_m0a, m1a, m2a] = utn_moments(0.5, 0.0001);
        expect_true(abs(log_m0a - 0.0) < tol);
        expect_true(abs(m1a - 0.5) < tol);
        expect_true(abs(m2a - 0.0001) < tol);

        auto [log_m0b, m1b, m2b] = utn_moments(0.5, 1);
        expect_true(abs(log_m0b - -0.959916333695623) < tol);
        expect_true(abs(m1b - 0.5) < tol);

        auto [log_m0c, m1c, m2c] = utn_moments(8, 1);
        expect_true(abs(log_m0c - -27.3847937007199) < 1e-4);
    }

    test_that("Gaussian moments are correctly calculated") {
        const double tol = 1.4901161e-8;
        vec mu = {0.5, 0.5};
        mat L_small = 0.01 * eye(2, 2);
        auto res = ep_moments(mu, L_small);
        expect_true(std::abs(res.log_Z) < tol);
        expect_true(vecnorm(res.m1 - mu) < tol);

        mat L_med = eye(2, 2);
        res = ep_moments(mu, L_med);
        // being in bounds in each dimension is independent,
        // so this should be square of the 1D case
        expect_true(std::abs(res.log_Z - 2*std::log(0.382924922548026)) < tol);
        expect_true(vecnorm(res.m1 - mu) < tol);

        int n = 20;
        vec mu_long = 0.5 * ones(n);
        mat L_long = 0.01 * eye(n, n);
        res = ep_moments(mu_long, L_long);
        expect_true(abs(res.log_Z) < tol);

        vec mu_corner = {1.0, 1.0, 1.0};
        res = ep_moments(mu_corner, 0.01 * eye(3, 3));
        expect_true(std::abs(res.log_Z - std::log(0.125)) < tol);
        expect_true(res.m1(0) < mu_corner(0));
        expect_true(res.m1(1) < mu_corner(1));
        expect_true(res.m1(2) < mu_corner(2));
    }
}
