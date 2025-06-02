// (c) 2025 Cory McCARTAN

#include <cmath>

#include "armadillo.hpp"

#include "likelihood.h"
#include "epmgp.h"
#include "tmvn.h"

using namespace arma;

// log unnormlized Normal density
inline double norm_lupdf(double x, double var) {
    return -0.5 * (x * x / var + std::log(var));
}


double llik_intonly(const vec& eta, const mat& L, const vec& y, const mat& X,
                    const vec& weights, double tol) {
    int k = eta.n_elem;
    vec eps_loc = y - X * eta;
    vec var_loc = diagvec(X * trimatl(L) * trimatl(L).t() * X.t());

    // 1 / R(eta, L) [norm. const. for overall TMVN]; false=no gradient
    double llik = -sum(weights) * ep_moments(eta, L, false, tol).log_Z;

    mat L_eta_proj(k, k);
    vec Lx(k);
    vec eta_ep(k - 1);
    mat L_ep(k - 1, k - 1);
    vec c_ep(k - 1);
    for (int i = 0; i < y.n_elem; ++i) {
        proj_mvn(eta, L, X.row(i).t(), eps_loc(i), Lx, L_eta_proj, tol);

        double log_S;
        if (k == 2) { // analytical
            int j = X(i, 0) > 0 ? 0 : 1;
            double lb = (y(i) - X(i, 1 - j)) / X(i, j);
            double ub = y(i) / X(i, j);
            double sd = j == 0 ? L_eta_proj(j, 0) : std::sqrt(sum(square(L_eta_proj.col(0))));
            log_S = utn_logZ(L_eta_proj(j, 1), sd, lb < 0 ? 0 : lb, ub > 1 ? 1 : ub);
        } else { // EP
            // find last nonzero entry of x
            int last_nz;
            for (last_nz = k - 1; last_nz >= 0; last_nz--) if (X(i, last_nz) != 0) break;

            // set up EP inputs
            for (int j = 0; j < k - 1; ++j) {
                int js = j + (j >= last_nz); // shifted
                eta_ep(j) = L_eta_proj(js, k - 1);
                L_ep.row(j) = L_eta_proj.row(js).head(k - 1);
                c_ep(j) = -X(i, js) / X(i, last_nz);
            }

            log_S = ep_moments(eta_ep, L_ep, c_ep, y(i)/X(i, last_nz), false, tol).log_Z;
        }
        // 2nd term is  y ~ N(x'eta, x'Sigma x)
        llik += weights(i) * (log_S + norm_lupdf(eps_loc(i), var_loc(i)));
    }

    return llik;
}

vec draw_local(int draws, const vec& eta, const mat& L,
               const vec& y, const mat& X, int warmup, double tol) {
    int n = y.n_elem;
    int k = eta.n_elem;
    vec eps_loc = y - X * eta;
    vec out(draws * n * k);

    mat L_eta_proj(k, k);
    vec Lx(k);
    vec init(k);
    for (int i = 0; i < n; ++i) {
        proj_mvn(eta, L, X.row(i).t(), eps_loc(i), Lx, L_eta_proj, tol);

        init.fill(y(i));
        mat draws_i = ess_tmvn(draws + warmup, L_eta_proj.col(k - 1),
                               L_eta_proj.head_cols(k - 1), init).tail_rows(draws);

        // strided copy :(
        for (int l = 0; l < k; ++l) {
            int first = draws*(i + n*l);
            out.subvec(first, first + draws - 1) = draws_i.col(l);
        }
    }

    return out;
}


void proj_mvn(const vec& eta, const mat& L, const vec& x, double eps, vec& Lx, mat& L_out, double tol) {
    int n = eta.n_elem;
    Lx = L * L.t() * x;
    double xLx = as_scalar(x.t() * Lx);

    L_out = L;
    L_out.col(n - 1) = eta + (eps / xLx) * Lx;

    for (int i = 0; i < n - 1; ++i) {
        double r = std::sqrt(L_out(i, i)*L_out(i, i) - Lx(i)*Lx(i) / xLx);
        if (r > tol) {
            L_out(span(i+1, n-1), i) = (L_out(i, i) * L_out(span(i+1, n-1), i)
                                         - Lx(i) * Lx(span(i+1, n-1)) / xLx) / r;
            Lx(span(i+1, n-1)) = (r * Lx(span(i+1, n-1))
                                      - Lx(i) * L_out(span(i+1, n-1), i)) / L_out(i, i);
        } else {
            r = 0.0;
            L_out(span(i+1, n-1), i).zeros();
            Lx(span(i+1, n-1)) = -Lx(i) * L_out(span(i+1, n-1), i) / L_out(i, i);
            // figure out last zero in x, which needs its own entry
            for (int j = n - 1; j >= 0; ++j) {
                if (x(j) == 0) {
                    L_out(j, i) = L(j, j);
                    break;
                }
            }
        }
        L_out(i, i) = r;
    }
}



#include <testthat.h>
#include "random.h"
context("Likelihood Tests") {
    double tol = 1e-8;

    test_that("MVN subspace projections calculated correctly") {
        int k = 5;
        mat Sigma = 3 + eye(k, k);
        mat L = chol(Sigma, "lower");
        vec eta = 0.5 * ones(k);
        vec Lx(k);
        mat L2(k, k);

        proj_mvn(eta, L, ones(k)/k, 0.0, Lx, L2);
        expect_true(max(abs(L2.col(k - 1) - eta)) < tol); // eps = 0 should give back eta

        auto err = [&](const vec &x) {
            mat Sigma2 = Sigma - (Sigma * x * x.t() * Sigma) / as_scalar(x.t() * Sigma * x);
            proj_mvn(eta, L, x, 0.0, Lx, L2);
            return max(abs(Sigma2 - L2.head_cols(k - 1) * L2.head_cols(k - 1).t()).as_col());
        };

        expect_true(err({1, 0, 0, 0, 0}) < tol);
        expect_true(err({0, 1, 0, 0, 0}) < tol);
        expect_true(err({0, 0, 1, 0, 0}) < tol);
        expect_true(err({0, 0, 0, 1, 0}) < tol);
        expect_true(err({0, 0, 0, 0, 1}) < tol);
        expect_true(err({0.5, 0.5, 0, 0, 0}) < tol);
        expect_true(err({0, 0.5, 0, 0.5, 0}) < tol);
        expect_true(err({0, 0.5, 0, 0, 0.5}) < tol);
        expect_true(err({0, 0.3, 0.4, 0.3, 0}) < tol);

        vec x(k);
        bool any = false;
        for (int i = 0; i < 1000; ++i) {
            // x ~ Dirichlet(1, 1, ... 1)
            x.imbue(r_unif);
            x = -log(x);
            x = x / sum(x);
            if (err(x) > 100 * tol) {
                any = true;
                break;
            }
        }
        expect_false(any);
    }
}
