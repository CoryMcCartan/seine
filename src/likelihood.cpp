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
double mvnorm_lupdf(const vec& x, const vec& mu, const mat& L) {
    vec diff = solve(trimatl(L), x - mu);
    return -sum(log(diagvec(L))) - 0.5 * dot(diff, diff);
}

// integral over tomography line, 2x2 case
double llik_S_2x2(const mat& L_eta_proj, double y, const mat& X, int i) {
    // pick variable to param. tomog. line: 0 by default, but 1 if at boundary
    int j = X(i, 0) == 0 || X(i, 0) == 1;
    // bounds and sd in tomog. line parametrization
    double lb = (y - X(i, 1 - j)) / X(i, j);
    double ub = y / X(i, j);
    double sd;
    if (j == 0) {
        // sd of 1st component is just L[0, 0]
        sd = L_eta_proj(j, 0);
    } else {
        // var. of 2nd component is L[0, 0]^2 + L[1, 0]^2
        sd = std::sqrt(sum(square(L_eta_proj.col(0))));
    }
    return utn_logZ(L_eta_proj(j, 1), sd, lb < 0 ? 0 : lb, ub > 1 ? 1 : ub);
}

// integral over tomography hyperplane
double llik_S_ep(const mat& L_eta_proj, double y, const mat& X, int i,
                  int k, vec& eta_ep, mat& L_ep, vec& c_ep, double tol) {
    // find last nonzero entry of x
    int last_nz;
    for (last_nz = k - 1; last_nz >= 0; last_nz--) {
      if (X(i, last_nz) != 0) {
          break;
      }
    }

    // set up EP inputs
    for (int j = 0; j < k - 1; ++j) {
        int js = j + (j >= last_nz); // shifted
        eta_ep(j) = L_eta_proj(js, k - 1);
        L_ep.row(j) = L_eta_proj.row(js).head(k - 1);
        c_ep(j) = -X(i, js) / X(i, last_nz);
    }

    return ep_moments(eta_ep, L_ep, c_ep, y / X(i, last_nz), false, tol).log_Z;
}

double llik(const mat& eta, const mat& L, const vec& y, const mat& X,
            const vec& weights, int p, double tol) {
    int n = y.n_elem;
    int k = eta.n_cols;
    vec eps_loc = y - sum(X % eta, 1); // row sums
    vec var_loc = diagvec(X * trimatl(L) * trimatl(L).t() * X.t());
    double L_det = sum(log(diagvec(L)));

    double llik = 0;
    if (p == 1) { // when no covariates, R term shared across all obs
        llik += -sum(weights) * ep_moments(eta.row(0).t(), L, false, tol).log_Z;
    }

    mat L_eta_proj(k, k);
    vec diff_unam(k);
    vec Lx(k);
    vec eta_ep(k - 1);
    mat L_ep(k - 1, k - 1);
    vec c_ep(k - 1);
    for (int i = 0; i < n; ++i) {
        // 1 / R(eta, L); false=no gradient
        double llik_i = 0.0;
        if (p > 1) {
            llik_i += -ep_moments(eta.row(i).t(), L, false, tol).log_Z;
        }

        // unanimous units
        // calculate MVN density at (0, 0, ... 0) and (1, 1, ... 1)
        if (y[i] == 0) {
            diff_unam = solve(trimatl(L), -eta.row(i).t());
            llik_i += -L_det - 0.5 * dot(diff_unam, diff_unam);
        } else if (y[i] == 1) {
            diff_unam = solve(trimatl(L), 1 - eta.row(i).t());
            llik_i += -L_det - 0.5 * dot(diff_unam, diff_unam);
        } else {
            proj_mvn(eta.row(i).t(), L, X.row(i).t(), eps_loc(i), Lx, L_eta_proj, tol);

            // log S
            if (k == 2) { // analytical
                llik_i += llik_S_2x2(L_eta_proj, y[i], X, i);
            } else { // EP
                llik_i += llik_S_ep(L_eta_proj, y[i], X, i, k, eta_ep, L_ep, c_ep, tol);
            }

            // y ~ N(x'eta, x'Sigma x)
            llik_i += norm_lupdf(eps_loc(i), var_loc(i));
        }

        llik += weights(i) * llik_i;
    }

    return llik;
}

vec draw_local(int draws, const mat& eta, const mat& L,
               const vec& y, const mat& X, int warmup, double tol) {
    int n = y.n_elem;
    int k = eta.n_cols;
    vec eps_loc = y - sum(X % eta, 1); // row sums
    vec out(draws * n * k);

    mat L_eta_proj(k, k);
    vec Lx(k);
    vec init(k);
    for (int i = 0; i < n; ++i) {
        proj_mvn(eta.row(i).t(), L, X.row(i).t(), eps_loc(i), Lx, L_eta_proj, tol);

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
            for (int j = n - 1; j >= 0; --j) {
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
    double tol = 1e-5;

    test_that("Normal distributions calculated correctly") {
        // dnorm(0.0, 0, 0.3, log=TRUE) - dnorm(0.5, 0, 0.3, log=TRUE)
        double diff_r = 1.3888888888889;
        double diff_act = norm_lupdf(0.0, 0.09) - norm_lupdf(0.5, 0.09);
        expect_true(abs(diff_r - diff_act) < tol);

        vec x = {0.5, 1.5};
        vec mu = {0.0, 1.0};
        mat L = eye(2, 2) * 0.3;
        double d_un = norm_lupdf(0.5, 0.09);
        double d_mv = mvnorm_lupdf(x, mu, L);
        expect_true(abs(2*d_un - d_mv) < tol);
    }

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
