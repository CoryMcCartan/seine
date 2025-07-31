#include <cassert>
#include "cpp11.hpp"
#include "cpp11armadillo.hpp"

#include "random.h"
#include "tmvn.h"
#include "epmgp.h"
#include "likelihood.h"

using namespace arma;
using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
void R_sync_rng() {
    GetRNGstate();
    int seed = R_unif_index(1000000);
    PutRNGstate();
    seed_rng(seed);
}

[[cpp11::register]]
doubles_matrix<> R_ess_tmvn(int N, const doubles& mu, const doubles_matrix<>& L, const doubles& init) {
    vec _mu = as_Col(mu);
    mat _L = as_Mat(L);
    vec _init = as_Col(init);
    return as_doubles_matrix(ess_tmvn(N, _mu, _L, _init));
}

[[cpp11::register]]
list R_ep_moments(const doubles& mu, const doubles_matrix<>& L,
                  const doubles& addl_c, double addl_shift, double tol) {
    vec _mu = as_Col(mu);
    mat _L = as_Mat(L);
    vec _addl_c = as_Col(addl_c);
    EPResult res = ep_moments(_mu, _L, _addl_c, addl_shift, 32, tol);

    return writable::list({
        as_sexp(res.log_Z), as_doubles(res.m1), as_doubles_matrix(res.m2),
        as_doubles(res.dlogZ_mu), as_doubles_matrix(res.dlogZ_L)
    });
}

[[cpp11::register]]
doubles R_utn_moments(double mu, double sigma2) {
    writable::doubles out(3);
    auto [m0, m1, m2] = utn_moments(mu, sigma2);
    out[0] = m0;
    out[1] = m1;
    out[2] = m2;
    return out;
}

[[cpp11::register]]
double R_llik_intonly(const doubles& eta, const doubles_matrix<>& L, const doubles& y,
                      const doubles_matrix<>& X, const doubles& weights, double tol) {
    vec _eta = as_Col(eta);
    mat _L = as_Mat(L);
    vec _y = as_Col(y);
    mat _X = as_Mat(X);
    vec _weights = as_Col(weights);
    return llik_intonly(_eta, _L, _y, _X, _weights, tol);
}

[[cpp11::register]]
doubles R_draw_local(int draws, const doubles& eta, const doubles_matrix<>& L,
                     const doubles& y, const doubles_matrix<>& X, int warmup, double tol) {
    vec _eta = as_Col(eta);
    mat _L = as_Mat(L);
    vec _y = as_Col(y);
    mat _X = as_Mat(X);
    return as_doubles(draw_local(draws, _eta, _L, _y, _X, warmup, tol));
}

[[cpp11::register]]
doubles_matrix<> r_proj_mvn(const doubles& eta, const doubles_matrix<>& l, const doubles& x, double eps) {
    vec _eta = as_col(eta);
    mat _l = as_mat(l);
    vec _x = as_col(x);
    vec _lx(_x.n_elem);
    mat _l_out(_x.n_elem, _x.n_elem);

    proj_mvn(_eta, _l, _x, eps, _lx, _l_out);
    return as_doubles_matrix(_l_out);
}

// Returns list; 1st entry is projected means, 2nd is projected variances
[[cpp11::register]]
list r_proj_local(const doubles_matrix<>& eta, const doubles_matrix<>& e_cov,
                  const doubles_matrix<>& r_cov,
                  const doubles_matrix<>& x, const doubles eps) {
    int n = eta.nrow();
    int c = x.ncol();
    vec _eta(c);
    mat _l(c, c);
    vec _x(c);
    vec _lx(c);
    mat _l_out(c, c);
    writable::doubles_matrix<> eta_proj(n, c);
    writable::doubles_matrix<> var_proj(n, c);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < c; j++) {
            _eta(j) = eta(i, j);
            _x(j) = x(i, j);
            var_proj(i, j) = 0;
        }

        proj_mvn(_eta, _l, _x, eps[i], _lx, _l_out);

        for (int j = 0; j < c; j++) {
            eta_proj(i, j) = _l_out(j, c - 1);
            for (int k = 0; k < c - 1; k++) {
                var_proj(i, j) += _l_out(j, k)*_l_out(j, k);
            }
        }
    }

    return writable::list({eta_proj, var_proj});
}
