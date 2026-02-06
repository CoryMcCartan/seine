#include <cassert>
#include "cpp11.hpp"
#include "cpp11armadillo.hpp"

#include "random.h"
#include "tmvn.h"
#include "epmgp.h"
#include "bounds.h"
#include "lp.h"

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
list R_bounds_lp(const doubles_matrix<>& x, const doubles_matrix<>& y, const doubles& bounds) {
    mat _x = as_Mat(x);
    mat _y = as_Mat(y);
    vec _bounds = as_Col(bounds);

    auto [min_mat, max_mat] = bounds_lp(_x, _y, _bounds);

    return writable::list({
        "min"_nm = as_doubles_matrix(min_mat),
        "max"_nm = as_doubles_matrix(max_mat)
    });
}

[[cpp11::register]]
list R_bounds_lp_contrast(const doubles_matrix<>& x, const doubles_matrix<>& y, 
                          const doubles_matrix<>& contr_m, double ub, 
                          double scale, double shift, bool sum_one, bool has_ub) {
    mat _x = as_Mat(x);
    mat _y = as_Mat(y);
    mat _contr_m = as_Mat(contr_m);

    auto [min_mat, max_mat] = bounds_lp_contrast_cpp(_x, _y, _contr_m, ub, 
                                                       scale, shift, sum_one, has_ub);

    return writable::list({
        "min"_nm = as_doubles_matrix(min_mat),
        "max"_nm = as_doubles_matrix(max_mat)
    });
}