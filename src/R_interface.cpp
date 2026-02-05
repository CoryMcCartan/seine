#include <cassert>
#include "cpp11.hpp"
#include "cpp11armadillo.hpp"

#include "random.h"
#include "tmvn.h"
#include "epmgp.h"
#include "simplex.h"
#include "bounds.h"

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

// Simplex solver
[[cpp11::register]]
list simplex_cpp(const doubles& a,
                 SEXP A1, SEXP b1,
                 SEXP A2, SEXP b2,
                 SEXP A3, SEXP b3,
                 bool maxi = false,
                 int n_iter = -1,
                 double eps = 1e-10) {

    vec obj_coef = as_Col(a);
    int n = obj_coef.n_elem;

    // Convert NULL inputs to empty matrices/vectors
    mat _A1, _A2, _A3;
    vec _b1, _b2, _b3;
    int m1 = 0, m2 = 0, m3 = 0;

    if (!Rf_isNull(A1)) {
        _A1 = as_Mat(as_cpp<doubles_matrix<>>(A1));
        m1 = _A1.n_rows;
        _b1 = as_Col(as_cpp<doubles>(b1));
    }
    if (!Rf_isNull(A2)) {
        _A2 = as_Mat(as_cpp<doubles_matrix<>>(A2));
        m2 = _A2.n_rows;
        _b2 = as_Col(as_cpp<doubles>(b2));
    }
    if (!Rf_isNull(A3)) {
        _A3 = as_Mat(as_cpp<doubles_matrix<>>(A3));
        m3 = _A3.n_rows;
        _b3 = as_Col(as_cpp<doubles>(b3));
    }

    int m = m1 + m2 + m3;
    if (n_iter == -1) {
        n_iter = n + 2 * m;
    }

    // Call C++ simplex
    simplex::SimplexResult result = simplex::simplex(
        obj_coef, _A1, _b1, _A2, _b2, _A3, _b3, maxi, n_iter, eps
    );

    // Build return list matching R structure
    writable::list out;

    // Extract solution (first n elements only)
    vec solution = result.soln.subvec(0, n - 1);
    out.push_back({"soln"_nm = as_doubles(solution)});
    out.push_back({"solved"_nm = result.solved});
    out.push_back({"value"_nm = result.value});
    out.push_back({"maxi"_nm = maxi});

    // Add slack variables if m1 > 0
    if (m1 > 0) {
        vec slack = result.soln.subvec(n, n + m1 - 1);
        out.push_back({"slack"_nm = as_doubles(slack)});
    }

    // Add surplus variables if m2 > 0
    if (m2 > 0) {
        vec surplus = result.soln.subvec(n + m1, n + m1 + m2 - 1);
        out.push_back({"surplus"_nm = as_doubles(surplus)});
    }

    // Add artificial variables if in stage 1
    if (result.solved == -1) {
        vec artificial = result.soln.subvec(n + m1 + m2, result.soln.n_elem - 1);
        out.push_back({"artificial"_nm = as_doubles(artificial)});
    }

    // Add objective coefficients
    out.push_back({"obj"_nm = as_doubles(obj_coef)});

    return out;
}