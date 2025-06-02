// (c) 2025 Cory McCARTAN

#include <cmath>
#include "armadillo.hpp"

#include "random.h"

#ifndef M_2PI
const double M_2PI = 6.28318530717959;
#endif

using namespace arma;

// helper functions
void ell_calc_active_ints(const mat &A, const vec &u, subview_col<double> zt, const vec &nu,
                          mat &ints, int &n_int, vec &p, vec &q, vec &ur, mat &ab);
inline double mod_2pi(double ab);


mat ess_tmvn(int N, const vec &mu, const mat &L, const vec &init) {
    int n = L.n_cols;
    // Have x ~ N(mu, LL') constrained to [0, 1]^d
    // Equivalent to z ~ N(0, I); A z <= u, where A and u come from this function
    mat A = join_vert(-L, L);
    vec u = join_vert(mu, 1 - mu);
    mat z(n, N);
    // solve(L, init - mu) but handle rank-deficient case while using triangularity of L
    z.col(0) = solve(trimatl(L.head_rows(n)), init.head(n) - mu.head(n));

    vec nu(n);
    mat ints(2*A.n_rows + 1, 2);
    vec cum_lens(ints.n_rows);
    uvec rlookup(ints.n_rows);
    int n_int;
    double theta;
    // preallocate memory for `ell_calc_active_ints`
    vec p(A.n_rows), q(A.n_rows), ur(u.n_elem);
    mat ab(A.n_rows, 2);

    for (int i = 1; i < N; ++i) {
        nu.imbue(r_norm);

        ell_calc_active_ints(A, u, z.col(i-1), nu, ints, n_int, p, q, ur, ab);

        // find angle
        int i_int;
        double u = r_unif();
        if (n_int == 0) {
            theta = u * M_2PI;
        } else { // sample a point in the union of the intervals
            // build vector of cumulative interval lengths
            int n_active = 0;
            for (int j = 0; j < n_int; ++j) {
                if (ints(j, 0) < ints(j, 1)) {
                    double prev = n_active > 0 ? cum_lens(n_active - 1) : 0;
                    cum_lens(n_active) = prev + (ints(j, 1) - ints(j, 0));
                    rlookup(n_active) = j;
                    ++n_active;
                }
            }

            // sample an interval and re-use
            i_int = rlookup(r_int_wgt(n_active, cum_lens));
            theta = ints(i_int, 0) + u * (ints(i_int, 1) - ints(i_int, 0));
        }

        z.col(i) = z.col(i-1) * std::cos(theta) + nu * std::sin(theta);
    }

    z = (L * z); // can't call .each_col() on an expression like this; also, allocation needed
    return (z.each_col() + mu).t();
}

// helper to calculate active intervals
void ell_calc_active_ints(const mat &A, const vec &u, subview_col<double> zt, const vec &nu,
                          mat &ints, int &n_int, vec &p, vec &q, vec &ur, mat &ab) {
    int m = A.n_rows;
    p = A * zt; // last value of z is passed in for p
    q = A * nu;
    ur = u / arma::sqrt(square(p) + square(q));

    n_int = 0;
    for (int i = 0; i < m; ++i) {
        if (std::abs(ur(i)) <= 1) {
            double tau = std::atan2(q(i), p(i));
            double acos_ur = std::acos(ur(i));
            ab(n_int, 0) = mod_2pi(tau - acos_ur);
            ab(n_int, 1) = mod_2pi(tau + acos_ur);
            if (ab(n_int, 0) > ab(n_int, 1)) {
                std::swap(ab(n_int, 0), ab(n_int, 1));
            }
            ++n_int;
        }
    }

    const auto ab_h = ab.head_rows(n_int);
    uvec idx = sort_index(ab_h.col(0));
    ints(0, 0) = 0;
    ints(n_int, 1) = M_2PI;
    double cmax = 0.0;
    for (int i = 0; i < n_int; ++i) {
        cmax = std::max(cmax, ab_h(idx(i), 1));
        ints(i + 1, 0) = cmax;
        ints(i, 1) = ab_h(idx(i), 0);
    }
    ++n_int;
}

// wrap-around so output lies in [0, 2*pi]
inline double mod_2pi(double a) {
    const double result = std::fmod(a, M_2PI);
    return result >= 0 ? result : result + M_2PI;
}
