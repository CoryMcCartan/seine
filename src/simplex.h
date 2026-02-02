// (c) 2025 Cory McCARTAN

#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/matrix.hpp>
#include <cpp11/list.hpp>
#include "armadillo.hpp"

namespace simplex {

// Helper structure for simplex output
struct SimplexResult {
    arma::vec soln;
    int solved;        // 1 = solved, 0 = max iterations, -1 = stage 1
    double value;
    arma::mat constr_mat;
    arma::vec obj_coef;
    arma::uvec basic;
    
    // Stage 1 specific
    double val_aux;
    arma::vec obj_aux;
    
    SimplexResult() : solved(0), value(0.0), val_aux(0.0) {}
};

// Helper functions
arma::mat zero_mat(int n, int m);
arma::mat identity_mat(int n);
arma::mat pivot(const arma::mat& tab, int pivot_row, int pivot_col);

// Core simplex solver (original simplex1)
SimplexResult simplex_core(
    const arma::vec& obj_coef,
    const arma::mat& constr_mat,
    const arma::vec& rhs,
    const arma::vec& init,
    arma::uvec basic,
    double obj_val = 0.0,
    int stage = 2,
    int n_orig = -1,
    double eps = 1e-10,
    int max_iter = -1
);

// Main simplex wrapper
SimplexResult simplex(
    const arma::vec& obj_coef,
    const arma::mat& A1,
    const arma::vec& b1,
    const arma::mat& A2,
    const arma::vec& b2,
    const arma::mat& A3,
    const arma::vec& b3,
    bool maximize = false,
    int max_iter = -1,
    double eps = 1e-10
);

} // namespace simplex

#endif
