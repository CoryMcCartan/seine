#ifndef LP_H
#define LP_H

#include "armadillo.hpp"
#include <tuple>

/**
 * Status codes for LP solver
 */
enum class LPStatus {
    OPTIMAL = 0,      // Optimal solution found
    INFEASIBLE = 1,   // Problem is infeasible
    UNBOUNDED = 2,    // Problem is unbounded
    ERROR = 3         // Numerical or other error
};

/**
 * Result of LP solver
 */
struct LPResult {
    LPStatus status;
    double objval;           // Optimal objective value (if status == OPTIMAL)
    arma::vec solution;      // Optimal solution (if status == OPTIMAL)
    
    LPResult() : status(LPStatus::ERROR), objval(0.0) {}
};

/**
 * Solve a linear program in standard form using two-phase simplex method
 * 
 * Solves: minimize c'x subject to Ax = b, x >= 0
 * 
 * @param c Objective coefficients (n_vars x 1)
 * @param A Constraint matrix (n_constraints x n_vars)
 * @param b Right-hand side (n_constraints x 1)
 * @param tol Numerical tolerance for zero comparisons (default: 1e-10)
 * @return LPResult containing status, objective value, and solution
 */
LPResult solve_lp_simplex(
    const arma::vec& c,
    const arma::mat& A,
    const arma::vec& b,
    double tol = 1e-10
);

/**
 * Compute bounds on contrasts using custom LP solver
 * 
 * For each observation i and contrast j, solves LP to find min/max of:
 *   contr_m[,j]' * B
 * subject to:
 *   B %*% x[i,] = y[i,]  (n_y constraints)
 *   sum(B[,k]) = scale*(1-shift) for k=1..n_x if sum_one (n_x constraints)
 *   B <= ub if has_ub (n_vars constraints)
 *   B >= 0 (implicit after transformation)
 * 
 * where B is n_y x n_x matrix (vectorized to n_vars = n_y*n_x)
 * 
 * @param x Predictor matrix (n x n_x)
 * @param y Outcome matrix (n x n_y), already transformed: scale*(y - shift)
 * @param contr_m Contrast matrix (n_vars x n_c)
 * @param ub Upper bound (after transformation)
 * @param scale Scaling factor (+1 or -1)
 * @param shift Shift amount
 * @param sum_one Whether columns of B should sum to scale*(1-shift)
 * @param has_ub Whether upper bound constraints should be enforced
 * @return Tuple of (min_mat, max_mat) each n x n_c
 */
std::tuple<arma::mat, arma::mat> bounds_lp_contrast_cpp(
    const arma::mat& x,
    const arma::mat& y,
    const arma::mat& contr_m,
    double ub,
    double scale,
    double shift,
    bool sum_one,
    bool has_ub
);

#endif
