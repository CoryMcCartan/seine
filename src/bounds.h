#ifndef BOUNDS_H
#define BOUNDS_H

#include "armadillo.hpp"

/**
 * Compute bounds on entries of matrix B satisfying B * x = y
 *
 * For each row of x and y, solves LP to minimize/maximize each entry of B
 * subject to: B %*% x = y, B in [bounds[0], bounds[1]]
 *
 * B is a q x p matrix (q = ncol(y), p = ncol(x)) satisfying B %*% x[i,] = y[i,]
 * for row i. Each entry of the output matrices represents the optimal value
 * when optimizing that specific entry of B independently.
 *
 * @param x Matrix (n x p) where each row is a probability vector (nonnegative, sums to 1)
 * @param y Matrix (n x q) with same number of rows as x
 * @param bounds Vector of length 2: [lower, upper] for entries of B
 * @return A tuple of two matrices (min, max), each n x (q*p), with entries in row-major order
 */
std::tuple<arma::mat, arma::mat> bounds_lp(
    const arma::mat& x,
    const arma::mat& y,
    const arma::vec& bounds
);

#endif
