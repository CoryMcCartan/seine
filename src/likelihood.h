#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include "armadillo.hpp"

/**
 * Log-likelihood of an intercept-only truncated normal EI model
 *
 * The parameters are `eta` is the global mean, and `L` the Cholesky factor of
 * the global covariance.
 *
 * All dimensions are assumed to be correct and are not checked. `tol` controls
 * the tolerance of the projection (see `proj_mvn()`).
 */
double llik(const arma::vec& eta, const arma::mat& L,
                    const arma::vec& y, const arma::mat& X,
                    const arma::vec& weights, double tol);

/**
 * Draw local parameters from their projected truncated normal distributions.
 *
 * Returns a column-major 3d array (first index varies fastest) of dimension
 * (draws, y.n_elem, eta.n_elem).
 */
arma::vec draw_local(int draws, const arma::vec& eta, const arma::mat& L,
                     const arma::vec& y, const arma::mat& X, int warmup, double tol);

/**
 * Project the MVN(eta, LL') onto the subspace of vectors b satisfying eps=x'(b-eta)
 *
 * Argument `Lx` is a vector that will be overwritten; used to avoid memory
 * allocation inside this funciton. It should be the same size as `eta`/`x`.
 * Argument `L_out` is a matrix that will be overwritten with the output; used
 * to allow memory reuse in loops that call this function. It should be square
 * of the same dimension as the length of `eta`/`x`.
 *
 * Returns a matrix where the final column is the projected mean and the other
 * columns are the projected Cholesky factor, which has one less column because
 * the projection operation lowers the rank by 1.
 */
void proj_mvn(const arma::vec& eta, const arma::mat& L, const arma::vec& x,
              double eps, arma::vec& Lx, arma::mat& L_out, double tol = 1e-7);

#endif
