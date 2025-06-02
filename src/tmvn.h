#ifndef TMVN_H
#define TMVN_H

#include "armadillo.hpp"

/**
 * Elliptical Slice Sampler for multivariate normal distribution truncated to
 * the unit hypercube
 *
 * Implements the following algorithm:
 *   Wu, Kaiwen, and Jacob R. Gardner. "A Fast, Robust Elliptical Slice Sampling
 *   Implementation for Linearly Truncated Multivariate Normal Distributions."
 *   arXiv preprint arXiv:2407.10449 (2024).
 *
 * Seems to converge extremely quickly (5 steps or so).
 *
 * @param N the number of samples
 * @param mu the mean of the distribution
 * @param L the Cholesky decomposition of the covariance matrix
 * @param init the initial value for the sampler
 */
arma::mat ess_tmvn(int N, const arma::vec &mu, const arma::mat &L, const arma::vec &init);

#endif
