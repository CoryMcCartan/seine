#ifndef EPMGP_H
#define EPMGP_H

#include "armadillo.hpp"

/**
 * Container for moments of a truncated univariate Normal distribution.
 */
struct Moments {
    double log_Z; /** log of the normalizing constant */
    double m1; /** mean */
    double m2; /** variance */
};

/**
 * Container for moments of a truncated multivariate Normal distribution.
 */
struct EPResult {
    double log_Z; /** log of the normalizing constant */
    arma::vec m1; /** mean */
    arma::mat m2; /** Covariance */
    arma::vec dlogZ_mu; /** derivative of log_Z wrt mean */
    arma::mat dlogZ_L; /** derivative of log_Z wrt Cholesky of covariance */
};

/**
 * Estimate moments of multivariate Normal distribution truncated to unit hypercube
 *
 * The MVN has mean `mu` and Cholesky decomposition of covariance `L`.
 * Optionally, an additional truncation constraint can be described by providing
 * a vector `addl_c` and a constraint value `addl_shift`. If not used, `addl_c`
 * should have length 0.
 *
 * Implements:
 *   Cunningham, J. P., Hennig, P., & Lacoste-Julien, S. (2011). Gaussian
 *   probabilities and expectation propagation. arXiv preprint arXiv:1111.6832.
 */
EPResult ep_moments(const arma::vec& mu, const arma::mat& L, bool gr=true, double tol=1e-7);
// with additional constraint
EPResult ep_moments(const arma::vec& mu, const arma::mat& L, const arma::vec& addl_c,
                    double addl_shift, bool gr=true, double tol=1e-7);


/**
 * Exact moments of univariate Normal distribution truncated to [0, 1]
 * `sigma2` here is the variance, not the standard deviation.
 */
Moments utn_moments(double mu, double sigma2);

/**
 * Exact log normalizing constant univariate Normal distribution truncated to [l, u]
 * `sigma` here is the standard deviation, not the variance.
 */
double utn_logZ(double mu, double sigma, double l, double u);


#endif
