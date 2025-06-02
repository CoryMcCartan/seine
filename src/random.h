#ifndef RANDOM_H
#define RANDOM_H

#include "armadillo.hpp"

/*
 * Set RNG seed
 */
void seed_rng(int seed);

/*
 * Generate a uniform random double in [0, 1). Very slightly biased.
 */
double r_unif();

/*
 * Generate a standard normal deviate via Marsaglia polar method.
 */
double r_norm();

/*
 * Generate a random integer in [0, max) according to weights.
 */
int r_int_wgt(int max, const arma::vec& cum_wgts);
// helper
int find_u(double u, int max, const arma::vec& cum_wgts);


#endif
