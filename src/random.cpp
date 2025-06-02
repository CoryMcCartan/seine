#include <cmath>
#include <vector>
#include <cstdint>
#include <random>

#include "armadillo.hpp"

std::random_device rd;


/* This is a fixed-increment version of Java 8's SplittableRandom generator
 See http://dx.doi.org/10.1145/2714064.2660195 and
 http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 Written in 2015 by Sebastiano Vigna (vigna@acm.org)
 [Public Domain]
 */
static uint64_t state_sr; /* The state can be seeded with any value. */
uint64_t next_sr() {
    uint64_t z = (state_sr += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}


/* This is xoshiro256+ 1.0
 Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 [Public domain]
 */
static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t state_xo[4] = {rd(), rd(), rd(), rd()};

uint64_t generator(void) {
    const uint64_t result = state_xo[0] + state_xo[3];

    const uint64_t t = state_xo[1] << 17;

    state_xo[2] ^= state_xo[0];
    state_xo[3] ^= state_xo[1];
    state_xo[1] ^= state_xo[2];
    state_xo[0] ^= state_xo[3];

    state_xo[2] ^= t;

    state_xo[3] = rotl(state_xo[3], 45);

    return result;
}
// Rest of file is original code --------------------------------


/*
 * Set RNG seed
 */
void seed_rng(int seed) {
    state_sr = seed;
    // seed xoshiro256+ with SplittableRandom, as recommended by authors
    state_xo[0] = next_sr();
    state_xo[1] = next_sr();
    state_xo[2] = next_sr();
    state_xo[3] = next_sr();
}


/*
 * Generate a uniform random double in [0, 1). Slightly biased.
 */
double r_unif() {
    return (generator() >> 11) * 0x1.0p-53;
}


static double _norm_z = 0.0;
static double _norm_valid = false;
/*
 * Generate a standard normal deviate via Marsaglia polar method.
 */
double r_norm() {
    if (_norm_valid) {
        _norm_valid = false;
        return _norm_z;
    }

    double u1, u2, s;
    do {
        u1 = r_unif()*2.0 - 1.0;
        u2 = r_unif()*2.0 - 1.0;
        s = u1*u1 + u2*u2;
    } while (s >= 1);

    double adj = std::sqrt(-2 * std::log(s) / s);
    _norm_z = u2 * adj;
    _norm_valid = true;
    return u1 * adj;
}

// helper
int find_u(double u, int max, const arma::vec& cum_wgts) {
    int low = 0, high = max - 1;
    double upper = cum_wgts[high];

    if (cum_wgts[0] > u*upper)
        return 0;

    while (high - low > 1) {
        int midpt = std::ceil((high + low) / 2.0);
        if (cum_wgts[midpt] <= u*upper)
            low = midpt;
        else
            high = midpt;
    }

    return high;
}

/*
 * Generate a random integer in [0, max) according to (cumulative) weights.
 */
int r_int_wgt(int max, const arma::vec& cum_wgts) {
    return find_u(r_unif(), max, cum_wgts);
}
