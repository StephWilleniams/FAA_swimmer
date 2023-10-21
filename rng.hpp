#ifndef rng
#define rng

#include "header.hpp"

// Manually set the seed value
unsigned int seed = 1;

// Create a random number generator with the manual seed in main
mt19937 gen(seed);

double nrand(double mean, double stddev, mt19937& gen) {
    normal_distribution<double> normal(mean, stddev);
    return normal(gen);
}

double urand(double a, double b, mt19937& ) {
    uniform_real_distribution<double> uniform(a,b);
    return uniform(gen);
}

#endif // 