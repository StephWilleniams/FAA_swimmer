#ifndef functions_rng
#define functions_rng

#include "header.hpp"

// Manually set the seed value
////unsigned int seed = 1;

mt19937 generate_seed(int seed){
// Create a random number generator with the manual seed in main
    mt19937 gen(seed);
    return gen;
}

double nrand(double mean, double stddev, mt19937& gen) {
    normal_distribution<double> normal(mean, stddev);
    return normal(gen);
}

double urand(double a, double b, mt19937& gen) {
    uniform_real_distribution<double> uniform(a,b);
    return uniform(gen);
}

#endif // 