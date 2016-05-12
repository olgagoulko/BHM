#include "basic.hpp"

using namespace std;

int main(int argc, char **argv) {
    cout << "Hello, world!" << std::endl;
    
    //Random Number generator
    long unsigned int seed=rdtsc();
    gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
    gsl_rng_set (RNG, seed);
    double randy=gsl_rng_uniform(RNG);
    cout << "My random number is: " << randy << endl;
    
    return 0;
}
