/*** LICENCE: ***
Bin histogram method for restoration of smooth functions from noisy integrals. Copyright (C) 2017 Olga Goulko

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.

*** END OF LICENCE ***/
//Includes parameters and basic definitions needed for all subroutines

#ifndef BASIC_HPP_cdce02947f3f4772ad67a5420fdedf42
#define BASIC_HPP_cdce02947f3f4772ad67a5420fdedf42

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include <limits.h>
#include <float.h>
#include <cstdlib>
#include <cassert>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf.h>
#include <time.h>

// Our logger class
#include "logger.hpp"

//basic parameters and constants
#define PI 3.14159265358979323846
#define VERY_SMALL_NUMBER 1e-10
#define ACCURACY 1e-9
#define SVD_THRESHOLD 1e-12

#define defaultMinNumberTimesSampled 10

//convenience functions
double get_option(int inputN,char *inputV[],char *was);		//Passing options when starting program, e.g. as ./sampling -option1 5.0 
bool acceptreject(double probability, gsl_rng* RNG);		//whether to accept or reject update according to probability
bool isAround(double whatWeHave, double whatItShouldBe, double accuracy = ACCURACY);
int whatsign(double a);
int rounding(double a);
long unsigned int rdtsc();					//generates RNG seed by measuring the total pseudo-cycles since the processor was powered on

/// verbosity level
enum verbosity_level_type { CONCISE=0, VERBOSE=1 };

/// return codes for OS
enum return_code_type {
    OK        = 0, ///< No problems
    BAD_ARGS  = 1, ///< Invalid/missing arguments to the program
    BAD_DATA  = 2, ///< Invalid data supplied
    ZERO_DATA = 3, ///< Data is consistent with zero, and the stop was requested
    BAD_FIT   = 4, ///< Failed to find good fit (the result may still be output)
    OTHER_ERROR =100 ///< Some other error
};


/// if n>0 && (n == 2^k), then return k, otherwise return -1.
inline int ilog2(int n) {
    // ref: https://en.wikipedia.org/wiki/Power_of_two#Fast_algorithm_to_check_if_a_positive_number_is_a_power_of_two
    if ((n<=0) || ((n & (n-1)) != 0)) return -1;
    int m=0;
    while ((n>>=1)) m++;
    return m;
}


#endif /* BASIC_HPP_cdce02947f3f4772ad67a5420fdedf42 */
