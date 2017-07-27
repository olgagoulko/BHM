//Includes parameters and basic definitions needed for all subroutines

#ifndef BASIC_HPP
#define BASIC_HPP


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

//basic parameters and constants
#define PI 3.14159265358979323846
#define VERY_SMALL_NUMBER 1e-10
#define ACCURACY 1e-9
#define SVD_THRESHOLD 1e-12

#define defaultMinNumberTimesSampled 10

//aliases


//global parameters that may be passed as options
//extern double someGlobalParameter;

//convenience functions
double get_option(int inputN,char *inputV[],char *was);		//Passing options when starting program, e.g. as ./sampling -option1 5.0 -option2 1.0 ...
bool acceptreject(double probability, gsl_rng* RNG);		//whether to accept or reject update according to probability
bool isAround(double whatWeHave, double whatItShouldBe, double accuracy = ACCURACY);
int whatsign(double a);
int rounding(double a);
long unsigned int rdtsc();					//generates RNG seed by measuring the total pseudo-cycles since the processor was powered on


#endif
