/** @file poly.cpp
    Implementation of polynomial function from coefficients.
**/


#include "poly.hpp"

double Poly::operator()(double x) const
{
    double y=0;
    // Straightforward Horner's method, https://en.wikipedia.org/wiki/Horner%27s_method#C_implementation
    for (container_type_::const_reverse_iterator cit=coeff_.rbegin();
         cit!=coeff_.rend(); ++cit) {
        y *= x;
        y += *cit;
    }
    return y;
}
