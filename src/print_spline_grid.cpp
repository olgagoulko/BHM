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
/** @file print_spline_grid.cpp
    Implements convenience function to print spline grid.
*/

#include "print_spline_grid.hpp"
#include "grid.hpp"

// FIXME: this can be done much simpler using boost::function (from Boost) or std::function (from C++11)
struct callable {
    virtual double operator()(double x) const =0;
    virtual ~callable() {};
};

struct spline_value : public callable {
    const splineArray& spline_;

    spline_value(const splineArray& s): spline_(s) {}
    virtual double operator()(double x) const { return spline_.splineValue(x); }
};

struct spline_error : public callable {
    const splineArray& spline_;

    spline_error(const splineArray& s): spline_(s) {}
    virtual double operator()(double x) const { return spline_.splineError(x); }
};


std::ostream& print_spline_grid(std::ostream& strm, const splineArray& aspline, unsigned int nrows)
{
    Grid<callable> grid;
    grid.set_rows(nrows)
        .set_xmin(aspline.getLowerBound())
        .set_xmax(aspline.getUpperBound());
    
    spline_value s_val(aspline);
    grid.add_column(&s_val);
    
    spline_error s_err(aspline);
    grid.add_column(&s_err);

    strm << grid;
    return strm;
}
