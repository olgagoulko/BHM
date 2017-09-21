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
