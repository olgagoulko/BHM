/** @file print_spline_grid.hpp
    Declares convenience function to print spline grid.
*/

#ifndef BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d
#define BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d

#include <iosfwd>
#include "spline.hpp"

std::ostream& print_spline_grid(std::ostream&, const splineArray&, unsigned int nrows=1024);

#endif /* BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d */

