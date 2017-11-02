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
/** @file print_spline_grid.hpp
    Declares convenience function to print spline grid.
*/

#ifndef BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d
#define BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d

#include <iosfwd>
#include "spline.hpp"

std::ostream& print_spline_grid(std::ostream&, const splineArray&, unsigned int nrows=1024);

#endif /* BHM_PRINT_SPLINE_GRID_HPP_INCLUDED_35a04bc33ac14ad28a24c217cd83025d */

