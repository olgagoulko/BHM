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
/** @file grid.hpp is a header for generating values of a function on a grid */

#ifndef BHM_GRID_HPP_INCLUDED_426f62d8250842fb9f11de8116019c61
#define BHM_GRID_HPP_INCLUDED_426f62d8250842fb9f11de8116019c61

#include <vector>
#include <iostream>
#include <stdexcept>

namespace detail {
    class uncopyable {
      private:
        uncopyable(const uncopyable&) {}
        uncopyable& operator=(const uncopyable&) { return *this; }
      protected:
        uncopyable() {}
        ~uncopyable() {}
    };
}

template <typename F>
class Grid : private detail::uncopyable {
  public:
    typedef F fn_type;
  private:
    typedef std::vector<const fn_type*> fn_vector;
    fn_vector columns_;

    double xmin_;
    double xmax_;
    unsigned int nrows_;

  public:
    Grid() : columns_(), xmin_(0), xmax_(0), nrows_(512) {}

    unsigned int columns() const { return columns_.size(); }
    unsigned int rows() const { return nrows_; }

    double xmax() const { return xmax_; }
    double xmin() const { return xmin_; }

    Grid& set_xmin(double xmin) { xmin_=xmin; return *this; }
    Grid& set_xmax(double xmax) { xmax_=xmax; return *this; }
    Grid& set_rows(unsigned int rows) { nrows_=rows; return *this; }
    
    /// Adds a column containing (borrowed!) pointer to a callable object
    /** @note The pointer is borrowed: this grid object
        must not outlive the pointed-to callable object. */
    Grid& add_column(const fn_type* f) { columns_.push_back(f); return *this; }

    std::ostream& print(std::ostream&) const;
    friend std::ostream& operator<<(std::ostream& s, const Grid& g) {
        return g.print(s);
    }
};

template <typename F>
std::ostream& Grid<F>::print(std::ostream& outs) const 
{
    if (nrows_==0) return outs;
    // double step=(xmax-xmin)/(rows-1);
    double x=xmin_;
    const double step=(nrows_==1)? 0 : (xmax_-xmin_)/(nrows_-1);
    for (unsigned int ir=0; ir<nrows_; ++ir, x+=step) {
        outs << x;
        for (typename fn_vector::const_iterator it=columns_.begin();
             it!=columns_.end(); ++it) {
            outs << " " << (**it)(x);
        }
        outs << "\n";
    }
    return outs;
}


#endif /* BHM_GRID_HPP_INCLUDED_426f62d8250842fb9f11de8116019c61 */
