/** @file poly.hpp
    Provides polynomial function from coefficients.
**/

#ifndef BHM_POLY_HPP_INCLUDED_6c2a3e627f1f4861a8c93a6a460873ff
#define BHM_POLY_HPP_INCLUDED_6c2a3e627f1f4861a8c93a6a460873ff

#include <vector>

class Poly {
    typedef std::vector<double> container_type_;
    container_type_ coeff_;
  public:
    /// default constructor, creates 0-valued monomial
    Poly() : coeff_() {}

    /// create a polynomial from sequence of coefficients
    template <typename S>
    Poly(const S& seq) { this->set(seq); }

    /// create a polynomial from pair of sequence iterators
    template <typename IT>
    Poly(IT start, IT end) { this->set(start, end); }

    /// sets coefficients from sequence
    template <typename S>
    Poly& set(const S& seq) { this->set(seq.begin(), seq.end()); return *this; }
    
    /// sets coefficients from pair of sequence iterators
    template <typename IT>
    Poly& set(IT start, IT end) { coeff_.assign(start,end); return *this; }
    
    /// return a polynomial at point
    double operator()(double x) const;
};

#endif /* BHM_POLY_HPP_INCLUDED_6c2a3e627f1f4861a8c93a6a460873ff */
