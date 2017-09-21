/**
   @file poly_test.cpp

   Tests polynomials
*/

#include <vector>

#include <cmath>
#include <float.h>

#include "sput.h"

#include <poly.hpp>

/// Compares two floating point values, see http://floating-point-gui.de/errors/comparison/
static bool is_near(const double a, const double b, const double eps=DBL_EPSILON)
{
    using std::fabs;
    if (a==b) return true;
    
    const double diff=fabs(a-b);
    if (a==0 || b==0 || diff<DBL_MIN) return diff < (eps*DBL_MIN);

    const double abs_a=fabs(a);
    const double abs_b=fabs(b);
    const double sum = (abs_b > DBL_MAX-abs_a)? DBL_MAX : (abs_a+abs_b);
    return diff/sum < eps;
}

static void ctor()
{
    Poly p;
    sput_fail_unless(p(0)==0, "Default ctor value 1");
    sput_fail_unless(p(43.25)==0, "Default ctor value 2");
}

// Global variables for testing
namespace Test {
    static const double coeff[]={1.1, 2.2, 3.3};
    static const std::size_t coeff_size=sizeof(coeff)/sizeof(*coeff);
    
    static const double xset[]={-2, -1, 0, 1, 2};
    static const std::size_t xset_size=sizeof(xset)/sizeof(*xset);
}

static void poly_value_check(const Poly& p)
{
    using Test::coeff;
    using Test::xset;
    using Test::xset_size;
    
    for (const double* px=xset; px != xset+xset_size; ++px) {
        const double y = coeff[0] + coeff[1]*(*px) + coeff[2]*(*px)*(*px);
        sput_fail_unless(is_near(y, p(*px)), "value check");
    }
}

static void from_sequence()
{
    std::vector<double> p_coeff(Test::coeff, Test::coeff+Test::coeff_size);
    Poly p;
    p.set(p_coeff);
    poly_value_check(p);
}

static void ctor_from_sequence()
{
    std::vector<double> p_coeff(Test::coeff, Test::coeff+Test::coeff_size);
    Poly p(p_coeff);
    poly_value_check(p);
}

static void from_iters()
{
    Poly p;
    p.set(Test::coeff, Test::coeff+Test::coeff_size);
    poly_value_check(p);
}

static void ctor_from_iters()
{
    Poly p(Test::coeff, Test::coeff+Test::coeff_size);
    poly_value_check(p);
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_poly_test");

    sput_run_test(ctor);
    sput_run_test(from_sequence);
    sput_run_test(ctor_from_sequence);
    sput_run_test(from_iters);
    sput_run_test(ctor_from_iters);

    sput_finish_testing();
    return sput_get_return_value();
}
