#include "basic.hpp"
#include "sput.h"

using namespace std; 

static void test_isAround()
{
	sput_fail_unless(isAround(1.,1.)  == true, "1 == 1");
	sput_fail_unless(isAround(5.,5.)  == true, "5 == 5");
	sput_fail_unless(isAround(1.6e121,1.6e121)  == true, "1.6e121 == 1.6e121");
	sput_fail_unless(isAround(0.,0.)  == true, "0 == 0");
	sput_fail_unless(isAround(-1.,-1.)  == true, "-1 == -1");
	sput_fail_unless(isAround(1.+0.9999*VERY_SMALL_NUMBER,1.)  == true, "1+EPS == 1");
	
	sput_fail_unless(isAround(1.,2.)  == false, "1 != 2");
	sput_fail_unless(isAround(1.,-1.)  == false, "1 != -1");
	sput_fail_unless(isAround(1.,1.+2*VERY_SMALL_NUMBER)  == false, "1 != 1+2EPS");
	sput_fail_unless(isAround(1.+2*VERY_SMALL_NUMBER,1.)  == false, "1+2EPS != 1");
	sput_fail_unless(isAround(1.,0.)  == false, "1 != 0");
	sput_fail_unless(isAround(0.,1.)  == false, "0 != 1");
}

static void test_whatsign()
{
	sput_fail_unless(whatsign(1.)  == 1, "1. is positive");
	sput_fail_unless(whatsign(1.6e121)  == 1, "1.6e121 is positive");
	sput_fail_unless(whatsign(123)  == 1, "123 is positive");
	sput_fail_unless(whatsign(0)  == 1, "0 is positive");
	sput_fail_unless(whatsign(-1.)  == -1, "-1. is negative");
	sput_fail_unless(whatsign(-1.6e121)  == -1, "-1.6e121 is negative");
	sput_fail_unless(whatsign(-24)  == -1, "-24 is negative");
}

static void test_rounding()
{
	sput_fail_unless(rounding(1.0)  == 1, "1.0 rounded is 1");
	sput_fail_unless(rounding(1.9)  == 2, "1.9 rounded is 2");
	sput_fail_unless(rounding(1.5)  == 2, "1.5 rounded is 2");
	sput_fail_unless(rounding(0.2)  == 0, "0.2 rounded is 0");
	sput_fail_unless(rounding(-0.2)  == 0, "-0.2 rounded is 0");
	sput_fail_unless(rounding(-0.7)  == -1, "-0.7 rounded is -1");
	sput_fail_unless(rounding(-2.5)  == -3, "-2.5 rounded is -3");
	sput_fail_unless(rounding(-4.)  == -4, "-4 rounded is -4");
}

int main(int argc, char **argv) {
	
	sput_start_testing();
	
	sput_enter_suite("test_isAround()");
	sput_run_test(test_isAround);
	sput_enter_suite("test_whatsign()");
	sput_run_test(test_whatsign);
	sput_enter_suite("test_rounding()");
	sput_run_test(test_rounding);
	
	sput_finish_testing();
	
	return sput_get_return_value();
	
}

