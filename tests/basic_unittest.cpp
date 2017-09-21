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
	sput_fail_unless(isAround(1.+0.9999*ACCURACY,1.)  == true, "1+EPS == 1");
	
	sput_fail_unless(isAround(1.,2.)  == false, "1 != 2");
	sput_fail_unless(isAround(1.,-1.)  == false, "1 != -1");
	sput_fail_unless(isAround(1.,1.+2*ACCURACY)  == false, "1 != 1+2EPS");
	sput_fail_unless(isAround(1.+2*ACCURACY,1.)  == false, "1+2EPS != 1");
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

