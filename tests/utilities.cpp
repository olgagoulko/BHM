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
/**
   @file utilities.cpp

   Tests utility functions
*/

#include <basic.hpp>

#include "sput.h"

#include <sstream>

// for diagnostics
static std::string toString(int i)
{
    std::ostringstream strbuf;
    strbuf << i;
    return strbuf.str();
}

static void test_ilog2()
{
    for (int i=0; i<sizeof(int)*8-1; ++i) {
        int n=1<<i;
        std::string msg=toString(n)+"=2^"+toString(i);
        sput_fail_unless(ilog2(n)==i, msg.c_str());
        if (i>1) {
            msg=toString(n+1)+"=2^"+toString(i)+"+1";
            sput_fail_unless(ilog2(n+1)==-1,msg.c_str());

            msg=toString(n-1)+"=2^"+toString(i)+"-1";
            sput_fail_unless(ilog2(n-1)==-1,msg.c_str());
        }
    }
    sput_fail_unless(ilog2(0)==-1,"zero is 2^n ?");
    sput_fail_unless(ilog2(-1)==-1,"negative 1 is 2^n ?");
    sput_fail_unless(ilog2(-2)==-1,"negative 2 is 2^n ?");
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_power2_check");

    sput_run_test(test_ilog2);

    sput_finish_testing();
    return sput_get_return_value();
}
