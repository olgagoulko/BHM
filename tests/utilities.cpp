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
