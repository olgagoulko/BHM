/**
   @file logger.cpp

   Tests the logger
*/

#include "sput.h"
#include <logger.hpp>
#include <sstream>

static void logger_macro()
{
    LOGGER_VERBOSITY(1);
    LOGGER << "Logger macro passed";
    sput_fail_unless(true, "Logger macro");
}

static void logger_quiet()
{
    LOGGER_VERBOSITY(0);
    std::ostringstream ostr;
    logger::LogLine(ostr) << "Something=" << 123;
    sput_fail_unless(ostr.str().empty(), "Quiet logging to a stream");
}

static void logger_verbose()
{
    LOGGER_VERBOSITY(1);
    std::ostringstream ostr;
    logger::LogLine(ostr) << "Something=" << 123;
    sput_fail_unless(ostr.str()=="Something=123\n", "Verbose logging to a stream");
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_logger");

    LOGGER << "Default verbosity: THIS LINE SHOULD NOT BE VISIBLE";

    sput_run_test(logger_macro);
    sput_run_test(logger_quiet);
    sput_run_test(logger_verbose);

    sput_finish_testing();
    return sput_get_return_value();
}
