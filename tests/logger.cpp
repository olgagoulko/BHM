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

    LOGGER << "Default stream output: THIS LINE SHOULD NOT BE VISIBLE";
    LOGGER_OUTPUT(std::cout);
    LOGGER << "Default verbosity: THIS LINE SHOULD NOT BE VISIBLE";

    sput_run_test(logger_macro);
    sput_run_test(logger_quiet);
    sput_run_test(logger_verbose);

    sput_finish_testing();
    return sput_get_return_value();
}
