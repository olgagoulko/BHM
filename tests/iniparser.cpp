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
   @file iniparser.cpp

   Tests the C++ fron-end to the INI file parser
*/

#include <iniparser_frontend.hpp>

#include "sput.h"

#include <iostream>

static void test_ctor()
{
    bool ok=0;
    try {
        iniparser::param p("nosuchfile");
    } catch (const iniparser::MissingFileError& err) {
        ok=1;
        std::cout << "Caught exception: " << err.what() << std::endl;
        sput_fail_unless(err.filename=="nosuchfile", "Filename that caused the error");
    }
    sput_fail_unless(ok, "Exception thrown");
    iniparser::param p("sample.param");
}

static void test_default_ctor()
{
    iniparser::param p;
    sput_fail_unless(p.get("int",1234)==1234, "Default constructor, int val");
    sput_fail_unless(p.get("dbl",12.5)==12.5, "Default constructor, double val");
    sput_fail_unless(p.get("bol",true)==true, "Default constructor, bool val");
    sput_fail_unless(p.get("str1","hi")=="hi", "Default constructor, char* val");
    sput_fail_unless(p.get("str2",std::string("hello"))=="hello", "Default constructor, string val");
}

static void test_load()
{
    iniparser::param p;
    bool ok=0;
    try {
        p.load("nosuchfile");
    } catch (const iniparser::MissingFileError& err) {
        ok=1;
        std::cout << "Caught exception: " << err.what() << std::endl;
        sput_fail_unless(err.filename=="nosuchfile", "Filename that caused the error");
    }
    sput_fail_unless(ok, "Exception thrown");
    sput_fail_unless(p.get(":datapointsmin", 123)==123, "Failed load");

    p.load("sample.param");
    int val=p.get(":datapointsmin", 123);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==100, "Load");
}

static void test_get_int_data()
{
    iniparser::param p("sample.param");
    int val=p.get(":datapointsmin", 123);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==100, "Getting int param");
}

static void test_get_double_data()
{
    iniparser::param p("sample.param");
    double val=p.get(":threshold", 14.25);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==2.1, "Getting double param");
}

static void test_get_bool_data()
{
    iniparser::param p("sample.param");
    bool val=p.get(":verbose", false);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==true, "Getting bool param");
}

static void test_get_string_data()
{
    iniparser::param p("sample.param");
    std::string val=p.get(":data", std::string("NONE"));
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val=="histogram.dat", "Getting string param");
}

static void test_get_string_chardef_data()
{
    iniparser::param p("sample.param");
    std::string val=p.get(":data", "NONE");
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val=="histogram.dat", "Getting string param with char* default");
}

static void test_get_int_nodata()
{
    iniparser::param p("sample.param");
    int val=p.get(":NOSUCH", 123);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==123, "Missing int param");
}

static void test_get_double_nodata()
{
    iniparser::param p("sample.param");
    double val=p.get(":NOSUCH", 14.25);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==14.25, "Missing double param");
}

static void test_get_bool_nodata()
{
    iniparser::param p("sample.param");
    bool val=p.get(":NOSUCH", false);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==false, "Missing bool param");
}

static void test_get_string_nodata()
{
    iniparser::param p("sample.param");
    std::string val=p.get(":NOSUCH", std::string("NONE"));
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val=="NONE", "Missing string param");
}

static void test_get_string_chardef_nodata()
{
    iniparser::param p("sample.param");
    std::string val=p.get(":NOSUCH", "NONE");
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val=="NONE", "Missing string param with char* default");
}

static void test_wrong_type_data()
{
    iniparser::param p("sample.param");
    int val=p.get(":jumpsuppression", 123);
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val==0, "Getting wrongly typed value");
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_iniparser");

    sput_run_test(test_ctor);
    sput_run_test(test_default_ctor);
    sput_run_test(test_load);

    sput_run_test(test_get_int_data);
    sput_run_test(test_get_double_data);
    sput_run_test(test_get_bool_data);
    sput_run_test(test_get_string_data);
    sput_run_test(test_get_string_chardef_data);

    sput_run_test(test_get_int_nodata);
    sput_run_test(test_get_double_nodata);
    sput_run_test(test_get_bool_nodata);
    sput_run_test(test_get_string_nodata);
    sput_run_test(test_get_string_chardef_nodata);
    
    sput_run_test(test_wrong_type_data);

    sput_finish_testing();
    return sput_get_return_value();
}
