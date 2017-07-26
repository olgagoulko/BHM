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
    } catch (const std::runtime_error& err) {
        ok=1;
        std::cout << "Caught exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Exception not thrown");
    iniparser::param p("sample.param");
    // int val=p.get<int>("DATAPOINTSMIN");
    // sput_fail_unless(val==100, "Parsing INI file");
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
    sput_fail_unless(val=="input.dat", "Getting string param");
}

static void test_get_string_chardef_data()
{
    iniparser::param p("sample.param");
    std::string val=p.get(":data", "NONE");
    std::cout << "value=" << val << std::endl;
    sput_fail_unless(val=="input.dat", "Getting string param with char* default");
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
