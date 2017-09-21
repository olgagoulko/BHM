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
   @file read_histogram.cpp

   Tests reading histogram from a stream
*/

#include "sput.h"
#include <sstream>

#include <histogram.hpp>

static void read_empty()
{
    const std::string input="";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Empty file exception");
}
    
static void read_incomplete()
{
    const std::string input=
        "5.75\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Incomplete file exception");
}
    
static void read_incomplete2()
{
    const std::string input=
        "5.75 0\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Incomplete file exception");
}
    
static void read_eof()
{
    const std::string input=
        "5.5   1\n"
        "1.25  10  1.  0.\n"
        "2.50  20  1.  0.\n"
        "3.50  30  1.  0.\n"
        "4.25  40  1.  0.\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Unexpected EOF exception");
}
    
static void read_garbage1()
{
    const std::string input=
        "something not a number\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Garbage in file exception");
}

static void read_garbage2()
{
    const std::string input=
        "5.0   1\n"
        "1.25  10  1.  0.\n"
        "2.50  20  1.  0.\n"
        "3.50  30  1.  0.\n"
        "4.25  40  1.  0.\n"
        "something not a number\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Garbage in file exception");
}

static void read_malformatted()
{
    const std::string input=
        "5.0   1\n"
        "1.25  10  1.  0.\n"
        "2.50  20  1.  0.\n"
        "3.50  30  1.  0.\n"
        "4.25  40  1.\n"     //< malformed
        "5.25\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::InvalidFileFormatError& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Invalid format exception");
}

static void read_overlapped()
{
    const std::string input=
        "5.0   1\n"
        "1.25  10  1.  0.\n"
        "2.50  20  1.  0.\n"
        "2.00  30  1.  0.\n" //< overlapping
        "4.25  40  1.  0.\n"
        "5.25\n";
    std::istringstream istrm(input);

    bool ok=false;
    try {
        histogramBasis hist(istrm);
    } catch (histogramBasis::OverlappingSlot& err) {
        ok=true;
        std::cout << "Exception: " << err.what() << std::endl;
    }
    sput_fail_unless(ok, "Overlapping slot exception");
}

static void read_bounds()
{
    const std::string input=
        " 5.0   1\n"
        "-8.75  10  1.  0.\n"
        "-7.5   20  0.5 0.\n"
        "-6.5   30  1.  0.25\n"
        "-5.75  40  1.5  0.75\n"
        "-4.25  45\n"
        "-3.25  5\n"
        "-3.0\n";

    const double slots[6][5]=
    {
        {-8.75, -7.5,  10,  1.0,  0.00},
        {-7.5,  -6.5,  20,  0.5,  0.00},
        {-6.5,  -5.75, 30,  1.0,  0.25},
        {-5.75, -4.25, 40,  1.5,  0.75},
        {-4.25, -3.25, 45,  1.0,  0.00},
        {-3.25, -3.0,   5,  1.0,  0.00}
    };
    const std::size_t nslots=sizeof(slots)/sizeof(*slots);
    
    std::istringstream istrm(input);
    histogramBasis hist(istrm);

    sput_fail_unless(hist.getLowerBound()==slots[0][0], "Histogram lower bound");
    sput_fail_unless(hist.getUpperBound()==slots[nslots-1][1], "Histogram upper bound");
    sput_fail_unless(hist.hasUpperBound(), "Histogram has upper bound");
    sput_fail_unless(hist.getSize()==nslots, "Histogram size");
}    

static void read()
{
    const std::string input=
        "5.0   2\n"
        "1.25  10  1.  0.\n"
        "2.50  20  0.5 0.\n"
        "3.50  30  1.  0.25\n"
        "4.25  40  1.5  0.75\n"
        "5.75  45\n"
        "6.75  5\n"
        "7.00\n";

    const double excess_val=5.0;
    const long excess_count=2;
    const double slots[6][5]=
    {
        {1.25, 2.50, 10,  1.0,  0.00},
        {2.50, 3.50, 20,  0.5,  0.00},
        {3.50, 4.25, 30,  1.0,  0.25},
        {4.25, 5.75, 40,  1.5,  0.75},
        {5.75, 6.75, 45,  1.0,  0.00},
        {6.75, 7.00,  5,  1.0,  0.00}
    };
    const std::size_t nslots=sizeof(slots)/sizeof(*slots);
    
    std::istringstream istrm(input);
    histogramBasis hist(istrm);

    sput_fail_unless(hist.getLowerBound()==slots[0][0], "Histogram lower bound");
    sput_fail_unless(hist.getUpperBound()==slots[nslots-1][1], "Histogram upper bound");
    sput_fail_unless(hist.hasUpperBound(), "Histogram has upper bound");
    sput_fail_unless(hist.getSize()==nslots, "Histogram size");
    if (nslots!=hist.getSize()) return;

    unsigned long nsamples=0;
    for (unsigned int i=0; i<hist.getSize(); ++i) {
        nsamples += slots[i][2];
        sput_fail_unless(hist.getSlot(i)->getBounds().getLowerBound() == slots[i][0],
                         "Slot lower bound");
        sput_fail_unless(hist.getSlot(i)->getBounds().getUpperBound() == slots[i][1],
                         "Slot upper bound");
        sput_fail_unless(hist.getSlot(i)->getNumberTimesSampled() == slots[i][2],
                         "Slot number of samples");
        sput_fail_unless(hist.getSlot(i)->sampledIntegral() == slots[i][3],
                         "Integral");
        sput_fail_unless(hist.getSlot(i)->getVariance() == slots[i][4],
                         "m2");
    }
    sput_fail_unless(hist.getNumberOfSamples()==nsamples+excess_count, "Histogram number of samples");
    sput_fail_unless(hist.getExcessCounter()==excess_count, "Histogram excess count");
    sput_fail_unless(hist.getExcessValues(1.0)==excess_val, "Histogram excess value");
}

int main(int argc, char **argv)
{
    sput_start_testing();
	
    sput_enter_suite("test_read_histogram");

    sput_run_test(read_empty);
    sput_run_test(read_incomplete);
    sput_run_test(read_incomplete2);
    sput_run_test(read_eof);
    sput_run_test(read_garbage1);
    sput_run_test(read_garbage2);
    sput_run_test(read_malformatted);
    sput_run_test(read_overlapped);
    sput_run_test(read_bounds);
    sput_run_test(read);

    sput_finish_testing();
    return sput_get_return_value();
}
