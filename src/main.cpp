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
#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"
#include "print_spline_grid.hpp"

using namespace std;

static std::ostream& print_help(const char* argv0, std::ostream& strm) {
        strm
            << "Usage:\n"
            << "1) " << argv0 << " param_file.ini\n"
            << "Reads parameters from file `param_file.ini` (that contains names\n"
               "of input sampling histogram and output spline files) and generates the spline\n"
               "approximating the sampled function.\n"
            << "2) " << argv0 << "'' <histogram.dat >spline.dat\n"
            << "Reads histogram from the standard input (here: `histogram.dat`)\n"
               "and writes the spline to the standard output (here: `spline.dat`)\n"
               "keeping all parameters at their default values."
            << std::endl;

        return strm;
}


int Main(int argc, char **argv) {
        if (argc!=2) {
            print_help(argv[0], std::cerr);
            return BAD_ARGS;
        }

        // provision for empty argument as a special case
        iniparser::param par;
        if (argv[1][0]!='\0') {
            par.load(argv[1]); // can throw
        }
        
        std::string outfile_name=par.get(":OutputName","");
	std::ofstream outfile_stream;
        bool default_verbose=false;

        if (!outfile_name.empty()) {
            outfile_stream.open(outfile_name.c_str());
            if (!outfile_stream) {
                std::cerr << "Cannot open file '" << outfile_name << "'" << std::endl;
                return BAD_ARGS;
            }
            default_verbose=true; // verbose by default only if output to a file
        }

        std::ostream& outfile=*(outfile_name.empty()? &std::cout : &outfile_stream);
        bool verbose=par.get(":verbose",default_verbose);

        LOGGER_VERBOSITY(verbose);

        LOGGER << "--------------------------- BHM fit -----------------------------";
        
	unsigned int dataPointsMin=par.get(":DataPointsMin", 100);
	if (dataPointsMin<10) {
		std::cerr << "Warning: DataPointsMin too small, resetting it to 10";
		dataPointsMin=10;
	}
	
        int splinePolynomialOrder=par.get(":SplineOrder", 3);
        if (splinePolynomialOrder<0) {
            std::cerr << "Polynomial order cannot be less than 0" << std::endl;
            return BAD_ARGS;
        }
        unsigned int splineOrder=splinePolynomialOrder+1; // number of polynomial coefficients
        
	unsigned int minLevel=par.get(":MinLevel", 2);
	if(splineOrder >= pow(2,minLevel+1)-1) {
		std::cerr << "Warning: Spline order too high for given MINLEVEL, resetting to defaults";
		minLevel=2; splineOrder=4;
		}
        if (minLevel<2) {
            std::cerr << "Warning: MINLEVEL must be at least 2, resetting it to 2";
            minLevel=2;
        }

        fitAcceptanceThreshold threshold;
	threshold.min=par.get(":THRESHOLD", 2.0);
	threshold.max=par.get(":THRESHOLDMAX", 2.0);	//if max<=min, use only min threshold value, steps set to zero
	threshold.steps=par.get(":THRESHOLDSTEPS", 0);	//if steps==0 only use min threshold value, ignore max
	if(threshold.max<=threshold.min) threshold.steps=0;
	if(threshold.steps<0) threshold.steps=0;
	
	double usableBinFraction=par.get(":UsableBinFraction",0.25);
	if (usableBinFraction>1) {
		std::cerr << "Warning: UsableBinFraction cannot be larger than 1, resetting it to default value 0.25";
		usableBinFraction=0.25;
	}

        bool enableJumpSuppression=par.get(":JumpSuppression", false);
        double jumpSuppression=enableJumpSuppression? 1.0 : 0.0;

        bool fail_if_bad=par.get(":FailOnBadFit", true);
        bool fail_if_zero=par.get(":FailOnZeroFit", false);
        bool print_fit=par.get(":PrintFitInfo", true);


        std::string infile_name=par.get(":Data","");
        // If no infile name is given, use std::cin as the file stream
        std::ifstream infile_stream;
        if (!infile_name.empty()) {
            infile_stream.open(infile_name.c_str());
            if (!infile_stream) {
                std::cerr << "Cannot open input file '" << infile_name << "'"
                          << std::endl;
                return BAD_ARGS;
            }
        }
        std::istream& infile = *(infile_name.empty()? &std::cin : &infile_stream);

        
        const std::string grid_name=par.get(":GridOutput", "");
	unsigned int grid_points=par.get(":GridPoints", 1024);
        std::ofstream grid_outfile;
        if (!grid_name.empty()) {
            grid_outfile.open(grid_name.c_str());
            if (!grid_outfile) {
                std::cerr << "Cannot open file '" << grid_name << "'" << std::endl;
                return BAD_ARGS;
            }
            if(grid_points<1) {
		    std::cerr << "Too few grid points '" << grid_points << "'" << std::endl;
		    return BAD_ARGS;
	    }
        }

        
        histogramBasis binHistogram(infile);
        
        if (!infile_name.empty()) infile_stream.close();
     
	unsigned int histogramPower=(unsigned int)ilog2(binHistogram.getSize());
        if (ilog2(binHistogram.getSize())<0) {
            std::cerr << "Number of bins (" << binHistogram.getSize()
                      << ") must be a power of 2"
                      << std::endl;
            return BAD_DATA;
        }
        else if (binHistogram.getSize()<4) { //this is still very few bins
		std::cerr << "Number of bins (" << binHistogram.getSize() << ") is too small" << std::endl;
		return BAD_DATA;
	}
	else if (minLevel+2>histogramPower) {
		std::cerr << "Warning: MINLEVEL too large for histogram size, resetting it to 2";
		minLevel=2;
	}
	
	if(binHistogram.getNumberOfSamples()<dataPointsMin) { //this is still very few data points
		std::cerr << "Not enough sampled data points (" << binHistogram.getNumberOfSamples() << ") in histogram" << std::endl;
		return BAD_DATA;
	}
        
        LOGGER << std::boolalpha
               << "Input parameters:\n"
	       << left << setw(20) << "DataPointsMin = " << left << setw(20) << dataPointsMin 	<< " # minimal number of data points per bin\n"
               << setw(20) << "SplineOrder = " 		<< left << setw(20) << splineOrder-1 	<< " # spline order\n"
               << setw(20) << "MinLevel = " 		<< left << setw(20) << minLevel 	<< " # minimual number of levels per interval\n"
               << setw(20) << "Threshold = " 		<< left << setw(20) << threshold.min 	<< " # minimal goodness-of-fit threshold\n"
               << setw(20) << "ThresholdMax = "		<< left << setw(20) << threshold.max 	<< " # maximal goodness-of-fit threshold (if applicable)\n"
               << setw(20) << "ThresholdSteps = " 	<< left << setw(20) << threshold.steps << " # number of steps for goodness-of-fit threshold increase\n"
	       << setw(20) << "UsableBinFraction = " 	<< left << setw(20) << usableBinFraction << " # minimal proportion of good bins for a level to be considered\n"
               << setw(20) << "JumpSuppression = " 	<< left << setw(20) << (jumpSuppression>0) << " # suppression of highest order derivative\n"
               << setw(20) << "Verbose = " 		<< left << setw(20) << verbose 		<< " # verbose output\n"
               << setw(20) << "FailOnZeroFit = " 	<< left << setw(20) << fail_if_zero 	<< " # do not proceed if the fit is consistent with 0\n"
               << setw(20) << "FailOnBadFit = " 	<< left << setw(20) << fail_if_bad 	<< " # do not proceed if the fit is bad\n"
               << setw(20) << "PrintFitInfo = " 	<< left << setw(20) << print_fit 	<< " # print the fit information\n"
               << setw(20) << "Data = " 		<< left << setw(20) << infile_name 	<< " # Input histogram file\n"
               << setw(20) << "OutputName = " 		<< left << setw(20) << outfile_name 	<< " # Output file to print results to\n"
               << setw(20) << "GridOutput = " 		<< left << setw(20) << grid_name 	<< " # Output file to print spline values on a grid\n"
	       << setw(20) << "GridPoints = " 		<< left << setw(20) << grid_points 	<< " # number of grid points for spline output\n"
               << "\nInput histogram:"
               << setw(20) << "\nNumber of bins: " << binHistogram.getSize()
               << setw(20) << "\nNumber of samples: " << binHistogram.getNumberOfSamples()
	       << setw(20) << "\nNormalization: " << binHistogram.getNorm()
               << setw(20) << "\nLower bound: " << binHistogram.getLowerBound()
               << setw(20) << "\nUpper bound: " << binHistogram.getUpperBound() << "\n";

        LOGGER << "BHM fit:";
	
	BHMparameters theParameters;
	theParameters.dataPointsMin=dataPointsMin;
	theParameters.splineOrder=splineOrder;
	theParameters.minLevel=minLevel;
	theParameters.threshold=threshold;
	theParameters.usableBinFraction=usableBinFraction;
	theParameters.jumpSuppression=jumpSuppression;
        
	histogramBasis normalizedBinHistogram = binHistogram.normalizedHistogram(binHistogram.getNorm());
	splineArray testBHMfit = normalizedBinHistogram.BHMfit(theParameters, binHistogram.getNumberOfSamples(), fail_if_zero);

        if (!testBHMfit.getAcceptance()) {
            LOGGER << "WARNING: no acceptable fit found";
            if (fail_if_bad) return BAD_FIT;
        }
        if (print_fit) {
            testBHMfit.printSplineArrayInfo(std::cout);
        }
        
	testBHMfit.printSplines(outfile);

        // Should we dump the grid?
        if (grid_outfile) print_spline_grid(grid_outfile, testBHMfit, grid_points);
        
	return OK;
	
}

int main(int argc, char** argv)
{
    try {
        return Main(argc, argv);
    } catch (const iniparser::Error& err) {
        std::cerr << "Cannot read parameters: " << err.what() << std::endl;
        return BAD_ARGS;
    } catch (const histogramBasis::ConsistentWithZero_Error& err) {
        std::cerr << err.what() << std::endl;
        return ZERO_DATA;
    } catch (const histogramBasis::NotEnoughData_Error& err) {
        std::cerr << err.what() << std::endl;
        return BAD_DATA;
    } catch (const std::runtime_error& err) {
        std::cerr << "Runtime error: " << err.what() << std::endl;
        return OTHER_ERROR;
    } catch (const std::logic_error& err) {
        std::cerr << "Internal logic error: " << err.what() << std::endl;
        return OTHER_ERROR;
    }
    // can never reach here
}
