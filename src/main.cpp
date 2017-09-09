#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"

using namespace std;

static std::ostream& print_help(const char* argv0, std::ostream& strm) {
        strm
            << "FIXME!!! Meaningful help should be printed here\n"
            << "Usage:\n" << argv0 << " param_file.ini"
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

        LOGGER << "----------------- Example BHM code ----------------";
        
        int splinePolynomialOrder=par.get(":SplineOrder", 4);
        if (splinePolynomialOrder<0) {
            std::cerr << "Polynomial order cannot be less than 0" << std::endl;
            return BAD_ARGS;
        }
        unsigned int splineOrder=splinePolynomialOrder+1; // number of polynomial coefficients
        
	unsigned int minLevel=par.get(":MinLevel", 2);
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

        histogramBasis binHistogram(infile);
        
        if (!infile_name.empty()) infile_stream.close();
     
        if (ilog2(binHistogram.getSize())<0) {
            std::cerr << "Number of bins (" << binHistogram.getSize()
                      << ") must be a power of 2"
                      << std::endl;
            return BAD_DATA;
        }
        
        LOGGER << std::boolalpha
               << "Input parameters:\n"
               << "SplineOrder = " << splineOrder-1 << " # spline order\n"
               << "MinLevel = " << minLevel << " # minimual number of levels per interval\n"
               << "Threshold = " << threshold.min << " # minimal goodness-of-fit threshold\n"
               << "ThresholdMax = " << threshold.max << " # maximal goodness-of-fit threshold (if applicable)\n"
               << "ThresholdSteps = " << threshold.steps << " # steps for goodness-of-fit threshold increase\n"
               << "JumpSuppression = " << (jumpSuppression>0) << " # suppression of highest order derivative\n"
               << "Verbose = " << verbose << " # verbose output\n"
               << "FailOnZeroFit = " << fail_if_zero << "# do not proceed if the fit is consistent with 0\n"
               << "FailOnBadFit = " << fail_if_bad << " # do not proceed if the fit is bad\n"
               << "PrintFitInfo = " << print_fit << " # print the fit information\n"
               << "Data = \"" << infile_name << "\" # Input histogram\n"
               << "OutputName = \"" << outfile_name << "\" # Output file to print results to\n"
               << "\nInput histogram:"
               << "\nNumber of bins: " << binHistogram.getSize()
               << "\nNumber of samples: " << binHistogram.getNumberOfSamples()
               << "\nLower bound: " << binHistogram.getLowerBound()
               << "\nUpper bound: " << binHistogram.getUpperBound();

        LOGGER << "BHM fit:";
        
	splineArray testBHMfit = binHistogram.BHMfit(splineOrder, minLevel, binHistogram.getNumberOfSamples(),
                                                     threshold, jumpSuppression, fail_if_zero);

        if (!testBHMfit.getAcceptance()) {
            LOGGER << "WARNING: no acceptable fit found";
            if (fail_if_bad) return BAD_FIT;
        }
        if (print_fit) {
            testBHMfit.printSplineArrayInfo(std::cout);
        }
        
	testBHMfit.printSplines(outfile);

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
