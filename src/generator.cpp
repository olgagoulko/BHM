#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"

using namespace std;

//definitions of several test functions with corresponding domains for sampling with rejection method

class testFunction {
	
private:
	
	int numberOfFunction;
	double minVar;
	double maxVar;
	double intervalSize;
	double testFunctionMax;
	
public:
	
	testFunction(int which)
		{
		numberOfFunction=which;
		if(which==0) {minVar=1.; maxVar=2.8; intervalSize=2.; testFunctionMax=0.6168917686383567;}
		else if(which==1) {minVar=-1.; maxVar=1.0; intervalSize=2.; testFunctionMax=0.2;}
		else if(which==2) {minVar=1.; maxVar=2.8; intervalSize=2.; testFunctionMax=3.1;}
		else if(which==3) {minVar=1.; maxVar=PI+0.6; intervalSize=PI; testFunctionMax=1.1/PI;}
		else if(which==4) {minVar=-5; maxVar=5.; intervalSize=10.; testFunctionMax=0.5;}	//don't use rejection method for this one
		}
	
	double getMinVar() const {return minVar;}
	double getMaxVar() const {return maxVar;}
	double getIntervalSize() const {return intervalSize;}
	double getTestFunctionMax() const {return testFunctionMax;}
	
	double theTestFunctionValue(double variable) const
		{
		double value=1;
		if(numberOfFunction==0) {value=3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;}
		else if(numberOfFunction==1) {value=pow(variable,4)-0.8*variable*variable;}
		else if(numberOfFunction==2) {value=exp(-3*variable)/(-1 + exp(6))*3*exp(9);}
		else if(numberOfFunction==3) {value=(10+cos(variable*10.))/10./PI;}
		return value;
		}
	
};


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
        iniparser::param par;
        if (argv[1][0]!='\0') {
            par.load(argv[1]); // can throw
        }

        std::string out_hist_name=par.get(":OUTPUT","");
	std::ofstream output_file_stream;
        bool default_verbose=false;
        if (!out_hist_name.empty()) {
            output_file_stream.open(out_hist_name.c_str());
            if (!output_file_stream) {
                std::cerr << "Cannot open file '" << out_hist_name << "'" << std::endl;
                return BAD_ARGS;
            }
            default_verbose=true; // verbose by default only if output to a file
        }
        std::ostream& output=*(out_hist_name.empty()? &std::cout : &output_file_stream);
        bool verbose=par.get(":VERBOSE",default_verbose);
        
	if (verbose) {
            cout << "Generating example histogram" << endl;
        }
	
	long unsigned int seed=par.get(":RANDOMSEED", 0);
	if(seed==0) seed=time(NULL);
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=par.get(":SAMPLESIZE", 1e4);
	if(samplingSteps<defaultMinNumberTimesSampled) {cerr << "Too few sampling steps" << endl; return BAD_ARGS;}

	int theFunctionChoice=par.get(":FUNCTION",0);
	if(theFunctionChoice<0 || theFunctionChoice>4) {cerr << "WARNING: No such function, will instead use function 0 (cubic polynomial)" << endl; theFunctionChoice=0;}
	testFunction myTestFunction(theFunctionChoice);
	
	double variable, random;
	double minVar=myTestFunction.getMinVar();
	double maxVar=myTestFunction.getMaxVar();
	double intervalSize=myTestFunction.getIntervalSize();
	double testFunctionMax=myTestFunction.getTestFunctionMax();
	
	unsigned int power=par.get(":POWERBINS",10);
	unsigned int numberHistogramBins=(1<<power);
	double slotWidth=(maxVar-minVar)/double(numberHistogramBins);

	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
	histogramBasis binHistogram(histogramVector);

	//rejection method to generate x with appropriate probabilities
	for(int i=0; i<samplingSteps;i++)
		{
		if(theFunctionChoice!=4)
			{
			bool accept=false;
			while(accept==false)
				{
				variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
				random=gsl_rng_uniform(RNG)*testFunctionMax;
				if(random<abs(myTestFunction.theTestFunctionValue(variable))) accept=true;
				}
			}
		else
			{
			if(i%5==0) {variable = gsl_ran_gaussian (RNG, 0.2);}
			else if(i%2==0) {variable = 2+ gsl_ran_gaussian (RNG, 1.0);}
			else {variable = -2+ gsl_ran_gaussian (RNG, 1.0);}
			}
		
		binHistogram.sampleUniform(variable,whatsign(myTestFunction.theTestFunctionValue(variable)));
		}


	for(unsigned int i=0; i<numberHistogramBins;i++) 
		{
		basisSlot* currentSlot = binHistogram.getSlot(i);
		output << setprecision(12) << currentSlot->getBounds().getLowerBound() << '\t' << currentSlot -> getNumberTimesSampled() << '\t' << currentSlot->sampledIntegral() << '\t' << currentSlot->getVariance() << endl;
		}
	output << maxVar << endl;

	return OK;
	
}

int main(int argc, char** argv)
{
    try {
        return Main(argc, argv);
    } catch (const iniparser::Error& e) {
        std::cerr << "Error parsing input file: " << e.what() << std::endl;
        return BAD_ARGS;
    } catch (const std::runtime_error& e) {
        std::cerr << "Runtime error: " << e.what() << std::endl;
        return OTHER_ERROR;
    }
}
