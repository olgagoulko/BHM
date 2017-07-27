#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"

using namespace std;

//definitions of several test functions with corresponding domains for sampling with rejection method

class testFunction {
	
protected:
	
	double minVar;
	double maxVar;
	double intervalSize;
	double testFunctionMax;
	
public:
	
	testFunction() {minVar=1.; maxVar=2.8; intervalSize=2.; testFunctionMax=100;}
	void setTestFunctionParameters(double theMin, double theMax, double theSize, double theFunMax) {minVar=theMin; maxVar=theMax; intervalSize=theSize; testFunctionMax=theFunMax;}
	
	double getMinVar() const {return minVar;}
	double getMaxVar() const {return maxVar;}
	double getIntervalSize() const {return intervalSize;}
	double getTestFunctionMax() const {return testFunctionMax;}
	
	virtual double theTestFunctionValue(double variable) const {return 1;}
	
};

class testFunctionCubicPolynomial: public testFunction {
	
private:
	
public:
	
	testFunctionCubicPolynomial() : testFunction()
		{
		setTestFunctionParameters(1., 2.8, 2., 0.6168917686383567);
		}
	double theTestFunctionValue(double variable) const {return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;}
	
};

class testFunctionQuatricPolynomial: public testFunction {
	
private:
	
public:
	
	testFunctionQuatricPolynomial() : testFunction()
		{
		setTestFunctionParameters(-1.,1,2.,0.2);
		}
	double theTestFunctionValue(double variable) const {return pow(variable,4)-0.8*variable*variable;}
	
};

class testFunctionExp: public testFunction {
	
private:
	
public:
	
	testFunctionExp() : testFunction()
		{
		setTestFunctionParameters(1.,2.8,2.,3.1);
		}
	double theTestFunctionValue(double variable) const {return exp(-3*variable)/(-1 + exp(6))*3*exp(9);}
	
};

class testFunctionCos: public testFunction {
	
private:
	
public:
	
	testFunctionCos() : testFunction()
		{
		setTestFunctionParameters(1.,PI+0.6,PI,1.1/PI);
		}
	double theTestFunctionValue(double variable) const {return (10+cos(variable*10.))/10./PI;}
	
};


static std::ostream& print_help(const char* argv0, std::ostream& strm) {
        strm
            << "FIXME!!! Meaningful help should be printed here\n"
            << "Usage:\n" << argv0 << " param_file.ini"
            << std::endl;

        return strm;
}
            


int main(int argc, char **argv) {

        if (argc!=2) {
            print_help(argv[0], std::cerr);
            return 1;
        }
        iniparser::param par(argv[1]); // FIXME: can throw!

	cout << "----------------- Example BHM code ----------------" << endl;
	
	long unsigned int seed=par.get(":RANDOMSEED", 0);
	if(seed==0) seed=time(NULL);
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=par.get(":SAMPLESIZE", 1e4);
	if(samplingSteps<defaultMinNumberTimesSampled) {cerr << "Too few sampling steps" << endl; return 1;}
	
	
	testFunctionExp myTestFunction;
	
	double variable, random;
	double minVar=myTestFunction.getMinVar();
	double maxVar=myTestFunction.getMaxVar();
	double intervalSize=myTestFunction.getIntervalSize();
	double testFunctionMax=myTestFunction.getTestFunctionMax();
	
	double slotWidth=(maxVar-minVar)/pow(2.,10);

	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
	histogramBasis binHistogram(histogramVector);

	//rejection method to generate x with appropriate probabilities
	for(int i=0; i<samplingSteps;i++)
		{
		bool accept=false;
		while(accept==false)
			{
			variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
			random=gsl_rng_uniform(RNG)*testFunctionMax;
			if(random<abs(myTestFunction.theTestFunctionValue(variable))) accept=true;
			}
		
		binHistogram.sampleUniform(variable,whatsign(myTestFunction.theTestFunctionValue(variable)));
		}

        
	unsigned int minLevel=par.get(":MinLevel", 2);
        if (minLevel<2) {
            std::cerr << "Warning: MINLEVEL must be at least 2, resetting it to 2";
            minLevel=2;
        }

	double threshold=par.get(":Threshold", 2.0);


            std::cout << std::boolalpha
                      << "Input parameters:\n"
                      << "SplineOrder = " << splineOrder-1 << " # spline order\n"
                      << "MinLevel = " << minLevel << " # minimual number of levels per interval\n"
                      << "Threshold = " << threshold << " # minimal goodness-of-fit threshold\n"

                      << std::endl;

        

	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; pair<double,double> basisResult;
	
	//normalization for output
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		variable=minVar+i*printStep;
		output << variable << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(variable).first << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(variable).second << '\t' 

		<< endl;
		}
	

	return 0;
	
}
