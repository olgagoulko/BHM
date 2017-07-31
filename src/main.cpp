#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"

using namespace std;

//definitions of several test functions with corresponding domains and possible basis for basis projection sampling

class testFunction {
	
protected:
	
	double minVar;
	double maxVar;
	double intervalSize;
	double testFunctionMax;
	
	vector<basisSlot*> basisVector;
	int numberBasisFunctions;
	
public:
	
	testFunction() {minVar=1.; maxVar=2.8; intervalSize=2.; testFunctionMax=100; numberBasisFunctions=4; basisVector.resize(0);}
	void setTestFunctionParameters(double theMin, double theMax, double theSize, double theFunMax) {minVar=theMin; maxVar=theMax; intervalSize=theSize; testFunctionMax=theFunMax;}
	
	double getMinVar() const {return minVar;}
	double getMaxVar() const {return maxVar;}
	double getIntervalSize() const {return intervalSize;}
	double getTestFunctionMax() const {return testFunctionMax;}
	vector<basisSlot*> getBasisVector() const {return basisVector;}
	
	virtual double theTestFunctionValue(double variable) const {return 1;}
	
};

class testFunctionCubicPolynomial: public testFunction {
	
private:
	
public:
	
	testFunctionCubicPolynomial() : testFunction()
		{
		setTestFunctionParameters(1., 2.8, 2., 0.6168917686383567);
		slotBounds bounds1(1, 2.8); basisVector.push_back(new taylorSlot(bounds1,numberBasisFunctions)); 
		}
	double theTestFunctionValue(double variable) const {return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;}
	
};

class testFunctionQuatricPolynomial: public testFunction {
	
private:
	
public:
	
	testFunctionQuatricPolynomial() : testFunction()
		{
		setTestFunctionParameters(-1.,1,2.,0.2);
		slotBounds bounds1(-1,-0.5); basisVector.push_back(new taylorSlot(bounds1,numberBasisFunctions)); 
		slotBounds bounds2(-0.5,0); basisVector.push_back(new taylorSlot(bounds2,numberBasisFunctions));
		slotBounds bounds3(0, 0.5); basisVector.push_back(new taylorSlot(bounds3,numberBasisFunctions));
		slotBounds bounds4(0.5, 1); basisVector.push_back(new taylorSlot(bounds4,numberBasisFunctions));
		}
	double theTestFunctionValue(double variable) const {return pow(variable,4)-0.8*variable*variable;}
	
};

class testFunctionExp: public testFunction {
	
private:
	
public:
	
	testFunctionExp() : testFunction()
		{
		setTestFunctionParameters(1.,2.8,2.,3.1);
		slotBounds bounds1(1, 1.9); basisVector.push_back(new taylorSlot(bounds1,numberBasisFunctions)); 
		slotBounds bounds2(1.9, 2.8); basisVector.push_back(new taylorSlot(bounds2,numberBasisFunctions));
		}
	double theTestFunctionValue(double variable) const {return exp(-3*variable)/(-1 + exp(6))*3*exp(9);}
	
};

class testFunctionCos: public testFunction {
	
private:
	
public:
	
	testFunctionCos() : testFunction()
		{
		setTestFunctionParameters(1.,PI+0.6,PI,1.1/PI);
		slotBounds bounds1(1,1.3); basisVector.push_back(new taylorSlot(bounds1,numberBasisFunctions)); 
		slotBounds bounds2(1.3,1.6); basisVector.push_back(new taylorSlot(bounds2,numberBasisFunctions));
		slotBounds bounds3(1.6,1.9); basisVector.push_back(new taylorSlot(bounds3,numberBasisFunctions));
		slotBounds bounds4(1.9,2.2); basisVector.push_back(new taylorSlot(bounds4,numberBasisFunctions));
		slotBounds bounds5(2.2,2.5); basisVector.push_back(new taylorSlot(bounds5,numberBasisFunctions));
		slotBounds bounds6(2.5,2.8); basisVector.push_back(new taylorSlot(bounds6,numberBasisFunctions));
		slotBounds bounds7(2.8,PI); basisVector.push_back(new taylorSlot(bounds7,numberBasisFunctions));
		slotBounds bounds8(PI,PI+0.3); basisVector.push_back(new taylorSlot(bounds8,numberBasisFunctions));
		slotBounds bounds9(PI+0.3,PI+0.6); basisVector.push_back(new taylorSlot(bounds9,numberBasisFunctions));
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


int Main(int argc, char **argv) {
        if (argc!=2) {
            print_help(argv[0], std::cerr);
            return BAD_ARGS;
        }
        iniparser::param par(argv[1]); // FIXME: can throw!

	cout << "----------------- Example BHM code ----------------" << endl;
	
	long unsigned int seed=956475;//time(NULL);
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	//adjust parameters ---------------------------------------------
	long samplingSteps=1e4;
	testFunctionExp myTestFunction;	//select which function to use here
	//---------------------------------------------------------------
	
	double variable, random;
	double minVar=myTestFunction.getMinVar();
	double maxVar=myTestFunction.getMaxVar();
	double intervalSize=myTestFunction.getIntervalSize();
	double testFunctionMax=myTestFunction.getTestFunctionMax();
	double slotWidth=(maxVar-minVar)/pow(2.,10);

	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
	histogramBasis binHistogram(histogramVector);
	histogramBasis basisHistogram(myTestFunction.getBasisVector());

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
		basisHistogram.sample(variable,whatsign(myTestFunction.theTestFunctionValue(variable)));
		}

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

	double threshold=par.get(":Threshold", 2.0);

        bool enableJumpSuppression=par.get(":JumpSuppression", false);
        double jumpSuppression=enableJumpSuppression? 0 : 1.0;

        bool verbose=par.get(":verbose", true);
        bool fail_if_bad=par.get(":FailOnBadFit", true);
        bool fail_if_zero=par.get(":FailOnZeroFit", false);
        bool print_fit=par.get(":PrintFitInfo", false);

        std::string outfile_name=par.get(":OutputName","");
        if (outfile_name.empty()) {
            std::cerr << "Output file name must be provided via 'OutputName' parameter"
                      << std::endl;
            return BAD_ARGS;
        }

        std::ofstream outfile(outfile_name.c_str());
        if (!outfile) {
            std::cerr << "Cannot open output file '" << outfile_name << "'"
                      << std::endl;
            return BAD_ARGS;
        }

        if (ilog2(binHistogram.getSize())<0) {
            std::cerr << "Number of bins (" << binHistogram.getSize()
                      << ") must be a power of 2"
                      << std::endl;
            return BAD_DATA;
        }
        
        if (verbose) {
            std::cout << std::boolalpha
                      << "Input parameters:\n"
                      << "SplineOrder = " << splineOrder-1 << " # spline order\n"
                      << "MinLevel = " << minLevel << " # minimual number of levels per interval\n"
                      << "Threshold = " << threshold << " # minimal goodness-of-fit threshold\n"
                      << "JumpSuppression = " << (jumpSuppression>0) << " # suppression of highest order derivative\n"
                      << "Verbose = " << verbose << "# verbose output\n"
                      << "FailOnZeroFit = " << fail_if_zero << "# do not proceed if the fit is consistent with 0\n"
                      << "FailOnBadFit = " << fail_if_bad << "# do not proceed if the fit is bad\n"
                      << "PrintFitInfo = " << print_fit << "# print the fit information\n"
                      << std::endl;
        }
        
	cout << endl << "BHM fit:" << endl;
	splineArray testBHMfit = binHistogram.BHMfit(splineOrder, minLevel, samplingSteps,
                                                     threshold, jumpSuppression,
                                                     verbose, fail_if_zero);
	cout << endl;
        if (!testBHMfit.getAcceptance()) {
            if (verbose) std::cout << "WARNING: no acceptable fit found" << std::endl;
            if (fail_if_bad) return BAD_FIT;
        }
        if (print_fit) {
            testBHMfit.printSplineArrayInfo(std::cout); cout << endl;
        }
        
	testBHMfit.printSplines(outfile);
        

#if 0
        
	cout << "BHM fit with highest derivative jump suppression:" << endl;
	jumpSuppression=1;
	splineArray testJumpSuppression = binHistogram.BHMfit(splineOrder, minLevel, samplingSteps, threshold, jumpSuppression);
	cout << "Printing spline:" << endl;
	testJumpSuppression.printSplineArrayInfo(); cout << endl;
	testJumpSuppression.printSplines();
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; pair<double,double> basisResult;
	
	//normalization for output
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);
	
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		variable=minVar+i*printStep;
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(variable);
		output << variable << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(variable).first << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(variable).second << '\t' 
		<< '\t' << testBHMfit.splineValue(variable) << '\t' << testBHMfit.splineError(variable)
		<< '\t' << testBHMfit.splineDerivative(variable,1) << '\t' << testBHMfit.splineDerivative(variable,2) << '\t' << testBHMfit.splineDerivative(variable,3)
		<< '\t' << testJumpSuppression.splineValue(variable) << '\t' << testJumpSuppression.splineError(variable)
		<< '\t' << testJumpSuppression.splineDerivative(variable,1) << '\t' << testJumpSuppression.splineDerivative(variable,2) << '\t' << testJumpSuppression.splineDerivative(variable,3)
		<< '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
#endif	
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
