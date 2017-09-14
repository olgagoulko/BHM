#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include "iniparser_frontend.hpp"

using namespace std;

//definitions of several test functions with corresponding domains for sampling with rejection method

/// Generic template for a type that holds value of type T, but has some specific name other than T
/** Objects of this type must be explicitly initialized by a value of type T,
    but they are easily convertible to type T
*/ 
template <int TAG, typename T=double> // TAG to make otherwise identical types differ
struct holder {
    T value;
    explicit holder(const T& v): value(v) {}
    operator T() const { return value; }
    T operator()() const { return value; }
};

typedef holder<0> min_var;
typedef holder<1> max_var;
typedef holder<2> interval_size;
typedef holder<3> test_function_max;
typedef holder<4,std::string> py_expr;

// This allows us to wrap each double (or string) argument into its own type,
// so that there will be no danger of accidently supplying `min_var` value where `max_var` was expected.

class testFunction_base {
private:
	double minVar;
	double maxVar;
	double intervalSize;
	double testFunctionMax;
        std::string pyExpr;

public:
        testFunction_base(min_var a_minVar, max_var a_maxVar, interval_size a_intervalSize,
                          test_function_max a_testFunctionMax, py_expr a_expr):
            minVar(a_minVar), maxVar(a_maxVar), intervalSize(a_intervalSize),
            testFunctionMax(a_testFunctionMax), pyExpr(a_expr)
        {}

        virtual ~testFunction_base() {}

	double getMinVar() const {return minVar;}
	double getMaxVar() const {return maxVar;}
	double getIntervalSize() const {return intervalSize;}
	double getTestFunctionMax() const {return testFunctionMax;}
        virtual std::string getExpr() const {
            return "def fn(x): return " + pyExpr;
        }
        
        virtual double theTestFunctionValue(double variable) const =0;
};

struct CubicPolynomial : public testFunction_base
{
    CubicPolynomial(): testFunction_base(min_var(1.),
                                         max_var(2.8),
                                         interval_size(2.),
                                         test_function_max(0.6168917686383567),
                                         py_expr("3*(1 - 3*x/2. + 2*x*x - x*x*x/2.)/10."))
    {}

    virtual double theTestFunctionValue(double x) const { return 3*(1 - 3*x/2. + 2*x*x - x*x*x/2.)/10.; }
};


struct QuarticPolynomial : public testFunction_base
{
    QuarticPolynomial(): testFunction_base(min_var(-1.0),
                                           max_var(1.0),
                                           interval_size(2.),
                                           test_function_max(0.2),
                                           py_expr("(x*x*x*x-0.8*x*x)/0.171964"))
    {}

    virtual double theTestFunctionValue(double x) const { return pow(x,4)-0.8*x*x; }
};


struct Exponential : public testFunction_base
{
    Exponential(): testFunction_base(min_var(1.),
                                     max_var(2.8),
                                     interval_size(2.),
                                     test_function_max(3.1),
                                     py_expr("np.exp(-3*x)/(-1 + np.exp(6.))*3.*np.exp(9.)"))
    {}

    virtual double theTestFunctionValue(double x) const { return exp(-3*x)/(-1 + exp(6.))*3*exp(9.); }
};

struct Cosine : public testFunction_base
{
    Cosine(): testFunction_base(min_var(1.),
                                max_var(PI+0.6),
                                interval_size(PI),
                                test_function_max(1.1/PI),
                                py_expr("(10.+np.cos(x*10.))/10./np.pi"))
    {}
    virtual double theTestFunctionValue(double x) const { return (10+cos(x*10.))/10./PI; }
};


// special case: don't use rejection sampling for this one
struct GaussianFromGSL : public testFunction_base
{
    GaussianFromGSL() : testFunction_base(min_var(-5),
                                          max_var(5.),
                                          interval_size(10.),
                                          test_function_max(0.5),
                                          py_expr(""))
    {}

    virtual double theTestFunctionValue(double x) const
    {
        throw std::logic_error("GaussianFromGSL::theTestFunctionValue() is not supposed to be called");
    }

    // this is more complex Python expression, we have to redefine the function from the base
    virtual std::string getExpr() const {
        return
            "def fn(x):\n"
            "    def gauss(m,s): return np.exp(-(x-m)**2/(2*s))/np.sqrt(2*np.pi*s)\n"
            "    return 0.2*gauss(0.,0.2) + 0.4*(gauss(2.,1.)+gauss(-2.,1))\n";
    }
};

static std::ostream& print_help(const char* argv0, std::ostream& strm) {
        strm
            << "Program to generate test data\n"
            << "Usage:\n"
            << "1) " << argv0 << " param_file.ini\n"
            << "Generates the histogram according to the parameters in `param_file.ini`\n"
            << "OR\n"
            << "2) " << argv0 << " -python n\n"
            << "(where n is a function number, 0-4) prints the Python expression corresponding to the generated function."
            << std::endl;

        return strm;
}
            


int Main(int argc, char **argv) {

        if (argc<2 || argc>3) {
            print_help(argv[0], std::cerr);
            return BAD_ARGS;
        }

        CubicPolynomial cubic_poly;
        QuarticPolynomial quartic_poly;
        Exponential exponential;
        Cosine cosine;
        GaussianFromGSL gauss;
        
        const testFunction_base* functions[]=
            {
                &cubic_poly,
                &quartic_poly,
                &exponential,
                &cosine,
                &gauss // must be the last!
            };
        const int nFunctions=sizeof(functions)/sizeof(*functions);
        const int nGaussFunctionChoice=nFunctions-1; // because gauss is the last
        
        if (std::string(argv[1])=="-python") {
            int fn=-1;
            if (argc!=3 ||
                sscanf(argv[2],"%d",&fn)!=1 ||
                fn<0 || fn>=nFunctions) {
                
                print_help(argv[0], std::cerr);
                return BAD_ARGS;
            }
            cout << functions[fn]->getExpr() << endl;
            return OK;
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
        LOGGER_VERBOSITY(verbose);
        
	long unsigned int seed=par.get(":RANDOMSEED", 0);
	if(seed==0) seed=time(NULL);
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=par.get(":SAMPLESIZE", 1e4);
	if(samplingSteps<defaultMinNumberTimesSampled) {cerr << "Too few sampling steps" << endl; return BAD_ARGS;}

	int theFunctionChoice=par.get(":FUNCTION",-1);
	if(theFunctionChoice<0 || theFunctionChoice>=nFunctions) {
            cerr << "WARNING: No such function, will instead use function 0 (cubic polynomial)" << endl;
            theFunctionChoice=0;
        }

	const testFunction_base& myTestFunction=*(functions[theFunctionChoice]);

	LOGGER << "Generating example histogram for function #" << theFunctionChoice << ":" << myTestFunction.getExpr();
       
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
		if(theFunctionChoice!=nGaussFunctionChoice)
			{
			bool accept=false;
			while(accept==false)
				{
				variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
				random=gsl_rng_uniform(RNG)*testFunctionMax;
				if(random<abs(myTestFunction.theTestFunctionValue(variable))) accept=true;
				}
			}
		else // special case: we do not use rejection sampling here
			{
			if(i%5==0) {variable = gsl_ran_gaussian (RNG, 0.2);}
			else if(i%2==0) {variable = 2+ gsl_ran_gaussian (RNG, 1.0);}
			else {variable = -2+ gsl_ran_gaussian (RNG, 1.0);}
			}
		
		binHistogram.sampleUniform(variable,whatsign(myTestFunction.theTestFunctionValue(variable)));
		}

        output << setprecision(12) << binHistogram;
        
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
