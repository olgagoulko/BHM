#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

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



int main(int argc, char **argv) {
	
	cout << "----------------- Example BHM code ----------------" << endl;
	
	long unsigned int seed=956475;//time(NULL);
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	//adjust parameters ---------------------------------------------
	long samplingSteps=1e4;
	unsigned int splineOrder=4;
	unsigned int minLevel=2;
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
	
	cout << endl; cout << "BHM fit:" << endl;
	double threshold=2; double jumpSuppression=0;
	splineArray testBHMfit = binHistogram.BHMfit(splineOrder, minLevel, samplingSteps, threshold, jumpSuppression);
	cout << endl;
	cout << "acceptable fit = " << testBHMfit.getAcceptance() << endl; cout << endl;
	cout << "Printing spline:" << endl;
	testBHMfit.printSplineArrayInfo(); cout << endl;
	testBHMfit.printSplines(); 
	cout << endl;
	
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
	
	return 0;
	
}