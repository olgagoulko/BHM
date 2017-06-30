#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	//return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;		//max 0.6168917686383567
	//return (2+variable)/8.;
	//return (10+cos(variable*10.))/10./PI;			//max 1.1/PI, maxVar=PI+0.6, intervalSize=PI;
	return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	//return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Production test: full default routine ----------------" << endl;
	
	long unsigned int seed=6482555;//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	
	double variable, random;
	double testFunctionMax=3.1;
	double intervalSize=2.0;
	double minVar=1.; double maxVar=2.8;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	double printStep=0.02; int numSteps=int((maxVar-minVar)/printStep);
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	double currentVar;
	histogramBasis binHistogram(histogramVector);
	
	ofstream output("histogram_testoutput.dat");
	
	/*for(int round=0;round<57;round++) for(int i=0; i<samplingSteps;i++) {bool accept=false;
		while(accept==false)
		{
			variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
			random=gsl_rng_uniform(RNG)*testFunctionMax;
			if(random<abs(testFunction(variable))) accept=true;
		}}*/

	for(int i=0; i<samplingSteps;i++)
		{
		bool accept=false;
		while(accept==false)
			{
			variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
			random=gsl_rng_uniform(RNG)*testFunctionMax;
			if(random<abs(testFunction(variable))) accept=true;
			}
		binHistogram.sample(variable,whatsign(testFunction(variable)));
		}
		
	int goodElementaryBins=0; for(unsigned int i=0;i<pow(2.,10);i++) {histogramVector[i] -> updateEnoughSampled(); goodElementaryBins+=histogramVector[i] -> enoughSampled();}
	cout << "Good elementary bins = " << goodElementaryBins << endl; cout << endl;
			
	double threshold=0;
	for(int round=0;round<30;round++)
		{
		splineArray testNewRoutine2 = binHistogram.splineProcedure2(4, 2, samplingSteps, threshold, 0);
		cout << "new procedure: " << testNewRoutine2.getAcceptance() << '\t' << threshold << '\t' << testNewRoutine2.numberKnots() << endl;
		//cout << endl; testNewRoutine2.printSplineArrayInfo(); cout << endl;
		//testNewRoutine2.printSplines(); cout << endl;
		//bool acceptance=testNewRoutine2.getAcceptance();
		splineArray testcurvature = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
		cout << "old procedure: " << testcurvature.getAcceptance() << '\t' << threshold << '\t' << testcurvature.numberKnots() << endl;
		//cout << endl; testcurvature.printSplineArrayInfo(); cout << endl;
		//testcurvature.printSplines(); cout << endl;
		
		for(int i=0; i<numSteps;i++) 
			{
			currentVar=minVar+i*printStep+printStep/2;

			output << currentVar << '\t' << testNewRoutine2.splineValue(currentVar) << '\t' << testNewRoutine2.splineError(currentVar) << '\t'
			<< testcurvature.splineValue(currentVar) << '\t' << testcurvature.splineError(currentVar) << " 0 0" << endl;
			}
		output << endl; output << endl;
				
		threshold+=0.5;
		}
		
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}