#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	//return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;		//max 0.6168917686383567
	//return (2+variable)/8.;
	//return (10+cos(variable*10.))/10./PI;			//max 1.1/pi, maxVar=Pi+1
	//return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Production test: full default routine ----------------" << endl;
	
	long unsigned int seed=324241209876;//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e4;
	double variable;
	
	vector<basisSlot*> basisVector; int bf=5;
	vector<slotBounds> intervalBounds;
	slotBounds bounds1(-1, 1); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	
	//slotBounds bounds1(1, 1.9); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	//slotBounds bounds2(1.9, 2.8); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);
	/*slotBounds bounds3(1.9, 2.35); taylorSlot slot3(bounds3,bf); basisVector.push_back(&slot3); intervalBounds.push_back(bounds3);
	slotBounds bounds4(2.35, 2.8); taylorSlot slot4(bounds4,bf); basisVector.push_back(&slot4); intervalBounds.push_back(bounds4);*/
	
	/*slotBounds bounds1(-1,-0.5); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	slotBounds bounds2(-0.5,0); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);
	slotBounds bounds3(0, 0.5); taylorSlot slot3(bounds3,bf); basisVector.push_back(&slot3); intervalBounds.push_back(bounds3);
	slotBounds bounds4(0.5, 1); taylorSlot slot4(bounds4,bf); basisVector.push_back(&slot4); intervalBounds.push_back(bounds4);*/
	
	/*slotBounds bounds1(-1,-0.75); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1);
	slotBounds bounds2(-0.75,-0.5); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2);
	slotBounds bounds3(-0.5,-0.25); taylorSlot slot3(bounds3,bf); basisVector.push_back(&slot3);
	slotBounds bounds4(-0.25,0); taylorSlot slot4(bounds4,bf); basisVector.push_back(&slot4);
	slotBounds bounds5(0,0.25); taylorSlot slot5(bounds5,bf); basisVector.push_back(&slot5);
	slotBounds bounds6(0.25,0.5); taylorSlot slot6(bounds6,bf); basisVector.push_back(&slot6);
	slotBounds bounds7(0.5,0.75); taylorSlot slot7(bounds7,bf); basisVector.push_back(&slot7);
	slotBounds bounds8(0.75,1); taylorSlot slot8(bounds8,bf); basisVector.push_back(&slot8);*/
	
	//initialize histogram and sample function
	//double minVar=1.; double maxVar=2.8;
	//double minVar=1; double maxVar=PI+0.6;
	double minVar=-1.; double maxVar=1.0;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	histogramBasis basisHistogram(basisVector);

	double random;
	double testFunctionMax=0.2;
	//double intervalSize=3.0-minVar;
	//double intervalSize=PI;
	double intervalSize=2;
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
		basisHistogram.sample(variable,whatsign(testFunction(variable)));
		}
		
	//make coarser histogram so that each bin has enough data
	int minNumberTimesSampled = 10;
	histogramBasis analysisHistogram = binHistogram.coarseGrainedHistogram(minNumberTimesSampled);
	cout << "Slots in coarse grained analysis histogram: " << analysisHistogram.getSize() << endl;
	
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledAnalysisHistogram = analysisHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);
	//cout << "Basis projection coefficients" << endl;
	//scaledBasisHistogram.getSlot(0) -> printSampledCoeffs();
	//cout << endl;
	
	splineArray testNewRoutine2 = binHistogram.splineProcedure(4, 2, samplingSteps, 2, 0);
	cout << "New spline procedure more big bins" << endl;
	testNewRoutine2.printSplineArrayInfo(); cout << endl;
	testNewRoutine2.printSplines();
	
	splineArray testcurvature = binHistogram.splineProcedure(4, 2, samplingSteps, 2, 1);
	cout << "Suppress curvature" << endl;
	testcurvature.printSplineArrayInfo(); cout << endl;
	testcurvature.printSplines();
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.02; double currentVar; double currentVar2; pair<double,double> coarseGrainedResult; pair<double,double> basisResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
	{
		currentVar=minVar+i*printStep+printStep/2;
		currentVar2=currentVar-printStep/2;
		coarseGrainedResult=scaledAnalysisHistogram.sampledFunctionValueWeightedAverage(currentVar);
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar2);
		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		<< '\t' << testNewRoutine2.splineValue(currentVar) << '\t' << testNewRoutine2.splineError(currentVar)
		<< '\t' << testNewRoutine2.splineDerivative(currentVar,1) << '\t' << testNewRoutine2.splineDerivative(currentVar,2) << '\t' << testNewRoutine2.splineDerivative(currentVar,3)
		<< '\t' << testcurvature.splineValue(currentVar) << '\t' << testcurvature.splineError(currentVar)
		<< '\t' << testcurvature.splineDerivative(currentVar,1) << '\t' << testcurvature.splineDerivative(currentVar,2) << '\t' << testcurvature.splineDerivative(currentVar,3)
		<< '\t' << currentVar2 << '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
	}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}