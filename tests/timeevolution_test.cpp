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
	
	long unsigned int seed=956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e4;
	
	double variable, random;
	double testFunctionMax=0.2;
	double intervalSize=2;
	double minVar=-1.; double maxVar=1.;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	double printStep=0.02; int numSteps=int((maxVar-minVar)/printStep);
	pair<double,double> basisResult;
	slotBounds bounds1(-1,-0.5);
	slotBounds bounds2(-0.5,0);
	slotBounds bounds3(0, 0.5);
	slotBounds bounds4(0.5, 1);
	vector<slotBounds> intervalBounds;
	intervalBounds.push_back(bounds1); intervalBounds.push_back(bounds2); intervalBounds.push_back(bounds3); intervalBounds.push_back(bounds4);
	int bf=4;
	double currentVar;
	int numberOfGoodFits=0;
	
	vector<histogramBasis> binHistogramVector; vector<histogramBasis> basisHistogramVector;
	
	ofstream output("histogram_testoutput.dat");
	
	int maxRounds=100;
	int samplingStepsCounter[maxRounds];
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> basisVector;
		samplingStepsCounter[round]=0;
	
		basisVector.push_back(new taylorSlot(bounds1,bf)); 
		basisVector.push_back(new taylorSlot(bounds2,bf)); 
		basisVector.push_back(new taylorSlot(bounds3,bf)); 
		basisVector.push_back(new taylorSlot(bounds4,bf)); 
		
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis binHistogram(histogramVector);
		histogramBasis basisHistogram(basisVector);
		
		binHistogramVector.push_back(binHistogram); basisHistogramVector.push_back(basisHistogram);
		}

	for(int i=0; i<samplingSteps*maxRounds;i++)
		{
		bool accept=false;
		while(accept==false)
			{
			variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
			random=gsl_rng_uniform(RNG)*testFunctionMax;
			if(random<abs(testFunction(variable))) accept=true;
			}
		for(int round=0;round<maxRounds;round++)
			{
			if(i<(round+1)*samplingSteps)
				{
				binHistogramVector[round].sample(variable,whatsign(testFunction(variable)));
				basisHistogramVector[round].sample(variable,whatsign(testFunction(variable)));
				samplingStepsCounter[round]++;
				}
			}
		}

	for(int round=0;round<maxRounds;round++)
		{
		histogramBasis scaledBasisHistogram = basisHistogramVector[round].scaledHistogram(samplingStepsCounter[round]);
		
		double threshold=2;
		splineArray testNewRoutine2 = binHistogramVector[round].splineProcedure(4, 2, samplingStepsCounter[round], threshold, 0);
		bool acceptance=testNewRoutine2.getAcceptance();

		while(acceptance==false)
			{
			threshold+=0.5;
			testNewRoutine2 = binHistogramVector[round].splineProcedure(4, 2, samplingStepsCounter[round], threshold, 0);
			acceptance=testNewRoutine2.getAcceptance();
			if(threshold>5) break;
			}
		
		splineArray testcurvature = binHistogramVector[round].splineProcedure(4, 2, samplingStepsCounter[round], threshold, 1);
		
		if( acceptance ) numberOfGoodFits++;

		for(int i=0; i<numSteps;i++) 
			{
			currentVar=minVar+i*printStep+printStep/2;
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);

			output << currentVar << '\t' << testNewRoutine2.splineValue(currentVar) << '\t' << testNewRoutine2.splineError(currentVar) << '\t'
			<< testcurvature.splineValue(currentVar) << '\t' << testcurvature.splineError(currentVar) << '\t' 
			<< basisResult.first << '\t' << basisResult.second << endl;
			}
		output << endl; output << endl;
		
		cout << "fit info: " << testNewRoutine2.getAcceptance() << '\t' << threshold << '\t' << testNewRoutine2.numberKnots() << endl;
		}
		
	cout << "number of good fits: " << numberOfGoodFits << endl;
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}