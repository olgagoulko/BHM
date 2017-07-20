#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	//return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;		//max 0.6168917686383567
	//return (2+variable)/8.;
	//return (10+cos(variable*10.))/10./PI;			//max 1.1/PI, maxVar=PI+0.6
	return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	//return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Evolution of BHM fit with MC time ----------------" << endl;
	
	long unsigned int seed=956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e3;
	
	double variable, random;
	double testFunctionMax=3.1;
	double intervalSize=2;
	double minVar=1.; double maxVar=2.8;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	double printStep=0.01; int numSteps=int((maxVar-minVar)/printStep)+1;
	pair<double,double> basisResult;
// 	slotBounds bounds1(-1,-0.5);
// 	slotBounds bounds2(-0.5,0);
// 	slotBounds bounds3(0, 0.5);
// 	slotBounds bounds4(0.5, 1);
	slotBounds bounds1(1, 2.8);
	//slotBounds bounds2(1.9, 2.8);
	int bf=4;
	double currentVar;
	int numberOfGoodFits=0;
	
	vector<basisSlot*> basisVector;
	basisVector.push_back(new taylorSlot(bounds1,bf)); 
	//basisVector.push_back(new taylorSlot(bounds2,bf)); 
	//basisVector.push_back(new taylorSlot(bounds3,bf)); 
	//basisVector.push_back(new taylorSlot(bounds4,bf)); 
	
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	histogramBasis basisHistogram(basisVector);
	
	ofstream output("histogram_testoutput.dat");
	
	int maxRounds=1000;
	long samplingStepsCounter=0;
	for(int round=0;round<maxRounds;round++)
		{
		for(int i=0; i<samplingSteps;i++)
			{
			bool accept=false;
			while(accept==false)
				{
				variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
				random=gsl_rng_uniform(RNG)*testFunctionMax;
				if(random<abs(testFunction(variable))) accept=true;
				}

			binHistogram.sampleUniform(variable,whatsign(testFunction(variable)));
			basisHistogram.sample(variable,whatsign(testFunction(variable)));
			samplingStepsCounter++;
			}
			
		histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingStepsCounter);
			
		double threshold=2;
		splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingStepsCounter, threshold, 0);
		bool acceptance=testBHMfit.getAcceptance();
			
		while(acceptance==false)
			{
			threshold+=0.5;
			testBHMfit = binHistogram.BHMfit(4, 2, samplingStepsCounter, threshold, 0);
			acceptance=testBHMfit.getAcceptance();
			if(threshold>5) break;
			}
			
		splineArray testJumpSuppression = binHistogram.BHMfit(4, 2, samplingStepsCounter, threshold, 1);
			
		if( acceptance ) numberOfGoodFits++;
			
		for(int i=0; i<numSteps;i++) 
			{
			currentVar=minVar+i*printStep;
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
				
			output << currentVar << '\t' << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar) << '\t'
			<< testJumpSuppression.splineValue(currentVar) << '\t' << testJumpSuppression.splineError(currentVar) << '\t' 
			<< basisResult.first << '\t' << basisResult.second << endl;
			}
		output << endl; output << endl;
			
		cout << "fit info: " << testBHMfit.getAcceptance() << '\t' << threshold << '\t' << testBHMfit.numberKnots() << endl;
		}

	cout << "number of good fits: " << numberOfGoodFits << endl;
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}