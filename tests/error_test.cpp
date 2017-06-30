#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	//return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;		//max 0.6168917686383567
	//return (2+variable)/8.;
	return (10+cos(variable*10.))/10./PI;			//max 1.1/pi, maxVar=Pi+1
	//return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	//return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Production test: full default routine ----------------" << endl;
	
	long unsigned int seed=76867;//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	
	double variable, random;
	double testFunctionMax=1.1/PI;
	double intervalSize=PI;
	double minVar=1.; double maxVar=PI+0.6;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	pair<double,double> basisResult;
	double currentVar;
	
	//double outputPoints[]={-1,-0.95,-0.5,-0.25,0,0.25,0.5,0.75,0.95,0.999};
	//double outputPoints[]={1,1.1,1.2,1.45,1.5,1.9,2.0,2.5,2.7,2.799};
	double outputPoints[]={1,1.2,1.5,1.9,2.0,2.5,2.7,3.0,3.2,3.5,3.74};
	
	ofstream output("histogram_testoutput.dat");
	
	for(int round=0;round<100;round++)
		{
		vector<basisSlot*> basisVector;
		vector<slotBounds> intervalBounds;
		//int bf=5;
		//slotBounds bounds1(-1, 1); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	
		int bf=4;
		//slotBounds bounds1(-1,-0.5); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
		//slotBounds bounds2(-0.5,0); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);
		//slotBounds bounds3(0, 0.5); taylorSlot slot3(bounds3,bf); basisVector.push_back(&slot3); intervalBounds.push_back(bounds3);
		//slotBounds bounds4(0.5, 1); taylorSlot slot4(bounds4,bf); basisVector.push_back(&slot4); intervalBounds.push_back(bounds4);
		
		slotBounds bounds1(1, 1.9); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
		//slotBounds bounds1b(1.45, 1.9); taylorSlot slot1b(bounds1b,bf); basisVector.push_back(&slot1b); intervalBounds.push_back(bounds1b);
		slotBounds bounds2(1.9, 3.8); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);
		
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis binHistogram(histogramVector);
		histogramBasis basisHistogram(basisVector);

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
		
		histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);
		
		double threshold=2;
		splineArray testNewRoutine2 = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
		bool acceptance=testNewRoutine2.getAcceptance();
		
		while(acceptance==false)
			{
			threshold+=0.5;
			testNewRoutine2 = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
			acceptance=testNewRoutine2.getAcceptance();
			if(threshold>4) break;
			}
		
		splineArray testcurvature = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 1);
		
		int numberSteps=10000; double stepWidth=(maxVar-minVar)/double(numberSteps);
		double differenceIntegral=0; double differenceCurvatureIntegral=0; double differenceBasisIntegral=0;
		int coverage[6]; for(int i=0; i<6;i++) coverage[i]=0;
		for(int i=0; i<numberSteps;i++) 
			{
			currentVar=minVar+(i+0.5)*stepWidth;
		
			//double difference = abs(testNewRoutine2.splineValue(currentVar)-testFunction(currentVar)/0.171964);
			double difference = abs(testNewRoutine2.splineValue(currentVar)-testFunction(currentVar));
			differenceIntegral+=difference/testNewRoutine2.splineError(currentVar);
			if(difference<testNewRoutine2.splineError(currentVar)) coverage[0]++;
			if(difference<2*testNewRoutine2.splineError(currentVar)) coverage[1]++;
			
			difference = abs(testcurvature.splineValue(currentVar)-testFunction(currentVar));//difference = abs(testcurvature.splineValue(currentVar)-testFunction(currentVar)/0.171964);
			differenceCurvatureIntegral+=difference/testcurvature.splineError(currentVar);
			if(difference<testcurvature.splineError(currentVar)) coverage[2]++;
			if(difference<2*testcurvature.splineError(currentVar)) coverage[3]++;
		
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
			difference = abs(basisResult.first-testFunction(currentVar)); //difference = abs(basisResult.first-testFunction(currentVar)/0.171964);
			differenceBasisIntegral+=difference/basisResult.second;
			if(difference<basisResult.second) coverage[4]++;
			if(difference<2*basisResult.second) coverage[5]++;
			}
		
		for(int i=0; i<10;i++) 
			{
			currentVar=outputPoints[i];
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
			output << testNewRoutine2.splineValue(currentVar) << '\t' << testNewRoutine2.splineError(currentVar) << '\t'
			<< testcurvature.splineValue(currentVar) << '\t' << testcurvature.splineError(currentVar) << '\t' 
			<< basisResult.first << '\t' << basisResult.second << '\t';
			}
		output << differenceIntegral*stepWidth << '\t' << differenceCurvatureIntegral*stepWidth << '\t' << differenceBasisIntegral*stepWidth << '\t';
		for(int i=0; i<6;i++) output << coverage[i]/double(numberSteps) << '\t';
		output << testNewRoutine2.getAcceptance() << '\t' << testcurvature.getAcceptance() << '\t' << threshold << '\t' << testNewRoutine2.numberKnots() << endl;
		}
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}