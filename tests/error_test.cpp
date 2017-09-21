/*** LICENCE: ***
Bin histogram method for restoration of smooth functions from noisy integrals. Copyright (C) 2017 Olga Goulko

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.

*** END OF LICENCE ***/
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
	
	cout << "----------------- Robust error analysis using repeated independent runs ----------------" << endl;
	
	long unsigned int seed=956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e6;
	
	double variable, random;
	double testFunctionMax=3.1;
	double intervalSize=2;
	double minVar=1.; double maxVar=2.8;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	pair<double,double> basisResult;
	double currentVar;
	
	//double outputPoints[]={-1,-0.95,-0.5,-0.25,0,0.25,0.5,0.75,0.95,0.999};
	double outputPoints[]={1,1.1,1.2,1.45,1.5,1.9,2.0,2.5,2.7,2.799};
	//double outputPoints[]={1,1.2,1.5,2.0,2.5,2.7,2.9,3.2,3.5,3.74};
	
	ofstream output("histogram_testoutput.dat");
	
	int maxRounds=100;
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> basisVector;
		//int bf=5;
		//slotBounds bounds1(-1, 1); basisVector.push_back(new taylorSlot(bounds1,bf));
	
		int bf=4;
		//slotBounds bounds1(-1,-0.5); basisVector.push_back(new taylorSlot(bounds1,bf));
		//slotBounds bounds2(-0.5,0); basisVector.push_back(new taylorSlot(bounds2,bf));
		//slotBounds bounds3(0, 0.5); basisVector.push_back(new taylorSlot(bounds3,bf));
		//slotBounds bounds4(0.5, 1); basisVector.push_back(new taylorSlot(bounds4,bf));
		
		slotBounds bounds1(1, 2.8); basisVector.push_back(new taylorSlot(bounds1,bf));
		//slotBounds bounds1b(1.45, 1.9); basisVector.push_back(new taylorSlot(bounds1b,bf));
		//slotBounds bounds2(1.9, 2.8); basisVector.push_back(new taylorSlot(bounds2,bf));
		
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
			binHistogram.sampleUniform(variable,whatsign(testFunction(variable)));
			basisHistogram.sample(variable,whatsign(testFunction(variable)));
			}
		
		histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);
		
		fitAcceptanceThreshold threshold; threshold.min=2; threshold.max=5.5; threshold.steps=7;
		splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0, true, false);
		bool acceptance=testBHMfit.getAcceptance();
		
		cout << "round = " << round << ", acceptance = " << acceptance << ", number knots = " << testBHMfit.numberKnots() << endl;

		splineArray testJumpSuppression = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 1, true, false);
		
		int numberSteps=10000; double stepWidth=(maxVar-minVar)/double(numberSteps);
		double differenceIntegral=0; double differenceCurvatureIntegral=0; double differenceBasisIntegral=0;
		int coverage[6]; for(int i=0; i<6;i++) coverage[i]=0;
		for(int i=0; i<numberSteps;i++) 
			{
			currentVar=minVar+(i+0.5)*stepWidth;
		
			double difference = abs(testBHMfit.splineValue(currentVar)-testFunction(currentVar));
			differenceIntegral+=difference/testBHMfit.splineError(currentVar);
			if(difference<testBHMfit.splineError(currentVar)) coverage[0]++;
			if(difference<2*testBHMfit.splineError(currentVar)) coverage[1]++;
			
			/**/difference = abs(testJumpSuppression.splineValue(currentVar)-testFunction(currentVar));
			differenceCurvatureIntegral+=difference/testJumpSuppression.splineError(currentVar);
			if(difference<testJumpSuppression.splineError(currentVar)) coverage[2]++;
			if(difference<2*testJumpSuppression.splineError(currentVar)) coverage[3]++;/**/
		
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
			difference = abs(basisResult.first-testFunction(currentVar));
			differenceBasisIntegral+=difference/basisResult.second;
			if(difference<basisResult.second) coverage[4]++;
			if(difference<2*basisResult.second) coverage[5]++;
			}
		
		for(int i=0; i<10;i++) 
			{
			currentVar=outputPoints[i];
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
			output << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar) << '\t'
			<< testJumpSuppression.splineValue(currentVar) << '\t' << testJumpSuppression.splineError(currentVar) << '\t' 
			<< basisResult.first << '\t' << basisResult.second << '\t';
			}
		output << differenceIntegral*stepWidth << '\t' << differenceCurvatureIntegral*stepWidth << '\t' << differenceBasisIntegral*stepWidth << '\t';
		for(int i=0; i<6;i++) output << coverage[i]/double(numberSteps) << '\t';
		output << testBHMfit.getAcceptance() << '\t' << testJumpSuppression.getAcceptance() << '\t' << testBHMfit.getThreshold() << '\t' << testBHMfit.numberKnots() << endl;
		
		}
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}