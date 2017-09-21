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
	//return (10+cos(variable*10.))/10./PI;			//max 1.1/PI, maxVar=PI+0.6, intervalSize=PI;
	return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	//return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Testing effect of goodness-of-fit threshold ----------------" << endl;
	
	long unsigned int seed=956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	
	double variable, random;
	double testFunctionMax=3.1;
	double intervalSize=2.0;
	double minVar=1.; double maxVar=2.8;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	double printStep=0.01; int numSteps=int((maxVar-minVar)/printStep);
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	
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
		}
		
	int goodElementaryBins=0; for(unsigned int i=0;i<pow(2.,10);i++) {histogramVector[i] -> updateEnoughSampled(); goodElementaryBins+=histogramVector[i] -> enoughSampled();}
	cout << "Good elementary bins = " << goodElementaryBins << endl; cout << endl;
	
	ofstream output("histogram_testoutput.dat");
			
	fitAcceptanceThreshold threshold;
	threshold.min=0; threshold.max=0; threshold.steps=0;
	for(int round=0;round<30;round++)
		{
		splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0, true, false);
		cout << "BHM fit: " << testBHMfit.getAcceptance() << '\t' << threshold.min << '\t' << testBHMfit.numberKnots() << endl;
		splineArray testJumpSuppression = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 1, true, false);
		cout << "Jump suppression: " << testJumpSuppression.getAcceptance() << '\t' << threshold.min << '\t' << testJumpSuppression.numberKnots() << endl;
		
		for(int i=0; i<numSteps;i++) 
			{
			variable=minVar+i*printStep;

			output << variable << '\t' << testBHMfit.splineValue(variable) << '\t' << testBHMfit.splineError(variable) << '\t'
			<< testJumpSuppression.splineValue(variable) << '\t' << testJumpSuppression.splineError(variable) << endl;
			}
		output << endl; output << endl;
				
		threshold.min+=0.5;
		}
		
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}