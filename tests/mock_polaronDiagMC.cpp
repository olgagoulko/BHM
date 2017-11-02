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

double valueToSample(double variable, int sector)
	{
	double value = 1.0;
	if(sector>0) {value*=exp(-variable); if(sector==1) value*=5; else if(sector==2) value*=-3;}	
	return value;
	}
	
double weight(double variable, int sector)
	{
	double theweight=1.;
	if(sector==0) theweight=1.5;
	else if(sector==1) theweight=exp(0.01*variable);
	else if(sector==3) theweight=exp(-variable);
	return theweight;
	}
	
bool update_insert(gsl_rng* RNG, int &currentSector, double &currentVar, double intervalSize)
	{
	bool accept=false; double newVar; int newSector;
	if(currentSector==0)
		{
		newVar=gsl_rng_uniform(RNG)*intervalSize;
		newSector=1+gsl_rng_uniform_int(RNG,3);
		double newWeight=weight(newVar,newSector);
		accept=acceptreject(3*intervalSize*newWeight/weight(0,0), RNG);
		}
	if(accept) {currentSector=newSector; currentVar=newVar;}
	return accept;
	}
	
bool update_delete(gsl_rng* RNG, int &currentSector, double &currentVar, double intervalSize)
	{
	bool accept=false;
	if(currentSector!=0)
		{
		accept=acceptreject(weight(0,0)/weight(currentVar,currentSector)/3./intervalSize, RNG);
		}
	if(accept) {currentSector=0; currentVar=0;}
	return accept;
	}
	
bool update_changevariable(gsl_rng* RNG, int currentSector, double &currentVar, double intervalSize)
	{
	bool accept=false; double newVar;
	if(currentSector!=0)
		{
		newVar=gsl_rng_uniform(RNG)*intervalSize;
		accept=acceptreject(weight(newVar,currentSector)/weight(currentVar,currentSector), RNG);
		}
	if(accept) {currentVar=newVar;}
	return accept;
	}
	
bool update_switchsector(gsl_rng* RNG, int &currentSector, double currentVar)
	{
	bool accept=false; int newSector=currentSector;
	if(currentSector!=0)
		{
		while(newSector==currentSector) newSector=1+gsl_rng_uniform_int(RNG,3);
		accept=acceptreject(weight(currentVar,newSector)/weight(currentVar,currentSector), RNG);
		}
	if(accept) {currentSector=newSector;}
	return accept;
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Very simplified version of polaron code ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e7;
	long MCsteps=1e8;
	
	vector<basisSlot*> basisVector; int bf=4;
	vector<slotBounds> intervalBounds;
	slotBounds bounds1(0,1); basisVector.push_back(new taylorSlot(bounds1,bf));
	slotBounds bounds2(1,3); basisVector.push_back(new taylorSlot(bounds2,bf));
	histogramBasis basisHistogram(basisVector);
	
	double minVar=0.; double maxVar=3.0; double maxVarTotal=3.0;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);

	double currentVar=0; double tosample;
	int currentSector=0; int norm=0;
	double intervalSize=maxVarTotal-minVar;
	int whichUpdate;
	int sectorCounter[4]; sectorCounter[0]=0; sectorCounter[1]=0; sectorCounter[2]=0; sectorCounter[3]=0;
	int updateCounter[4][2]; for(int i=0;i<4;i++) for(int j=0;j<2;j++) updateCounter[i][j]=0;
	for(int i=0; i<MCsteps;i++)
		{
		whichUpdate=gsl_rng_uniform_int(RNG,4);
		if(whichUpdate==0)
			{
			updateCounter[0][0]++;
			updateCounter[0][1]+=update_insert(RNG, currentSector, currentVar, intervalSize);
			}
		else if(whichUpdate==1)
			{
			updateCounter[1][0]++;
			updateCounter[1][1]+=update_delete(RNG, currentSector, currentVar, intervalSize);
			}
		else if(whichUpdate==2)
			{
			updateCounter[2][0]++;
			updateCounter[2][1]+=update_changevariable(RNG, currentSector, currentVar, intervalSize);
			}
		else if(whichUpdate==3)
			{
			updateCounter[3][0]++;
			updateCounter[3][1]+=update_switchsector(RNG, currentSector, currentVar);
			}

		sectorCounter[currentSector]++;
		
		//if(i%10==0)
			//{
			if(currentSector==0) norm++;
			else if(currentSector<3)
				{
				tosample=valueToSample(currentVar,currentSector);
				binHistogram.sampleUniform(currentVar,tosample); basisHistogram.sample(currentVar,tosample);
				}
			//}
		}
	
	cout << "Sector counters: " << sectorCounter[0]/double(samplingSteps) << '\t' << sectorCounter[1]/double(samplingSteps) << '\t' << sectorCounter[2]/double(samplingSteps) << '\t' << sectorCounter[3]/double(samplingSteps) << endl;
	
	cout << "Update acceptance ratios" << endl;
	for(int i=0;i<4;i++) cout << updateCounter[i][1]/double(updateCounter[i][0]) << endl;
	
	cout << "Norm = " << norm << endl;
	double fullnorm=norm/weight(0,0);
	cout << "Full norm = " << fullnorm << endl;
	cout << "Relative full norm = " << fullnorm/double(samplingSteps) << endl;
	
	ofstream bhmoutput("histogram.dat");
	bhmoutput << setprecision(12) << binHistogram;
	LOGGER_VERBOSITY(true);
	
	fitAcceptanceThreshold threshold; threshold.min=2; threshold.max=5; threshold.steps=4;
	BHMparameters theParameters;
	theParameters.dataPointsMin=100;
	theParameters.splineOrder=4;
	theParameters.minLevel=2;
	theParameters.threshold=threshold;
	theParameters.usableBinFraction=0.25;
	theParameters.jumpSuppression=0;
	
	splineArray testBHMfit = binHistogram.normalizedHistogram(fullnorm/double(samplingSteps)).BHMfit(theParameters, samplingSteps, true);
	testBHMfit.printSplineArrayInfo(cout); cout << endl;
	testBHMfit.printSplines(cout);
	
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; pair<double,double> coarseGrainedResult; pair<double,double> basisResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
	coarseGrainedResult=scaledBinHistogram.normalizedHistogram(fullnorm/double(samplingSteps)).sampledFunctionValueWeightedAverage(currentVar);
	basisResult=scaledBasisHistogram.normalizedHistogram(fullnorm/double(samplingSteps)).sampledFunctionValueWeightedAverage(currentVar);
	
		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		<< '\t' << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar)
		<< '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}

	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}