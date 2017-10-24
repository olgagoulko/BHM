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

double testFunction(double variable, int sector)
	{
	double coefficient=1; double prefactor=1;
	if(sector==1) coefficient=0.99;
	else if(sector==0) {coefficient=0; prefactor=1./3.;}
	
	return prefactor*exp(-variable*coefficient);
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Test with sign problem and Markov chain updates, normalization sector ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e7;
	long MCsteps=1e8;
	
	vector<basisSlot*> basisVector; int bf=4;
	vector<slotBounds> intervalBounds;
	slotBounds bounds1(0,3); basisVector.push_back(new taylorSlot(bounds1,bf));
	histogramBasis basisHistogram(basisVector);
	
	taylorSlot normSlot(bounds1,0);

	double minVar=0.; double maxVar=3.0; double maxVarTotal=3.0;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	
	/**/int maxRounds=100; int bootstrapSamples=100; vector<histogramBasis> binHistogramVector;
	vector< vector<double> > averageCoeffs; vector< vector<double> > averageCov;
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis binHistogram(histogramVector);
		binHistogramVector.push_back(binHistogram);
		}/**/

	double currentVar=0; double newVar;
	int currentSector=0;
	double currentWeight=testFunction(currentVar,currentSector); double newWeight;
	double intervalSize=maxVarTotal-minVar;
	int whichUpdate; bool accept;
	int sectorCounter[3]; sectorCounter[0]=0; sectorCounter[1]=0; sectorCounter[2]=0;
	for(int i=0; i<MCsteps;i++)
		{
		whichUpdate=gsl_rng_uniform_int(RNG,2);
		if(whichUpdate==0)
			{
			//switch between sectors
			int newSector=currentSector+1+gsl_rng_uniform_int(RNG,2);
			if(newSector>1) newSector-=3;
			newWeight=testFunction(currentVar,newSector);
			accept=acceptreject(newWeight/currentWeight, RNG);
			if(accept==true) {currentSector=newSector; currentWeight=newWeight;}
			}
		else
			{
			newVar=gsl_rng_uniform(RNG)*intervalSize;
			newWeight=testFunction(newVar,currentSector);
			accept=acceptreject(newWeight/currentWeight, RNG);
			if(accept==true)
				{
				currentVar=newVar;
				currentWeight=newWeight;
				}
			}
		if(i%10==0) 
			{
			binHistogram.sampleUniform(currentVar,currentSector);
			basisHistogram.sample(currentVar,currentSector);
		
			binHistogramVector[(i/10)%maxRounds].sampleUniform(currentVar,currentSector);
		
			sectorCounter[currentSector+1]++;
			if(currentSector==0) normSlot.sample(currentVar,1.);
			}
		}
	
	cout << "Sector counters: " << sectorCounter[0]/double(samplingSteps) << '\t' << sectorCounter[1]/double(samplingSteps) << '\t' << sectorCounter[2]/double(samplingSteps) << endl;
	
	cout << "Normalization slot " << endl;
	normSlot.printSampledCoeffs();
	cout << normSlot.sampledFunctionValue(0) << " " << normSlot.sampledFunctionError(0) << endl;
	
	normSlot.scale(samplingSteps);
	normSlot.printSampledCoeffs();
	cout << normSlot.sampledIntegral() << " " << normSlot.sampledIntegralError() << endl;
		
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	fitAcceptanceThreshold threshold; threshold.min=2;
	BHMparameters theParameters;
	theParameters.dataPointsMin=100;
	theParameters.splineOrder=4;
	theParameters.minLevel=2;
	theParameters.threshold=threshold;
	theParameters.usableBinFraction=0.25;
	theParameters.jumpSuppression=0;
	
	splineArray testBHMfit = binHistogram.normalizedHistogram(normSlot.sampledIntegral()).BHMfit(theParameters, samplingSteps, true);
	testBHMfit.printSplineArrayInfo(cout); cout << endl;
	testBHMfit.printSplines(cout);
	
	/**/intervalBounds=testBHMfit.getBounds();
	for(int round=0;round<bootstrapSamples;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis combinedHistogram(histogramVector);
		for(int i=0;i<maxRounds;i++)
			{
			unsigned int random=gsl_rng_uniform_int (RNG, maxRounds);
			combinedHistogram.addAnotherHistogram(binHistogramVector[random]);
			}

		vector< double > dummy;
		splineArray testBHMfit;

		for(unsigned int j=0;j<intervalBounds.size();j++)
			{
			vector<double> theCoeffVec(4, 0.0); vector<double> theCovVec(10, 0.0);
			averageCoeffs.push_back(theCoeffVec);
			averageCov.push_back(theCovVec);
			}

		//vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.binHierarchy(samplingSteps);
		vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.normalizedHistogram(normSlot.sampledIntegral()).binHierarchy(samplingSteps, theParameters.dataPointsMin, theParameters.usableBinFraction);
		testBHMfit = matchedSplineFit(currentAnalysisBins, intervalBounds, 4, 0, dummy, dummy,threshold.min);
			
		for(unsigned int j=0;j<intervalBounds.size();j++)
			{
			vector<double> currentCoeffs=testBHMfit.getSplinePiece(j) -> getCoefficients();
			vector<double> delta;
			for(unsigned int i=0;i<4;i++)
				{
				delta.push_back(currentCoeffs[i] - averageCoeffs[j][i]);
				averageCoeffs[j][i] += delta[i] / (round+1);
				}
				for(unsigned int i=0;i<4;i++)
					for(unsigned int k=i;k<4;k++)
					{
					averageCov[j][i*(7-i)/2+k] += delta[i]*delta[k]*round/(round+1);
					}
			}
		}
		
	vector<splinePiece*> theBootstrapSplines;
	for(unsigned int i=0;i<intervalBounds.size();i++) 
		{
		theBootstrapSplines.push_back(new splinePiece(intervalBounds[i]));
			
		vector<double> theErrorCoefficients; for(unsigned int j=0;j < 7;j++) theErrorCoefficients.push_back(0); 
		for(unsigned int j=0;j < 4;j++)
			for(unsigned int k=j;k < 4;k++)
				{
				theErrorCoefficients[j+k]+=averageCov[i][j*(7-j)/2+k]/(bootstrapSamples - 1);
				if(j!=k) theErrorCoefficients[j+k]+=averageCov[i][j*(7-j)/2+k]/(bootstrapSamples - 1);
				}
		theBootstrapSplines[i] -> setSplinePiece(averageCoeffs[i], theErrorCoefficients);
		}
		
	splineArray averageBootstrapSpline(theBootstrapSplines);
	cout << "Average Bootstrap Spline" << endl;
	averageBootstrapSpline.printSplineArrayInfo(cout); averageBootstrapSpline.printSplines(cout); cout << endl;
	/**/
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.05; pair<double,double> coarseGrainedResult; pair<double,double> basisResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
		coarseGrainedResult=scaledBinHistogram.normalizedHistogram(normSlot.sampledIntegral()).sampledFunctionValueWeightedAverage(currentVar);
		basisResult=scaledBasisHistogram.normalizedHistogram(normSlot.sampledIntegral()).sampledFunctionValueWeightedAverage(currentVar);
	
		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		<< '\t' << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar)
		<< '\t' << averageBootstrapSpline.splineValue(currentVar) << '\t' << averageBootstrapSpline.splineError(currentVar)
		<< '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}