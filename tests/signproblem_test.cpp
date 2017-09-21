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
	double coefficient=1;
	if(sector==1) coefficient=0.99;
	return exp(-variable*coefficient);
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Test with sign problem and Markov chain updates ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e7;
	long MCsteps=1e8;
	
	vector<basisSlot*> basisVector; int bf=4;
	vector<slotBounds> intervalBounds;
	slotBounds bounds1(0,3); basisVector.push_back(new taylorSlot(bounds1,bf));
	histogramBasis basisHistogram(basisVector);

	double minVar=0.; double maxVar=3.0; double maxVarTotal=3.0;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	
	int maxRounds=100; int bootstrapSamples=100; vector<histogramBasis> binHistogramVector;
	vector< vector<double> > averageCoeffs; vector< vector<double> > averageCov;
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis binHistogram(histogramVector);
		binHistogramVector.push_back(binHistogram);
		}

	double currentVar=0; double newVar;
	int currentSector=-1;
	double currentWeight=testFunction(currentVar,currentSector); double newWeight;
	double intervalSize=maxVarTotal-minVar;
	int whichUpdate; bool accept;
	//int sectorCounter;
	for(int i=0; i<MCsteps;i++)
		{
		whichUpdate=gsl_rng_uniform_int(RNG,2);
		if(whichUpdate==0)
			{
			//switch between branches
			newWeight=testFunction(currentVar,-currentSector);
			accept=acceptreject(newWeight/currentWeight, RNG);
			if(accept==true) {currentSector=-currentSector; currentWeight=newWeight;}
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
		
			binHistogramVector[(i/10)%maxRounds].sample(currentVar,currentSector);
			}
		//if(currentSector==1) sectorCounter++;
		}
		
	//cout << "Sector counter = " << sectorCounter/double(samplingSteps) << endl;
		
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	fitAcceptanceThreshold threshold; threshold.min=2;
	splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0, true, false);
	testBHMfit.printSplineArrayInfo(cout); cout << endl;
	testBHMfit.printSplines(cout);
	
	intervalBounds=testBHMfit.getBounds();
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

		vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.binHierarchy(samplingSteps);
		testBHMfit = matchedSplineFit(currentAnalysisBins, intervalBounds, 4, 0, dummy, dummy);
			
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
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.05; double currentVar2; pair<double,double> coarseGrainedResult; pair<double,double> basisResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2; currentVar2=currentVar-printStep/2;
		coarseGrainedResult=scaledBinHistogram.sampledFunctionValueWeightedAverage(currentVar);
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar2);
	
		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		<< '\t' << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar)
		<< '\t' << averageBootstrapSpline.splineValue(currentVar) << '\t' << averageBootstrapSpline.splineError(currentVar)
		<< '\t' << currentVar2 << '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}