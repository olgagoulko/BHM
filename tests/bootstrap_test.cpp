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
	
	cout << "----------------- Test of the bootstrap error estimate ----------------" << endl;
	
	long unsigned int seed=956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e4;
	
	double testFunctionMax=3.1;
	double intervalSize=2;
	double minVar=1.; double maxVar=2.8;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	
	double printStep=0.01; int numSteps=int((maxVar-minVar)/printStep)+1;
	pair<double,double> basisResult;
	slotBounds bounds1(1,2.8);
	//slotBounds bounds1(1,1.45);
	//slotBounds bounds1(-1,-0.5);
	//slotBounds bounds1b(1.45,1.9);
	//slotBounds bounds2(-0.5,0);
	//slotBounds bounds2(1.9,2.8);
	//slotBounds bounds3(0, 0.5);
	//slotBounds bounds4(0.5, 1);
	int bf=4;
	double variable, random;
	vector< slotBounds > intervalBounds; 
	int numberOfGoodFits=0;

	vector<histogramBasis> binHistogramVector; vector<histogramBasis> basisHistogramVector;

	pair<double,double> total[numSteps][2]; for(int i=0; i<numSteps;i++) for(int j=0;j<2;j++) total[i][j]=make_pair(0,0);
	
	vector< vector<double> > averageCoeffs; vector< vector<double> > averageCov;
	
	ofstream output("histogram_testoutput.dat");

	int maxRounds=1000; int bootstrapSamples=1000;
	for(int round=0;round<maxRounds+1;round++)
		{
		vector<basisSlot*> basisVector;
		basisVector.push_back(new taylorSlot(bounds1,bf)); 
		//basisVector.push_back(new taylorSlot(bounds1b,bf)); 
		//basisVector.push_back(new taylorSlot(bounds2,bf)); 
		//basisVector.push_back(new taylorSlot(bounds4,bf)); 
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
			
		histogramBasis binHistogram(histogramVector);
		histogramBasis basisHistogram(basisVector);
		
		binHistogramVector.push_back(binHistogram); basisHistogramVector.push_back(basisHistogram);
		}

	for(int i=0; i<samplingSteps;i++)
		{
		bool accept=false;
		while(accept==false)
			{
			variable=gsl_rng_uniform(RNG)*intervalSize+minVar;
			random=gsl_rng_uniform(RNG)*testFunctionMax;
			if(random<abs(testFunction(variable))) accept=true;
			}
		binHistogramVector[i%maxRounds].sampleUniform(variable,whatsign(testFunction(variable)));
		basisHistogramVector[i%maxRounds].sample(variable,whatsign(testFunction(variable)));
		binHistogramVector[maxRounds].sampleUniform(variable,whatsign(testFunction(variable)));
		basisHistogramVector[maxRounds].sample(variable,whatsign(testFunction(variable)));
		}

	for(int round=bootstrapSamples;round>=0;round--)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis combinedHistogram(histogramVector);
		vector<basisSlot*> basisVector;
		basisVector.push_back(new taylorSlot(bounds1,bf)); 
		//basisVector.push_back(new taylorSlot(bounds1b,bf)); 
		//basisVector.push_back(new taylorSlot(bounds2,bf)); 
		//basisVector.push_back(new taylorSlot(bounds4,bf)); 
		histogramBasis combinedBasisHistogram(basisVector);
		if(round==bootstrapSamples)
			{
			combinedHistogram.addAnotherHistogram(binHistogramVector[maxRounds]);
			combinedBasisHistogram.addAnotherHistogram(basisHistogramVector[maxRounds]);
			}
		else 
			{
			for(int i=0;i<maxRounds;i++)
				{
				unsigned int random=gsl_rng_uniform_int (RNG, maxRounds);
				combinedHistogram.addAnotherHistogram(binHistogramVector[random]);
				combinedBasisHistogram.addAnotherHistogram(basisHistogramVector[random]);
				}
			}
		histogramBasis scaledBasisHistogram = combinedBasisHistogram.scaledHistogram(samplingSteps);
		
		fitAcceptanceThreshold threshold; threshold.min=2; threshold.max=6; threshold.steps=4;
		bool acceptance;
		vector< double > dummy;
		splineArray theBHMfit;
		
		if(round==bootstrapSamples)
			{
			theBHMfit = combinedHistogram.BHMfit(4, 2, samplingSteps, threshold, 0, false);
			acceptance=theBHMfit.getAcceptance();
			//if(theBHMfit.numberKnots()>10) acceptance=false;
			
			intervalBounds=theBHMfit.getBounds();
			
			for(unsigned int j=0;j<intervalBounds.size();j++)
				{
				vector<double> theCoeffVec(4, 0.0); vector<double> theCovVec(10, 0.0);
				averageCoeffs.push_back(theCoeffVec);
				averageCov.push_back(theCovVec);
				}
			}
		else
			{
			vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.binHierarchy(samplingSteps);
			theBHMfit = matchedSplineFit(currentAnalysisBins, intervalBounds, 4, 0, dummy, dummy);
			acceptance=true;//theBHMfit.checkOverallAcceptance(threshold);
			if( acceptance ) numberOfGoodFits++;
			
			for(unsigned int j=0;j<intervalBounds.size();j++)
				{
				vector<double> currentCoeffs=theBHMfit.getSplinePiece(j) -> getCoefficients();
				vector<double> delta;
				for(unsigned int i=0;i<4;i++)
					{
					delta.push_back(currentCoeffs[i] - averageCoeffs[j][i]);
					averageCoeffs[j][i] += delta[i] / numberOfGoodFits;
					}
				for(unsigned int i=0;i<4;i++)
					for(unsigned int k=i;k<4;k++)
						{
						averageCov[j][i*(7-i)/2+k] += delta[i]*delta[k]*(numberOfGoodFits-1)/numberOfGoodFits;
						}
				}
			}

		for(int i=0; i<numSteps;i++) 
			{
			variable=minVar+i*printStep;
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(variable);
		
			if( acceptance && round<bootstrapSamples)
				{
				total[i][0].first+=theBHMfit.splineValue(variable); total[i][0].second+=theBHMfit.splineValue(variable)*theBHMfit.splineValue(variable);
				total[i][1].first+=basisResult.first; total[i][1].second+=basisResult.first*basisResult.first;
				}

			output << variable << '\t' << theBHMfit.splineValue(variable) << '\t' << theBHMfit.splineError(variable) << '\t'
			<< basisResult.first << '\t' << basisResult.second << endl;
			}
		output << endl; output << endl;
		
		cout << "fit info: " << theBHMfit.getAcceptance() << '\t' << theBHMfit.getThreshold() << '\t' << theBHMfit.numberKnots() << endl;
		}
	
	vector<splinePiece*> theBootstrapSplines;
	for(unsigned int i=0;i<intervalBounds.size();i++) 
		{
		theBootstrapSplines.push_back(new splinePiece(intervalBounds[i]));
	
		vector<double> theErrorCoefficients; for(unsigned int j=0;j < 7;j++) theErrorCoefficients.push_back(0); 
		for(unsigned int j=0;j < 4;j++)
			for(unsigned int k=j;k < 4;k++)
				{
				theErrorCoefficients[j+k]+=averageCov[i][j*(7-j)/2+k]/(numberOfGoodFits - 1);
				if(j!=k) theErrorCoefficients[j+k]+=averageCov[i][j*(7-j)/2+k]/(numberOfGoodFits - 1);
				}
		theBootstrapSplines[i] -> setSplinePiece(averageCoeffs[i], theErrorCoefficients);
		}

	splineArray averageBootstrapSpline(theBootstrapSplines);
	cout << "Average Bootstrap Spline" << endl;
	averageBootstrapSpline.printSplineArrayInfo(cout); averageBootstrapSpline.printSplines(cout); cout << endl;
		
	for(int i=0; i<numSteps;i++) 
		{
		output << minVar+i*printStep;
		for(int j=0; j<2;j++) 
			{
			total[i][j].first*=1./double(numberOfGoodFits); total[i][j].second*=1./double(numberOfGoodFits); total[i][j].second-=total[i][j].first*total[i][j].first;
			if(total[i][j].second<0) total[i][j].second=0;
			total[i][j].second=sqrt(total[i][j].second);
			output << '\t' << total[i][j].first << '\t' << total[i][j].second;
			}
		output << endl;
		}
	output << endl; output << endl;
	for(int i=0; i<numSteps;i++) 
		{
		variable=minVar+i*printStep;
		output << variable << '\t' << averageBootstrapSpline.splineValue(variable) << '\t' << averageBootstrapSpline.splineError(variable) << endl;
		}
		
	cout << "number of good fits: " << numberOfGoodFits << endl;
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}