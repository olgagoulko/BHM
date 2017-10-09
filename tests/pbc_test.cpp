#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
{
	return pow(variable,4)-0.8*variable*variable; 	//max 0.2
}


int main(int argc, char **argv) {
	
	cout << "----------------- Test distribution with periodic boundary conditions ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e6;

	vector<slotBounds> intervalBounds;

	double testFunctionMax=0.2;
	double intervalSize=2;
	double minVar=-1.; double maxVar=1.0;
	double slotWidth=(maxVar-minVar)/pow(2.,8);
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
	
	histogramBasis binHistogram(histogramVector);
	
	int maxRounds=100; int bootstrapSamples=100; vector<histogramBasis> binHistogramVector;
	vector< vector<double> > averageCoeffs; vector< vector<double> > averageCov;
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
		histogramBasis binHistogram(histogramVector);
		binHistogramVector.push_back(binHistogram);
		}

	double variable, random;
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
		binHistogramVector[i%maxRounds].sampleUniform(variable,whatsign(testFunction(variable)));
		}

	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);

	fitAcceptanceThreshold threshold; threshold.min=2; threshold.max=2; threshold.steps=0;
	splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0, true, false);	//type of boundary condition to be entered here
	testBHMfit.printSplineArrayInfo(cout); cout << endl;
	testBHMfit.printSplines(cout);
	
	intervalBounds=testBHMfit.getBounds();
	for(int round=0;round<bootstrapSamples;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
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
	double printStep=0.01; pair<double,double> coarseGrainedResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		variable=minVar+i*printStep+printStep/2;
		coarseGrainedResult=scaledBinHistogram.sampledFunctionValueWeightedAverage(variable);
	
		output << variable << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		<< '\t' << testBHMfit.splineValue(variable) << '\t' << testBHMfit.splineError(variable)
		<< '\t' << averageBootstrapSpline.splineValue(variable) << '\t' << averageBootstrapSpline.splineError(variable)
		<< endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}