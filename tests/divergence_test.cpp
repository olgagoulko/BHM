#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	return tan(variable*PI/2.)*tan(variable*PI/2.);
	}
	
double variableTransform(double variable)
	{
	return 2*atan(variable)/PI;
	}
	
double transformDerivative(double variableTransformed)
	{
	return PI/2./cos(variableTransformed*PI/2.)/cos(variableTransformed*PI/2.);
	}
	
double divergenceWeight(double variable)
	{
	return sqrt(variable);
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Divergent function on a semi-infinite domain ----------------" << endl;
	
	long unsigned int seed=time(NULL);//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	double variable; double variableTransformed;
	
	//make basis with square root
	vector<basisSlot*> basisVector; int bf=4;
	slotBounds bounds1(0, 0.5); basisVector.push_back(new sqrtSlot(bounds1,bf));
	slotBounds bounds2(0.5, 1.0); basisVector.push_back(new taylorSlot(bounds2,bf));

	//initialize histogram and sample function
	//double minVar=0; double maxVar=10.0;
	double minVar=0; double maxVar=1;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	histogramBasis basisHistogram(basisVector);
	
	int maxRounds=100; int bootstrapSamples=100; vector<histogramBasis> binHistogramVector;
	vector< vector<double> > averageCoeffs; vector< vector<double> > averageCov;
	for(int round=0;round<maxRounds;round++)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis binHistogram(histogramVector);
		binHistogramVector.push_back(binHistogram);
		}

	for(int i=0; i<samplingSteps;i++)
		{
		variable=gsl_rng_uniform(RNG);
		variable=testFunction(variable);
	
		//plain sampling
		//basisHistogram.sample(variable,1);
		//binHistogram.sampleUniform(variable,1);
	
		//just compensate for divergence
		//basisHistogram.sample(variable,divergenceWeight(variable));
		//binHistogram.sampleUniform(variable,divergenceWeight(variable));

		//interval transform
		variableTransformed=variableTransform(variable);
		
		//without compensation for divergence
		//basisHistogram.sample(variableTransformed,1/transformDerivative(variableTransformed));
		//binHistogram.sampleUniform(variableTransformed,1/transformDerivative(variableTransformed));
		
		//with compensation for divergence
		basisHistogram.sample(variableTransformed,divergenceWeight(variable)/transformDerivative(variableTransformed));
		binHistogram.sampleUniform(variableTransformed,divergenceWeight(variable)/transformDerivative(variableTransformed));
		binHistogramVector[i%maxRounds].sampleUniform(variableTransformed,divergenceWeight(variable)/transformDerivative(variableTransformed));
		}

	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	double threshold=2;
	splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0);
	bool acceptance=testBHMfit.getAcceptance();
	while(acceptance==false)
		{
		threshold+=1;
		testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, threshold, 0);
		acceptance=testBHMfit.getAcceptance();
		if(threshold>5) break;
		}
	testBHMfit.printSplines();
	
	vector<slotBounds> intervalBounds=testBHMfit.getBounds();
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
		splineArray testfit2;
		
		for(unsigned int j=0;j<intervalBounds.size();j++)
			{
			vector<double> theCoeffVec(4, 0.0); vector<double> theCovVec(10, 0.0);
			averageCoeffs.push_back(theCoeffVec);
			averageCov.push_back(theCovVec);
			}
		
		vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.binHierarchy(samplingSteps);
		testfit2 = matchedSplineFit(currentAnalysisBins, intervalBounds, 4, 0, dummy, dummy);
		
		for(unsigned int j=0;j<intervalBounds.size();j++)
			{
			vector<double> currentCoeffs=testfit2.getSplinePiece(j) -> getCoefficients();
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

	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; double currentVar; pair<double,double> basisResult;
	
	//print without any transform
	/*for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
		
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
		output << currentVar << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(currentVar).first << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(currentVar).second << '\t' 
		<< '\t' << testNewRoutine.splineValue(currentVar) << '\t' << testNewRoutine.splineError(currentVar)
		//<< '\t' << testcurvature.splineValue(currentVar) << '\t' << testcurvature.splineError(currentVar)
		<< '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}*/
	
	//variable and divergence transform
	/**/maxVar=10;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
		variableTransformed=variableTransform(currentVar);
		
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(variableTransformed);
		output << currentVar << '\t' << scaledBinHistogram.sampledFunctionValueWeightedAverage(variableTransformed).first/divergenceWeight(currentVar) << '\t' 
		<< scaledBinHistogram.sampledFunctionValueWeightedAverage(variableTransformed).second/divergenceWeight(currentVar) << '\t' 
		<< '\t' << testBHMfit.splineValue(variableTransformed)/divergenceWeight(currentVar) << '\t' << testBHMfit.splineError(variableTransformed)/divergenceWeight(currentVar)
		<< '\t' << basisResult.first/divergenceWeight(currentVar) << '\t' << basisResult.second/divergenceWeight(currentVar)
		<< '\t' << averageBootstrapSpline.splineValue(variableTransformed)/divergenceWeight(currentVar) << '\t' << averageBootstrapSpline.splineError(variableTransformed)/divergenceWeight(currentVar)
		<< endl;
		}/**/
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}