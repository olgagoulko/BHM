#include "histogram.hpp"
#include "matrix.hpp"
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
	
	cout << "----------------- Production test: full default routine ----------------" << endl;
	
	long unsigned int seed=time(NULL);//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	double variable; double variableTransformed;
	
	//make basis with square roots
	vector<basisSlot*> basisVector; int bf=4;
	vector<slotBounds> intervalBounds;
	
	/*slotBounds bounds1(0, 1.25); sqrtSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	slotBounds bounds2(1.25, 2.5); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);
	slotBounds bounds3(2.5, 5); taylorSlot slot3(bounds3,bf); basisVector.push_back(&slot3); intervalBounds.push_back(bounds3);
	slotBounds bounds4(5, 7.5); taylorSlot slot4(bounds4,bf); basisVector.push_back(&slot4); intervalBounds.push_back(bounds4);
	slotBounds bounds5(7.5, 10); taylorSlot slot5(bounds5,bf); basisVector.push_back(&slot5); intervalBounds.push_back(bounds5);*/
	
	slotBounds bounds1(0, 0.5); sqrtSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	slotBounds bounds2(0.5, 1.0); taylorSlot slot2(bounds2,bf); basisVector.push_back(&slot2); intervalBounds.push_back(bounds2);

	//initialize histogram and sample function
	//double minVar=0; double maxVar=10.0;
	double minVar=0; double maxVar=1;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);
	histogramBasis basisHistogram(basisVector);

	for(int i=0; i<samplingSteps;i++)
		{
		variable=gsl_rng_uniform(RNG);
		variable=testFunction(variable);
	
		//plain sampling
		//basisHistogram.sample(variable,1);
		//binHistogram.sample(variable,1);
	
		//just compensate for divergence
		//basisHistogram.sample(variable,divergenceWeight(variable));
		//binHistogram.sample(variable,divergenceWeight(variable));

		//interval transform
		variableTransformed=variableTransform(variable);
		
		//without compensation for divergence
		//basisHistogram.sample(variableTransformed,1/transformDerivative(variableTransformed));
		//binHistogram.sample(variableTransformed,1/transformDerivative(variableTransformed));
		
		//with compensation for divergence
		basisHistogram.sample(variableTransformed,divergenceWeight(variable)/transformDerivative(variableTransformed));
		binHistogram.sample(variableTransformed,divergenceWeight(variable)/transformDerivative(variableTransformed));
		}

	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	double threshold=2;
	splineArray testNewRoutine = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
	bool acceptance=testNewRoutine.getAcceptance();
	while(acceptance==false)
		{
		threshold+=1;
		testNewRoutine = binHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
		acceptance=testNewRoutine.getAcceptance();
		if(threshold>5) break;
		}
	testNewRoutine.printSplines();
	
	/*splineArray testcurvature = binHistogram.splineProcedure(4, 2, samplingSteps, 2, 1);
	cout << "Suppress curvature" << endl;
	testcurvature.printSplines();*/

	
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
		<< '\t' << testNewRoutine.splineValue(variableTransformed)/divergenceWeight(currentVar) << '\t' << testNewRoutine.splineError(variableTransformed)/divergenceWeight(currentVar)
		<< '\t' << basisResult.first/divergenceWeight(currentVar) << '\t' << basisResult.second/divergenceWeight(currentVar)
		<< endl;
		}/**/
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}