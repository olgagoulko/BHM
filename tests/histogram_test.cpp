#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

int main(int argc, char **argv) {
	
	cout << "---------- Production test: sampling realistic function with taylorSlot -----------" << endl;
	
	long unsigned int seed=456475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e7;
	double variable;
	
	double minVar=0.5; double maxVar=6.5; double slotWidth=0.5; int numberOverlaps=10; int totalNumOfBasisFn=4;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram(histogramVector);
	
	for(int i=0; i<samplingSteps;i++)
		{
		variable=minVar+(maxVar-minVar)*gsl_rng_uniform(RNG);
		myTestHistogram.sample(variable,1./sqrt(variable));
		}
		
	ofstream output("histogram_testoutput.dat");
	ofstream outputwav("histogram_testoutput_weightedav.dat");
	double printStep=0.01; double currentVar; pair<double,double> sampledResult; pair<double,double> sampledWeightedResult;
	unsigned int vectorSize;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
	{
		currentVar=minVar+i*printStep;
		sampledResult=myTestHistogram.sampledFunctionValueAverage(currentVar);
		sampledWeightedResult=myTestHistogram.sampledFunctionValueWeightedAverage(currentVar);
		output << currentVar << '\t' << sampledResult.first << '\t' << sampledResult.second << endl;
		outputwav << currentVar << '\t' << sampledWeightedResult.first << '\t' << sampledWeightedResult.second << endl;
		vectorSize=i;
	}
	
	vectorSize++;
	double scalingErrorMatrix=1e-9; double scaledElement;
	double* matrix=secondDerivativeMatrix(vectorSize);
	double vector[vectorSize];
	for(unsigned int i=0; i<vectorSize;i++) 
	{
		currentVar=minVar+i*printStep;
		sampledResult=myTestHistogram.sampledFunctionValueAverage(currentVar);
		scaledElement=scalingErrorMatrix/sampledResult.second/sampledResult.second;
		matrix[i*vectorSize+i]+=scaledElement;
		vector[i]=sampledResult.first*scaledElement;
	}
	
	solveLinearEquation(matrix,vector,vectorSize);
	
	ofstream output2("histogram_testoutput2.dat");
	for(unsigned int i=0; i<vectorSize;i++) 
	{
		currentVar=minVar+i*printStep;
		output2 << currentVar << '\t' << vector[i] << endl;
	}
	
	scalingErrorMatrix=1e-7;
	for(unsigned int i=0; i<vectorSize;i++) 
	{
		for(unsigned int j=0; j<vectorSize;j++) matrix[i*vectorSize+j]=0;
		currentVar=minVar+i*printStep;
		sampledResult=myTestHistogram.sampledFunctionValueAverage(currentVar);
		scaledElement=scalingErrorMatrix/sampledResult.second/sampledResult.second;
		matrix[i*vectorSize+i]=scaledElement+1.;
		vector[i]=(sampledResult.first*scaledElement)+(1./sqrt(currentVar));
	}
	
	solveLinearEquation(matrix,vector,vectorSize);
	
	ofstream output3("histogram_testoutput3.dat");
	for(unsigned int i=0; i<vectorSize;i++) 
	{
		currentVar=minVar+i*printStep;
		output3 << currentVar << '\t' << vector[i] << endl;
	}

	cout << "----------------------------- Testing finished ------------------------------------" << endl;
	
	return 0;
	
}