#include "histogram.hpp"
#include "matrix.hpp"
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
	
	cout << "----------------- Production test with sign problem and Markov chain updates ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	long MCsteps=1e6;
	
	vector<basisSlot*> basisVector; int bf=4;
	vector<slotBounds> intervalBounds;
	slotBounds bounds1(0,3); taylorSlot slot1(bounds1,bf); basisVector.push_back(&slot1); intervalBounds.push_back(bounds1);
	histogramBasis basisHistogram(basisVector);

	double minVar=0.; double maxVar=3.0; double maxVarTotal=3.0;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis binHistogram(histogramVector);

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
			binHistogram.sample(currentVar,currentSector);
			basisHistogram.sample(currentVar,currentSector);
			}
		//if(currentSector==1) sectorCounter++;
		}
		
	//cout << "Sector counter = " << sectorCounter/double(samplingSteps) << endl;
		
	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);
	histogramBasis scaledBasisHistogram = basisHistogram.scaledHistogram(samplingSteps);

	splineArray testNewRoutine = binHistogram.splineProcedure(4, 2, samplingSteps, 2);
	cout << "New spline procedure more big bins" << endl;
	testNewRoutine.printSplineArrayInfo(); cout << endl;
	testNewRoutine.printSplines();
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.05; double currentVar2; pair<double,double> coarseGrainedResult; pair<double,double> basisResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2; currentVar2=currentVar-printStep/2;
		coarseGrainedResult=scaledBinHistogram.sampledFunctionValueWeightedAverage(currentVar);
		basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar2);
	
		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
		//<< '\t' << bestFit.splineValue(currentVar) << '\t' << bestFit.splineError(currentVar)
		//<< '\t' << testNewRoutine.splineValue(currentVar) << '\t' << testNewRoutine.splineError(currentVar)
		<< '\t' << testNewRoutine.splineValue(currentVar) << '\t' << testNewRoutine.splineError(currentVar)
		<< '\t' << currentVar2 << '\t' << basisResult.first << '\t' << basisResult.second
		<< endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}