#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double testFunction(double variable)
	{
	//return 3*(1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.)/10.;		//max 0.6168917686383567
	//return (2+variable)/8.;
	return (10+cos(variable*10.))/10./PI;			//max 1.1/PI, maxVar=PI+0.6
	//return exp(-3*variable)/(-1 + exp(6))*3*exp(9);	//max 3.1
	//return pow(variable,4)-0.8*variable*variable; 	//max 0.2
	}

int main(int argc, char **argv) {
	
	cout << "----------------- Production test: full default routine ----------------" << endl;
	
	long unsigned int seed=time(NULL);//11145792458;//956475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e5;
	
	double variable, random;
	double testFunctionMax=1.1/PI;
	double intervalSize=PI;
	double minVar=1.; double maxVar=PI+0.6;
	double slotWidth=(maxVar-minVar)/pow(2.,10); int numberOverlaps=1; int totalNumOfBasisFn=0;
	double printStep=0.02; int numSteps=int((maxVar-minVar)/printStep);
	pair<double,double> basisResult;
	slotBounds bounds1(1,2.8);
	//slotBounds bounds1(1,1.45);
	//slotBounds bounds1(-1,-0.5);
	//slotBounds bounds2(1.45,1.9);
	//slotBounds bounds2(-0.5,0);
	//slotBounds bounds3(1.9,2.8);
	//slotBounds bounds3(0, 0.5);
	//slotBounds bounds4(0.5, 1);
	int bf=4;
	double currentVar;
	vector< slotBounds > intervalBounds; 
	int numberOfGoodFits=0;
	
	vector<histogramBasis> binHistogramVector; vector<histogramBasis> basisHistogramVector;
	
	pair<double,double> total[numSteps][2]; for(int i=0; i<numSteps;i++) for(int j=0;j<2;j++) total[i][j]=make_pair(0,0);
	
	ofstream output("histogram_testoutput.dat");
	
	int maxRounds=100; int bootstrapSamples=100;
	for(int round=0;round<maxRounds+1;round++)
		{
		vector<basisSlot*> basisVector;
		basisVector.push_back(new taylorSlot(bounds1,bf)); 
		//basisVector.push_back(new taylorSlot(bounds2,bf)); 
		//basisVector.push_back(new taylorSlot(bounds3,bf)); 
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
		binHistogramVector[i%maxRounds].sample(variable,whatsign(testFunction(variable)));
		basisHistogramVector[i%maxRounds].sample(variable,whatsign(testFunction(variable)));
		binHistogramVector[maxRounds].sample(variable,whatsign(testFunction(variable)));
		basisHistogramVector[maxRounds].sample(variable,whatsign(testFunction(variable)));
		}

	for(int round=bootstrapSamples;round>=0;round--)
		{
		vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
		histogramBasis combinedHistogram(histogramVector);
		vector<basisSlot*> basisVector;
		basisVector.push_back(new taylorSlot(bounds1,bf)); 
		//basisVector.push_back(new taylorSlot(bounds2,bf)); 
		//basisVector.push_back(new taylorSlot(bounds3,bf)); 
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
		
		double threshold=2;
		bool acceptance;
		vector< double > dummy;
		splineArray testNewRoutine2;
		if(round==bootstrapSamples)
			{
			testNewRoutine2 = combinedHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
			acceptance=testNewRoutine2.getAcceptance();
			//if(testNewRoutine2.numberKnots()>10) acceptance=false;

			while(acceptance==false)
				{
				threshold+=1;
				testNewRoutine2 = combinedHistogram.splineProcedure(4, 2, samplingSteps, threshold, 0);
				acceptance=testNewRoutine2.getAcceptance();
				//if(testNewRoutine2.numberKnots()>10) acceptance=false;
				if(threshold>6) break;
				}
			
			intervalBounds=testNewRoutine2.getBounds();
			}
		else
			{
			vector< vector< basisSlot* > > currentAnalysisBins=combinedHistogram.binHierarchy(samplingSteps);
			testNewRoutine2 = matchedSplineFit(currentAnalysisBins, intervalBounds, 4, 0, dummy, dummy);
			acceptance=true;//testNewRoutine2.checkOverallAcceptance(threshold);
			}
				
		if( acceptance && round<bootstrapSamples) numberOfGoodFits++;

		for(int i=0; i<numSteps;i++) 
			{
			currentVar=minVar+i*printStep+printStep/2;
			basisResult=scaledBasisHistogram.sampledFunctionValueWeightedAverage(currentVar);
		
			if( acceptance && round<bootstrapSamples)
				{
				total[i][0].first+=testNewRoutine2.splineValue(currentVar); total[i][0].second+=testNewRoutine2.splineValue(currentVar)*testNewRoutine2.splineValue(currentVar);
				total[i][1].first+=basisResult.first; total[i][1].second+=basisResult.first*basisResult.first;
				}

			output << currentVar << '\t' << testNewRoutine2.splineValue(currentVar) << '\t' << testNewRoutine2.splineError(currentVar) << '\t'
			<< basisResult.first << '\t' << basisResult.second << endl;
			}
		output << endl; output << endl;
		
		cout << "fit info: " << testNewRoutine2.getAcceptance() << '\t' << threshold << '\t' << testNewRoutine2.numberKnots() << endl;
		}
		
	for(int i=0; i<numSteps;i++) 
		{
		output << minVar+i*printStep+printStep/2;
		for(int j=0; j<2;j++) 
			{
			total[i][j].first*=1./double(numberOfGoodFits); total[i][j].second*=1./double(numberOfGoodFits); total[i][j].second-=total[i][j].first*total[i][j].first;
			if(total[i][j].second<0) total[i][j].second=0;
			total[i][j].second=sqrt(total[i][j].second);
			output << '\t' << total[i][j].first << '\t' << total[i][j].second;
			}
		output << endl;
		}
		
	cout << "number of good fits: " << numberOfGoodFits << endl;
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}