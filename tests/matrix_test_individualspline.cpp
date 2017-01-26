#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double polynomialTestFunction(double variable)
	{
	//return 1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.;
	//return 2+variable;
	return 2*cos(variable*10.);
	}

int main(int argc, char **argv) {
	
	cout << "---------- Production test: sampling realistic function and fitting one spline -----------" << endl;
	
	long unsigned int seed=456475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e6;
	double variable;
	
	//sampling 3rd degree polynomial, which can be exactly fitted
	double minVar=1.; double maxVar=3.; double slotWidth=0.01; int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram(histogramVector);
	
	for(int i=0; i<samplingSteps;i++)
		{
		variable=minVar+(maxVar-minVar)*gsl_rng_uniform(RNG);
		myTestHistogram.sample(variable,polynomialTestFunction(variable));
		}
		
	//basisSlot * checkSlot = myTestHistogram.getSlot(0);
	//cout << checkSlot -> sampledFunctionValue(1.005) << '\t' << checkSlot -> sampledFunctionError(1.005) << endl;
		
	int defaultSplineOrder=4;
	//very good spline, more accurate than just raw data
	spline fineSpline=myTestHistogram.oneSplineFit(defaultSplineOrder);
	fineSpline.printSpline(); cout << endl;
	
	//worse accuracy in cubic due to overfitting, linear data ok
	spline higherOrderSpline=myTestHistogram.oneSplineFit(defaultSplineOrder+2);
	higherOrderSpline.printSpline(); cout << endl;
	
	//worse accuracy in cubic due to underfitting, linear data ok
	spline lowerOrderSpline=myTestHistogram.oneSplineFit(defaultSplineOrder-2);
	lowerOrderSpline.printSpline(); cout << endl;
	
	vector<basisSlot*> analysisSlots; unsigned int totalNumberSlots=histogramVector.size();
	analysisSlots.push_back(myTestHistogram.combinedSlot(0,totalNumberSlots-1));
	for(unsigned int i=0;i<2;i++) {analysisSlots.push_back(myTestHistogram.combinedSlot(totalNumberSlots*i/2,totalNumberSlots*(i+1)/2-1));}
	for(unsigned int i=0;i<4;i++) {analysisSlots.push_back(myTestHistogram.combinedSlot(totalNumberSlots*i/4,totalNumberSlots*(i+1)/4-1));}
	
	//rough spline less accurate, does not catch very bad fit of unsuitable function via chi^2
	histogramBasis analysisHistogram(analysisSlots);
	spline roughSpline=analysisHistogram.oneSplineFit(defaultSplineOrder);
	roughSpline.printSpline(); cout << endl;
	
	for(unsigned int i=0;i<8;i++) {analysisSlots.push_back(myTestHistogram.combinedSlot(totalNumberSlots*i/8,totalNumberSlots*(i+1)/8-1));}
	
	//still too rough, but better than before because more slots
	histogramBasis analysisHistogram2(analysisSlots);
	spline roughSpline2=analysisHistogram2.oneSplineFit(defaultSplineOrder);
	roughSpline2.printSpline(); cout << endl;
	
	//adding large combined slots to fine grid doesn't change outcome at all
	histogramBasis myTestHistogram2=myTestHistogram;
	for(unsigned int i=0;i<analysisSlots.size();i++) {myTestHistogram2.appendSlot(analysisSlots[i]);}
	spline allSpline=myTestHistogram2.oneSplineFit(defaultSplineOrder);
	allSpline.printSpline(); cout << endl;
	
	//doubling the slot width doubles the error bars and errors
	vector<basisSlot*> doubleSlots;
	for(unsigned int i=0;i<totalNumberSlots/2;i++) {doubleSlots.push_back(myTestHistogram.combinedSlot(2*i,2*i+1));}
	histogramBasis doubleHistogram(doubleSlots);
	spline doubleSpline=doubleHistogram.oneSplineFit(defaultSplineOrder);
	doubleSpline.printSpline(); cout << endl;
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; double currentVar; pair<double,double> sampledResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
		sampledResult=myTestHistogram.sampledFunctionValueWeightedAverage(currentVar);
		output << currentVar << '\t' << sampledResult.first << '\t' << sampledResult.second << '\t' << fineSpline.splineValue(currentVar) << '\t' << fineSpline.splineError(currentVar) << '\t' << roughSpline.splineValue(currentVar) << '\t' << roughSpline.splineError(currentVar) << '\t' << roughSpline2.splineValue(currentVar) << '\t' << roughSpline2.splineError(currentVar) << '\t' << allSpline.splineValue(currentVar) << '\t' << allSpline.splineError(currentVar) << '\t' << doubleSpline.splineValue(currentVar) << '\t' << doubleSpline.splineError(currentVar) << endl;
		}


	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}