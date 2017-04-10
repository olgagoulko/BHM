#include "histogram.hpp"
#include "matrix.hpp"
#include "basic.hpp"
#include "slot.hpp"

using namespace std;  

double polynomialTestFunction(double variable)
	{
	//return 1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.;
	//return 2+variable;
	//return 2*cos(variable*10.);
	return 5*exp(-3*variable);
	}

int main(int argc, char **argv) {
	
	cout << "---------- Production test: sampling realistic function and fitting several splines -----------" << endl;
	
	long unsigned int seed=456475;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e6;
	double variable;
	
	//sampling 3rd degree polynomial, which can be exactly fitted
	double minVar=1.; double maxVar=3.; double slotWidth=0.01; int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram(histogramVector);
	unsigned int totalNumberSlots=histogramVector.size();
	cout << "Total number of bins is " << totalNumberSlots << endl; cout << endl;
	
	for(int i=0; i<samplingSteps;i++)
		{
		variable=minVar+(maxVar-minVar)*gsl_rng_uniform(RNG);
		myTestHistogram.sample(variable,polynomialTestFunction(variable));
		}
		
	int defaultSplineOrder=4;
	
	spline allSpline=myTestHistogram.oneSplineFit(defaultSplineOrder);
	//cout << "Spline over all data" << endl; allSpline.printSpline(); cout << endl;
		
	vector<basisSlot*> slots1; vector<basisSlot*> slots2; vector<basisSlot*> slots3; vector<basisSlot*> slots4; 
	for(unsigned int i=0;i<totalNumberSlots/4;i++) {slots1.push_back(myTestHistogram.getSlot(i));}
	for(unsigned int i=totalNumberSlots/4;i<totalNumberSlots/2;i++) {slots2.push_back(myTestHistogram.getSlot(i));}
	for(unsigned int i=totalNumberSlots/2;i<3*totalNumberSlots/4;i++) {slots3.push_back(myTestHistogram.getSlot(i));}
	for(unsigned int i=3*totalNumberSlots/4;i<totalNumberSlots;i++) {slots4.push_back(myTestHistogram.getSlot(i));}
	histogramBasis analysisHistogram1(slots1);
	histogramBasis analysisHistogram2(slots2);
	histogramBasis analysisHistogram3(slots3);
	histogramBasis analysisHistogram4(slots4);

	spline spline1=analysisHistogram1.oneSplineFit(defaultSplineOrder);
	//cout << "Spline over interval 0" << endl; spline1.printSpline(); cout << endl;
	spline spline2=analysisHistogram2.oneSplineFit(defaultSplineOrder);
	//cout << "Spline over interval 1" << endl; spline2.printSpline(); cout << endl;
	spline spline3=analysisHistogram3.oneSplineFit(defaultSplineOrder);
	//cout << "Spline over interval 2" << endl; spline3.printSpline(); cout << endl;
	spline spline4=analysisHistogram4.oneSplineFit(defaultSplineOrder);
	//cout << "Spline over interval 3" << endl; spline4.printSpline(); cout << endl;
	
	vector<unsigned int> intervalBoundaries; intervalBoundaries.push_back(totalNumberSlots/4); intervalBoundaries.push_back(totalNumberSlots/2); intervalBoundaries.push_back(3*totalNumberSlots/4);
	splineArray theFit = myTestHistogram.splineFit(intervalBoundaries, 0, defaultSplineOrder);
	//cout << "Spline over intervals" << endl; theFit.printSplineArrayInfo(); theFit.printSplines();
	
	splineArray theMatchedFit = myTestHistogram.splineFit(intervalBoundaries, 1e10, defaultSplineOrder);
	cout << "Smoothed spline over intervals" << endl; theMatchedFit.printSplineArrayInfo(); theMatchedFit.printSplines();
	
	splineArray exactlyMatchedFit = myTestHistogram.matchedSplineFitSmallBins(intervalBoundaries, defaultSplineOrder);
	cout << "Exactly matched spline over intervals" << endl; exactlyMatchedFit.printSplineArrayInfo(); exactlyMatchedFit.printSplines();
	
	ofstream output("histogram_testoutput.dat");
	double printStep=0.01; double currentVar; pair<double,double> sampledResult;
	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
		{
		currentVar=minVar+i*printStep+printStep/2;
		sampledResult=myTestHistogram.sampledFunctionValueWeightedAverage(currentVar);
		output << currentVar << '\t' << sampledResult.first << '\t' << sampledResult.second << '\t' << allSpline.splineValue(currentVar) << '\t' << allSpline.splineError(currentVar) << '\t' << theFit.splineValue(currentVar) << '\t' << theFit.splineError(currentVar) << '\t' << theMatchedFit.splineValue(currentVar) << '\t' << theMatchedFit.splineError(currentVar) << '\t' << exactlyMatchedFit.splineValue(currentVar) << '\t' << exactlyMatchedFit.splineError(currentVar) << '\t';
		if(currentVar<=1.5) output << spline1.splineValue(currentVar) << '\t' << spline1.splineError(currentVar) << endl;
		else if(currentVar<=2.0) output << spline2.splineValue(currentVar) << '\t' << spline2.splineError(currentVar) << endl;
		else if(currentVar<=2.5) output << spline3.splineValue(currentVar) << '\t' << spline3.splineError(currentVar) << endl;
		else output << spline4.splineValue(currentVar) << '\t' << spline4.splineError(currentVar) << endl;
		}


	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}