#include "basic.hpp"
#include "spline.hpp"
#include "histogram.hpp"
#include "sput.h"

using namespace std; 

static void test_spline()
{
	double lowerBound=1.; double upperBound=3.;
	slotBounds myTestBounds(lowerBound, upperBound);
	unsigned int defaultSplineOrder=4;
	vector<double> mySplineCoeffs;
	for(unsigned int i=0;i<defaultSplineOrder;i++) mySplineCoeffs.push_back(i*0.5+1);
	splinePiece myTestSpline(myTestBounds,defaultSplineOrder);
	vector<double> covarianceDummy; covarianceDummy.resize(defaultSplineOrder*defaultSplineOrder);
	myTestSpline.setSplinePiece(mySplineCoeffs,covarianceDummy);
	
	sput_fail_unless(myTestSpline.splineValue(1.)==7, "correct spline value at 1");
	sput_fail_unless(myTestSpline.splineValue(1.5)==16.1875, "correct spline value at 1.5");
	
	lowerBound=1.5; upperBound=2.;
	slotBounds smallerBounds(lowerBound, upperBound);
	sput_fail_unless(myTestSpline.splineIntegral(smallerBounds)==4505./384., "correct spline integral between 1.5 and 2");
	
	mySplineCoeffs[0]=20; mySplineCoeffs[3]=-10;
	myTestSpline.setSplinePiece(mySplineCoeffs,covarianceDummy);
	sput_fail_unless(myTestSpline.splineValue(2.5)==-120., "correct updated spline value at 2.5");
	
}

static void test_matrix()
{
	int binMatSize=6;
	double* binomMatrix=binomialMatrix(binMatSize);

	sput_fail_unless(isAround(binomMatrix[0],1), "binomial matrix element 0,0");
	sput_fail_unless(isAround(binomMatrix[5],1), "binomial matrix element 0,5");
	sput_fail_unless(isAround(binomMatrix[6],0), "binomial matrix element 1,0");
	sput_fail_unless(isAround(binomMatrix[9],3), "binomial matrix element 1,3");
	sput_fail_unless(isAround(binomMatrix[16],6), "binomial matrix element 2,4");
	
	double* binomVec=binomialVector(binMatSize-1);
	
	sput_fail_unless(binomVec[0]==1&&binomVec[1]==4&&binomVec[2]==6&&binomVec[3]==4&&binomVec[4]==1, "binomial vector correct");

} 

static void test_chisq()
{
	double lowerBound=1.; double upperBound=4.2;
	slotBounds myTestBounds(lowerBound, upperBound);
	unsigned int defaultSplineOrder=4;
	
	double slotWidth=0.1; int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(lowerBound, upperBound, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram(histogramVector);
	
	double stepWidth=0.001; double variable=lowerBound+stepWidth/2.;
	while(variable<upperBound)
		{
		myTestHistogram.sample(variable,1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.);
		variable+=stepWidth;
		}

	vector<basisSlot*> allAnalysisSlots;
	allAnalysisSlots.push_back(myTestHistogram.combinedSlot(0,19));
	for(unsigned int i=0;i<2;i++) {allAnalysisSlots.push_back(myTestHistogram.combinedSlot(10*i,10*i+9));}
	for(unsigned int i=0;i<4;i++) {allAnalysisSlots.push_back(myTestHistogram.combinedSlot(5*i,5*i+4));}
	
	sput_fail_unless(isAround(allAnalysisSlots[1] -> sampledIntegral(),24666667./16000000.), "check value of one of the combined slots");
	sput_fail_unless(isAround( (allAnalysisSlots[4] -> sampledIntegralError()) ,sqrt(8488771217029311.)/16000000000.), "check error of one of the combined slots");

	double* myDesignMatrix=designMatrix(allAnalysisSlots, defaultSplineOrder);
	sput_fail_unless(isAround(myDesignMatrix[0],32000000000./(3.*sqrt(1535863353193979))), "check design matrix element");
	sput_fail_unless(isAround(myDesignMatrix[23],30750000000*sqrt(3./157587733676437.)), "check design matrix element");

	double* matrixProduct=designMatrixProduct(myDesignMatrix, allAnalysisSlots.size(), defaultSplineOrder);
	sput_fail_unless(isAround(matrixProduct[7],5675789.13893), "check A^T A element");
	sput_fail_unless(isAround(matrixProduct[13],5675789.13893), "check A^T A element");

	double* bvec = integralVector(allAnalysisSlots);
	double* vectorProduct=integralVectorProduct(myDesignMatrix, bvec, allAnalysisSlots.size(), defaultSplineOrder);
	sput_fail_unless(isAround(vectorProduct[2],3360322.181223839), "check A^T b element");
		
	//gsl example
	double a_data[] = { 0.18, 0.60, 0.57, 0.96,
		0.41, 0.24, 0.99, 0.58,
		0.14, 0.30, 0.97, 0.66,
		0.51, 0.13, 0.19, 0.85 };
	double b_data[] = { 1.0, 2.0, 3.0, 4.0 };
	solveLinearEquation(a_data, b_data, 4);
	sput_fail_unless(isAround(b_data[0],-4.05205,0.00001)&&isAround(b_data[3],8.69377,0.00001), "linear equation solver");
	
	double a_tonullify[]= {6.881737852316666e9, 1.804944376798533e8, -3.987145745998773e7, 6.611225244750259e8, 2.812625188625582e8,
		1.863234860715052e8, 4.947682846340233e6, -1.0813733676412262e6, 1.7904516220312424e7, 7.610831122786341e6,
		-9.293121263062338e8, -2.4406545147881947e7, 5.385261254688273e6, -8.928065050991334e7, -3.797944454921767e7,
		-9.670404462932392e8, -2.5378126931406867e7, 5.603280940338399e6, -9.290382520446008e7, -3.952272533408071e7,
		-4.835202230948418e8, -1.2689063397262562e7, 2.801640521947063e6, -4.64519132593478e7, -1.976136242589578e7,
		6.626954978308027e9, 1.7382944231461877e8, -3.839582860663368e7, 6.366470721982002e8, 2.708480698863848e8,
		-6.066708590964211e8, -1.584794310773435e7, 3.5129784269692334e6, -5.827752661106431e7, -2.479974675252765e7,
		-1.81680924938357e9, -4.7630685210806735e7, 1.0525620669933341e7, -1.7453771627848065e8, -7.42560330606024e7,
		1.2133374481887615e9, 3.1836434592198167e7, -7.030226858019163e6, 1.1656524571006072e8, 4.9589207537878215e7,
		6.066687241461586e8, 1.5918217364539955e7, -3.515113377231718e6, 5.82826221979126e7, 2.4794604010083687e7};
	double b_tonullify[]={12, 14, -5, 10, 2, 1, 0, 12, 6, 1};
		
	gsl_matrix_view gslU = gsl_matrix_view_array (a_tonullify, 10, 5);
	gsl_vector_view gslb = gsl_vector_view_array (b_tonullify, 10);
	gsl_matrix * V = gsl_matrix_alloc (5, 5);
	gsl_vector * diagonal = gsl_vector_alloc (5);
	gsl_vector *x = gsl_vector_alloc (5);
	solveSVD(&gslU.matrix, V, &gslb.vector, diagonal, x);
	
	sput_fail_unless(isAround(gsl_vector_get(diagonal,0),1e10), "SVD singular value matches");
	sput_fail_unless(isAround(gsl_vector_get(diagonal,2),25.), "SVD singular value matches");
	sput_fail_unless(isAround(gsl_vector_get(diagonal,4),0), "smallest singular value was nullified");
	sput_fail_unless(isAround(gsl_vector_get(x,0),-0.24185932309093,1e-6), "SVD solution vector matches");
	sput_fail_unless(isAround(gsl_vector_get(x,4),-1.11168398015165,1e-6), "SVD solution vector matches");
	sput_fail_unless(isAround(gsl_matrix_get(V,0,0),-0.99423953219004,1e-7), "SVD V matrix matches");
	sput_fail_unless(isAround(gsl_matrix_get(V,1,3),0.09675019612512,1e-7), "SVD V matrix matches");
	sput_fail_unless(isAround(gsl_matrix_get(V,3,2),-0.09551581016585,1e-7), "SVD V matrix matches");
	sput_fail_unless(isAround(gsl_matrix_get(V,4,4),-0.93581193011763,1e-7), "SVD V matrix matches");
	
}

int main(int argc, char **argv) {
	
	sput_start_testing();
	
	sput_enter_suite("test_spline()");
	sput_run_test(test_spline);
	
	sput_enter_suite("test_matrix()");
	sput_run_test(test_matrix);
	
	sput_enter_suite("test_chisq()");
	sput_run_test(test_chisq);
	
	sput_finish_testing();
	
	return sput_get_return_value();
	
}

