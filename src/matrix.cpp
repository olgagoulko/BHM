#include "basic.hpp"
#include "slot.hpp"
#include "matrix.hpp"

using namespace std; 


spline::spline(slotBounds theBounds, unsigned int theSplineOrder){
	
bounds=theBounds;
splineOrder=theSplineOrder;
splineCoefficients.resize(splineOrder);
splineErrorCoefficients.resize(2*splineOrder-1);
chiSquared=0;
degreesOfFreedom=1;
	
}

void spline::setSpline(vector< double > theCoefficients, vector< double > theErrorCoefficients, double theChiSquared, int theDOF)
{
	splineCoefficients=theCoefficients;
	splineErrorCoefficients=theErrorCoefficients;
	chiSquared=theChiSquared;
	degreesOfFreedom=theDOF;
}

vector< double > spline::getCoefficients() const
{
return splineCoefficients;
}

double spline::goodnessOfFitQ() const
{
	return gsl_sf_gamma_inc_Q (double(degreesOfFreedom)/2., chiSquared/2.);
}


double spline::splineValue(double variable) const{
	
//think about what to do when out of bounds: error probably not warranted; might return zero or actual spline value
double result=0;
double currentVariablePower=1.;
for(unsigned int i=0;i<splineOrder;i++)
	{
	result+=splineCoefficients[i]*currentVariablePower;
	currentVariablePower*=variable;
	}
return result;

}

double spline::splineError(double variable) const{
	
	//think about what to do when out of bounds: error probably not warranted; might return zero or actual spline value
	double result=0;
	double currentVariablePower=1.;
	for(unsigned int i=0;i<2*splineOrder-1;i++)
		{
		result+=splineErrorCoefficients[i]*currentVariablePower;
		currentVariablePower*=variable;
		}
	return sqrt(result);
	
}


double spline::splineIntegral(slotBounds theBounds) const{
	
	//think about what to do when out of bounds: error probably not warranted; might return zero or actual spline value
	double result=0;
	for(unsigned int i=0;i<2*splineOrder-1;i++)
	{
		result+=splineCoefficients[i]*splineBasisFunction(theBounds, i);
	}
	return result;
	
}

void spline::printSpline() const
{
	bounds.printBoundsInfo();
	cout << left << setw(12) << "spline order" << '\t' << setw(12) << "chi^2" << '\t'  << setw(12) << "dof" << '\t' << setw(12) << "chi^2/dof" << '\t' << setw(12) << "goodness of fit Q" << endl;
	cout << left << setw(12) << setprecision(10) << splineOrder << '\t' << setw(12) << chiSquared << '\t' << setw(12) << degreesOfFreedom << '\t' << setw(12) << getChisquared() << '\t' << setw(12) << goodnessOfFitQ() << endl;
	cout << "spline coefficients and error coefficients" << endl;
	for(unsigned int i=0;i<splineOrder;i++) cout << setprecision(10) << splineCoefficients[i] << '\t';
	cout << endl;
	for(unsigned int i=0;i<2*splineOrder-1;i++) cout << setprecision(10) << splineErrorCoefficients[i] << '\t';
	cout << endl;
}

splineArray::splineArray(vector< spline* > theSplines)
{
	
	lowerBound=(theSplines[0] -> getBounds()).getLowerBound();
	upperBound=(theSplines[theSplines.size()-1] -> getBounds()).getUpperBound();
	splineOrder=(theSplines[0] -> getSplineOrder()); //should maybe be total of parameters
	totalChiSquared=0;
	totalDegreesOfFreedom=0;
	splines.push_back(theSplines[0]);
	for(unsigned int i=1;i<theSplines.size();i++)
	{
		splines.push_back(theSplines[i]);
		intervalBoundaries.push_back((theSplines[i] -> getBounds()).getLowerBound());
	}

}

void splineArray::updateProperties(double theTotalChisq, int theTotalDOF)
{
totalChiSquared=theTotalChisq;
totalDegreesOfFreedom=theTotalDOF;
}

spline* splineArray::getSpline(unsigned int whichSpline) const
{
if(whichSpline>=splines.size()) {cout << "ERROR in getSpline, " << whichSpline << " does not exist" << endl; exit(EXIT_FAILURE);}
return splines[whichSpline];
}

double splineArray::goodnessOfFitQ() const
{
	return gsl_sf_gamma_inc_Q (double(totalDegreesOfFreedom)/2., totalChiSquared/2.);
}


double splineArray::splineValue(double variable) const
{
	bool foundSpline; double result;
	for(unsigned int i=0;i<splines.size();i++)
		{
		foundSpline=(splines[i] -> getBounds()).checkIfInBounds(variable);
		if(foundSpline) {result = (splines[i] -> splineValue(variable)); break;}
		}
	return result;
}

double splineArray::splineError(double variable) const
{
	bool foundSpline; double result;
	for(unsigned int i=0;i<splines.size();i++)
	{
		foundSpline=(splines[i] -> getBounds()).checkIfInBounds(variable);
		if(foundSpline) {result = (splines[i] -> splineError(variable)); break;}
	}
	return result;
}

void splineArray::printSplineArrayInfo() const
{
	cout << "Lower Bound: " << lowerBound << endl;
	cout << "Upper Bound: " << upperBound << endl;
	cout << left << setw(12) << "spline order" << '\t' << setw(12) << "chi^2" << '\t'  << setw(12) << "dof" << '\t' << setw(12) << "chi^2/dof" << '\t' << setw(12) << "goodness of fit Q" << endl;
	cout << left << setw(12) << setprecision(10) << splineOrder << '\t' << setw(12) << totalChiSquared << '\t' << setw(12) << totalDegreesOfFreedom << '\t' << setw(12) << getChisquared() << '\t' << setw(12) << goodnessOfFitQ() << endl;
}

void splineArray::printSplines() const
{
	for(unsigned int i=0;i<splines.size();i++)
		{
		cout << "Spline " << i << "  -------------------------------------------------------------------------" << endl;
		splines[i] -> printSpline();
		cout << "-------------------------------------------------------------------------------------------" << endl;
		cout << endl;
		}
}


double splineBasisFunction(slotBounds theBounds, unsigned int functionOrder)
{
	return (pow(theBounds.getUpperBound(),functionOrder+1)-pow(theBounds.getLowerBound(),functionOrder+1))/double(functionOrder+1);
}

double* binomialMatrix(int vectorSize)
{
	double* matrix = new double[vectorSize*(vectorSize-1)];
	for(int k=0;k<vectorSize-1;k++)
		for(int n=0;n<vectorSize;n++)
			{
			if(k>n) matrix[k*vectorSize+n]=0;
			else {
				matrix[k*vectorSize+n]=1.;
				for(int i=n-k+1;i<=n;i++) matrix[k*vectorSize+n]*=double(i);
				for(int i=1;i<=k;i++) matrix[k*vectorSize+n]*=1./double(i);
				}
			}

return matrix;
}

void solveLinearEquation(double* matrix, double* vec, unsigned int vectorSize)
{
	gsl_matrix_view m = gsl_matrix_view_array (matrix, vectorSize, vectorSize);
	gsl_vector_view b = gsl_vector_view_array (vec, vectorSize);
	gsl_vector *x = gsl_vector_alloc (vectorSize);
	int s;
	gsl_permutation * p = gsl_permutation_alloc (vectorSize);
	
	gsl_linalg_LU_decomp (&m.matrix, p, &s);
	gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);
	
	gsl_permutation_free (p);
	
	for(unsigned int i=0;i<vectorSize;i++) vec[i] = gsl_vector_get (x, i);

	gsl_vector_free (x); 
}

double solveSVD(gsl_matrix * U, gsl_matrix * V, gsl_vector* b, gsl_vector* diagonal, gsl_vector* x)
	{
	unsigned int matrixRows = U->size1;
	unsigned int matrixCols = U->size2;

	gsl_matrix * helper = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * work = gsl_vector_alloc (matrixCols);

	gsl_matrix * A = gsl_matrix_alloc (matrixRows, matrixCols);
	gsl_matrix_memcpy (A, U);

	//SVD, original matrix is replaced by U
	if (matrixRows<10*matrixCols) gsl_linalg_SV_decomp (U, V, diagonal, work);		//Golub-Reinsch SVD algorithm
	else gsl_linalg_SV_decomp_mod (U, helper, V, diagonal, work);				//modified Golub-Reinsch algorithm, faster for M>>N

	//gsl_linalg_SV_decomp_jacobi (U, V, diagonal);						//another method: one-sided Jacobi orthogonalization
	//The Jacobi method can compute singular values to higher relative accuracy than Golub-Reinsch algorithms

	//nullify small diagonal values (they are sorted in decreasing order)
	double maximalValue = gsl_vector_get(diagonal,0);
	double threshold=maximalValue*SVD_THRESHOLD;
	for(unsigned int i=1;i<matrixCols;i++)
	{
		if(gsl_vector_get(diagonal,i)<threshold) gsl_vector_set(diagonal,i,0);
	}
	
	gsl_linalg_SV_solve (U, V, diagonal, b, x);

	//compute chi^2 as |Ax-b|^2
	gsl_blas_dgemv (CblasNoTrans, 1, A, x, -1, b);
	double chisq;
	gsl_blas_ddot (b, b, &chisq);
	
	return chisq;
	}



double* secondDerivativeMatrix(unsigned int vectorSize)
{
	double* matrix = new double[vectorSize*vectorSize];
	for(int i=0;i<vectorSize;i++)
		for(int j=0;j<vectorSize;j++)
		{
			if(i==j)
				{
				if((i==0)||(i==vectorSize-1))matrix[i*vectorSize+j]=1;
				else if((i==1)||(i==vectorSize-2))matrix[i*vectorSize+j]=5;
				else matrix[i*vectorSize+j]=6;
				}
			else if(abs(i-j)==1)
				{
				if((i+j==1)||(i+j==2*vectorSize-3)) matrix[i*vectorSize+j]=-2;
				else matrix[i*vectorSize+j]=-4;
				}
			else if(abs(i-j)==2) matrix[i*vectorSize+j]=1;
			else matrix[i*vectorSize+j]=0;
		}
		
	return matrix;
}

double* designMatrix(vector<basisSlot*> slotArray, unsigned int splineOrder){
	
	unsigned int numberOfSlots=slotArray.size();
	double* matrix = new double[numberOfSlots*splineOrder];
	
	for(unsigned int i=0;i<numberOfSlots;i++)
		for(unsigned int j=0;j<splineOrder;j++)
		{
			matrix[i*splineOrder+j]=splineBasisFunction(slotArray[i] -> getBounds(),j)/(slotArray[i] -> sampledIntegralError());
		}

	return matrix;
	
}

double* designMatrixProduct(double* matrix, unsigned int numberOfSlots, unsigned int splineOrder)
{
	gsl_matrix_view A = gsl_matrix_view_array (matrix, numberOfSlots, splineOrder);
	gsl_matrix* C = gsl_matrix_calloc( splineOrder, splineOrder);
	
	//result symmetric, so C only stores half the values
	gsl_blas_dsyrk (CblasUpper, CblasTrans, 1, &A.matrix, 0, C);
	
	double* result = new double[splineOrder*splineOrder];
	
	for(unsigned int i=0;i<splineOrder;i++)
		for(unsigned int j=i;j<splineOrder;j++)
		{
			result[i*splineOrder+j]=gsl_matrix_get(C,i,j);
			if(i!=j) result[j*splineOrder+i]=gsl_matrix_get(C,i,j);
		}
	
	return result;
}


double* integralVector(std::vector< basisSlot* > slotArray){
	
	unsigned int numberOfSlots=slotArray.size();
	double* vec = new double[numberOfSlots];
	
	for(unsigned int i=0;i<numberOfSlots;i++) vec[i]=(slotArray[i] -> sampledIntegral())/(slotArray[i] -> sampledIntegralError());
	
	return vec;
}

double* integralVectorProduct(double* matrix, double* vec, unsigned int numberOfSlots, unsigned int splineOrder){
	
	gsl_vector_view b = gsl_vector_view_array (vec, numberOfSlots);
	gsl_matrix_view A = gsl_matrix_view_array (matrix, numberOfSlots, splineOrder);
	gsl_vector * y = gsl_vector_calloc (splineOrder);
	
	gsl_blas_dgemv (CblasTrans, 1, &A.matrix, &b.vector, 0, y);
	
	double* result = new double[splineOrder];
	for(unsigned int i=0;i<splineOrder;i++) result[i]=gsl_vector_get(y,i);
	
	return result;
	
}





