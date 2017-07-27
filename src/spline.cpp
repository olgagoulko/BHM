#include "basic.hpp"
#include "slot.hpp"
#include "spline.hpp"

using namespace std; 


splinePiece::splinePiece(slotBounds theBounds, unsigned int theSplineOrder)
	{
	bounds=theBounds;
	splineOrder=theSplineOrder;
	splineCoefficients.resize(splineOrder);
	splineErrorCoefficients.resize(2*splineOrder-1);
	}

void splinePiece::setSplinePiece(vector< double > theCoefficients, vector< double > theErrorCoefficients)
	{
	splineCoefficients=theCoefficients;
	splineErrorCoefficients=theErrorCoefficients;
	}


bool splinePiece::checkIntervalAcceptance(vector< vector< basisSlot* > > currentAnalysisBins, double fitAcceptanceThreshold, double chisqArrayElement, unsigned int intervalOrder, bool checkIntervals) const
	{
	bool currentSplineGood=true;
	if(intervalOrder>currentAnalysisBins.size())
		{cout << "ERROR in checkIntervalAcceptance, intervalOrder " << intervalOrder << " is larger than size of analysis bins vector " << currentAnalysisBins.size() << endl; exit(EXIT_FAILURE);}
	for(unsigned int j=0;j<currentAnalysisBins.size()-intervalOrder;j++)
		{
		//test all bins fully within interval, order by order
		double currentChisq=0;
		unsigned int numberSlotsAtCurrentLevel=0;
		for(unsigned int k=0;k<currentAnalysisBins[j+intervalOrder].size();k++)
			{
			basisSlot* currentSlot = currentAnalysisBins[j+intervalOrder][k];
			if((currentSlot -> getBounds()).overlapping(this -> getBounds()))
				{
				numberSlotsAtCurrentLevel++;
				double currentSplineIntegral = this -> splineIntegral(currentSlot -> getBounds());
				currentSplineIntegral-=currentSlot -> sampledIntegral();
				currentSplineIntegral*=1./(currentSlot -> sampledIntegralError());
				currentChisq+=currentSplineIntegral*currentSplineIntegral;
				}
			}
		
		if(numberSlotsAtCurrentLevel>pow(2,j)/2.)
			{
			currentChisq*=1./double(numberSlotsAtCurrentLevel);
			double delta=1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))-currentChisq; if(delta<0) delta=0;
			if(delta<chisqArrayElement) chisqArrayElement=delta;
			if(checkIntervals) cout << intervalOrder+j << '\t' << numberSlotsAtCurrentLevel << '\t' << currentChisq << '\t' << 1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel)) << endl;
			if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {currentSplineGood=false; break;}
			}
		else 
			{
			if(checkIntervals) cout << "Only " << numberSlotsAtCurrentLevel << " good bins, not enough for evaluation" << endl;
			break;
			}
		}
		
	if(checkIntervals) {cout << "This interval fit is "; if(!currentSplineGood) cout << "not "; cout << "good" << endl;}
	return currentSplineGood;
	}


double splinePiece::splineValue(double variable) const
	{
	//think about what to do when out of bounds: return zero or actual spline value?
	double result=0;
	double currentVariablePower=1.;
	for(unsigned int i=0;i<splineOrder;i++)
		{
		result+=splineCoefficients[i]*currentVariablePower;
		currentVariablePower*=variable;
		}
	return result;
	}

double splinePiece::splineError(double variable) const
	{
	//think about what to do when out of bounds: error probably not warranted; return zero or actual value?
	double result=0;
	double currentVariablePower=1.;
	for(unsigned int i=0;i<2*splineOrder-1;i++)
		{
		result+=splineErrorCoefficients[i]*currentVariablePower;
		currentVariablePower*=variable;
		}
	if(result<0) result=0;
	return sqrt(result);
	}


double splinePiece::splineIntegral(slotBounds theBounds) const
	{
	//think about what to do when out of bounds: return zero or actual value?
	double result=0;
	for(unsigned int i=0;i<splineOrder;i++)
		{
		result+=splineCoefficients[i]*splineBasisFunction(theBounds, i);
		}
	return result;
	}

double splinePiece::splineDerivative(double variable, unsigned int derivativeOrder) const
	{
	double result=0;
	double currentVariablePower=1.;
	for(unsigned int i=derivativeOrder;i<splineOrder;i++)
		{
		int powerCoeff=1; unsigned int j=i; for(unsigned int k=0;k<derivativeOrder;k++) {powerCoeff*=j;j--;}
		result+=splineCoefficients[i]*currentVariablePower*powerCoeff;
		currentVariablePower*=variable;
		}
	return result;
	}


std::ostream& splinePiece::printSplinePiece(std::ostream& strm) const
	{
	bounds.printBoundsInfo(strm, CONCISE);
	// strm << "spline piece order: " << splineOrder << endl;
	// strm << "spline piece coefficients and error coefficients" << endl;
	for(unsigned int i=0;i<splineOrder;i++) strm << setprecision(10) << splineCoefficients[i] << " ";
	strm << endl;
	for(unsigned int i=0;i<2*splineOrder-1;i++) strm << setprecision(10) << splineErrorCoefficients[i] << " ";
	strm << endl;
        return strm;
	}

splineArray::splineArray()
	{
	lowerBound=0;
	upperBound=0;
	splineOrder=4;
	splines.resize(0);
	intervalBoundaries.resize(0);
	levelsChiSquared.resize(0);
	levelsDegreesOfFreedom.resize(0);
	acceptableSpline=false;
	}

splineArray::splineArray(vector< splinePiece* > theSplines)
	{
	lowerBound=(theSplines[0] -> getBounds()).getLowerBound();
	upperBound=(theSplines[theSplines.size()-1] -> getBounds()).getUpperBound();
	splineOrder=(theSplines[0] -> getSplineOrder());
	intervalBoundaries.resize(0);
	levelsChiSquared.resize(0);
	levelsDegreesOfFreedom.resize(0);
	acceptableSpline=false;
	splines.resize(0);
	splines.push_back(theSplines[0]);
	for(unsigned int i=1;i<theSplines.size();i++)
		{
		splines.push_back(theSplines[i]);
		intervalBoundaries.push_back((theSplines[i] -> getBounds()).getLowerBound());
		}
	}

splineArray::~splineArray()
	{
	for (vector<splinePiece*>::iterator i = splines.begin(); i != splines.end(); i++ ) delete *i;
	splines.clear();
	}
	
splineArray& splineArray::operator=(const splineArray& toBeAssigned)
	{
	if (this != &toBeAssigned)
		{
		for (vector<splinePiece*>::iterator i = splines.begin(); i != splines.end(); i++ ) delete *i;
		splines.clear();
		
		for (unsigned int i=0;i<toBeAssigned.splines.size();i++)
			{
			splinePiece* theSpline=toBeAssigned.splines[i] -> Clone();
			splines.push_back(theSpline);
			}
		lowerBound=(splines[0] -> getBounds()).getLowerBound();
		upperBound=(splines[splines.size()-1] -> getBounds()).getUpperBound();
		splineOrder=(splines[0] -> getSplineOrder());
		updateLevelProperties(toBeAssigned.levelsChiSquared, toBeAssigned.levelsDegreesOfFreedom);
		acceptableSpline=toBeAssigned.acceptableSpline;
		intervalBoundaries=toBeAssigned.intervalBoundaries;
		}
	return *this;
	}

splineArray::splineArray(const splineArray& toBeCopied)
	{
	for (unsigned int i=0;i<toBeCopied.splines.size();i++)
		{
		splinePiece* theSpline=toBeCopied.splines[i] -> Clone();
		splines.push_back(theSpline);
		}
	
	lowerBound=(splines[0] -> getBounds()).getLowerBound();
	upperBound=(splines[splines.size()-1] -> getBounds()).getUpperBound();
	splineOrder=(splines[0] -> getSplineOrder());
	updateLevelProperties(toBeCopied.levelsChiSquared, toBeCopied.levelsDegreesOfFreedom);
	acceptableSpline=toBeCopied.acceptableSpline;
	intervalBoundaries=toBeCopied.intervalBoundaries;
	}


void splineArray::updateLevelProperties(std::vector<double> theChisq, std::vector<int> theDOF)
	{
	levelsChiSquared.resize(0);
	levelsDegreesOfFreedom.resize(0);
	for(unsigned int i=0;i<theChisq.size();i++)
		{
		levelsChiSquared.push_back(theChisq[i]);
		levelsDegreesOfFreedom.push_back(theDOF[i]);
		}
	}

void splineArray::updateGoodness(bool acceptable)
	{
	acceptableSpline=acceptable;
	}


splinePiece* splineArray::getSplinePiece(unsigned int whichPiece) const
	{
	if(whichPiece>=splines.size()) {cout << "ERROR in getSplinePiece, piece number " << whichPiece << " does not exist" << endl; exit(EXIT_FAILURE);}
	return splines[whichPiece];
	}


bool splineArray::checkOverallAcceptance(double fitAcceptanceThreshold) const
	{
	bool overallAcceptance=true;
	for(unsigned int i=0;i<levelsChiSquared.size();i++)
		if( (levelsChiSquared[i]-1)/sqrt(2./double(levelsDegreesOfFreedom[i]))>fitAcceptanceThreshold ) overallAcceptance=false;
	return overallAcceptance;
	}


vector< slotBounds > splineArray::getBounds() const
	{
	vector<slotBounds> theBounds;
	for(unsigned int i=0;i<splines.size();i++) theBounds.push_back(splines[i] -> getBounds());
	return theBounds;
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

double splineArray::splineDerivative(double variable, unsigned int derivativeOrder) const
	{
	bool foundSpline;
        double result=0;
	for(unsigned int i=0;i<splines.size();i++)
		{
		foundSpline=(splines[i] -> getBounds()).checkIfInBounds(variable);
		if(foundSpline) {result = (splines[i] -> splineDerivative(variable,derivativeOrder)); break;}
		}
	return result;
	}

std::ostream& splineArray::printSplineArrayInfo(std::ostream& strm) const
	{
	strm << "##Lower Bound: " << lowerBound << endl;
	strm << "##Upper Bound: " << upperBound << endl;
	strm << "##Spline order: " << splineOrder << endl;
	strm << left << setw(4) << "## level" << '\t' << setw(8) << "n" << '\t'  << setw(12) << "chi_n^2/n" << '\t' << setw(12) << "sqrt(2/n)" << '\t' << setw(12) << "deviation in sqrt(2/n) units" << endl;
	for(unsigned int i=0;i<levelsChiSquared.size();i++)
		strm << "#" << left << setw(4) << setprecision(10) << i << '\t' << setw(8) <<  levelsDegreesOfFreedom[i] << '\t' << setw(12) << levelsChiSquared[i] << '\t' << setw(12) << sqrt(2./double(levelsDegreesOfFreedom[i])) << '\t' << setw(12) << max(0.0,(levelsChiSquared[i]-1)/sqrt(2./double(levelsDegreesOfFreedom[i]))) << endl;

        return strm;
	}

std::ostream& splineArray::printSplines(std::ostream& strm) const
	{
	strm << lowerBound << " " << upperBound << " " << (splineOrder-1) << endl;
        strm << splines.size() << endl;
	for(unsigned int i=0;i<splines.size();i++)
		{
		strm << "# spline piece " << i << endl;
		splines[i] -> printSplinePiece(strm);
		}
        return strm;
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

double* binomialVector(int vectorSize)
	{
	double* vector = new double[vectorSize];
	vector[0]=1;
	for(int k=1;k<vectorSize;k++)
		{
		vector[k]=vector[k-1]*(vectorSize-k)/double(k);
		}
		
	return vector;
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

void solveSVD(gsl_matrix * A, gsl_matrix * V, gsl_vector* b, gsl_vector* diagonal, gsl_vector* x)
	{
	unsigned int matrixRows = A->size1;
	unsigned int matrixCols = A->size2;

	gsl_matrix * helper = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * work = gsl_vector_alloc (matrixCols);

	gsl_matrix * U = gsl_matrix_alloc (matrixRows, matrixCols);
	gsl_matrix_memcpy (U, A);

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

	//compute chi^2 as |Ax-b|^2, b is replaced by Ax-b
	gsl_blas_dgemv (CblasNoTrans, 1, A, x, -1, b);
	//double chisq;
	//gsl_blas_ddot (b, b, &chisq);
	
	gsl_matrix_free(helper); gsl_vector_free(work); gsl_matrix_free(U);
	}


double* designMatrix(vector<basisSlot*> slotArray, unsigned int splineOrder)
	{
	unsigned int numberOfSlots=slotArray.size();
	double* matrix = new double[numberOfSlots*splineOrder];
	
	for(unsigned int i=0;i<numberOfSlots;i++)
		for(unsigned int j=0;j<splineOrder;j++)
			{
			matrix[i*splineOrder+j]=splineBasisFunction(slotArray[i] -> getBounds(),j)/(slotArray[i] -> sampledIntegralError());
			}

	return matrix;
	}


double* designMatrix(vector< vector<basisSlot*> > slotArray, unsigned int splineOrder)
	{
	unsigned int numberOfSlots=0;	
	for(unsigned int i=0;i<slotArray.size();i++) numberOfSlots+=slotArray[i].size();
	double* matrix = new double[numberOfSlots*splineOrder];
	
	unsigned int counter=0; double weightFactor;
	for(unsigned int i=0;i<slotArray.size();i++)
		{
		weightFactor=sqrt(double(pow(2,i)));
		for(unsigned int j=0;j<slotArray[i].size();j++)
			{
			double currentError=slotArray[i][j] -> sampledIntegralError();
			for(unsigned int k=0;k<splineOrder;k++)
				{
				if( currentError < VERY_SMALL_NUMBER) matrix[counter*splineOrder+k]=0;
				else matrix[counter*splineOrder+k]=splineBasisFunction(slotArray[i][j] -> getBounds(),k)/currentError/weightFactor;
				}
			counter++;
			}
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
	
	gsl_matrix_free(C);
	
	return result;
	}


double* integralVector(vector< basisSlot* > slotArray)
	{
	unsigned int numberOfSlots=slotArray.size();
	double* vec = new double[numberOfSlots];
	
	for(unsigned int i=0;i<numberOfSlots;i++) vec[i]=(slotArray[i] -> sampledIntegral())/(slotArray[i] -> sampledIntegralError());
	
	return vec;
	}

double* integralVector(vector< vector<basisSlot*> > slotArray)
	{
	unsigned int numberOfSlots=0;	
	for(unsigned int i=0;i<slotArray.size();i++) numberOfSlots+=slotArray[i].size();
	double* vec = new double[numberOfSlots];
	
	unsigned int counter=0; double weightFactor;
	for(unsigned int i=0;i<slotArray.size();i++)
		{
		weightFactor=sqrt(double(pow(2,i)));
		for(unsigned int j=0;j<slotArray[i].size();j++)
			{
			if( (slotArray[i][j] -> sampledIntegralError()) < VERY_SMALL_NUMBER) vec[counter]=0;
			else vec[counter]=(slotArray[i][j] -> sampledIntegral())/(slotArray[i][j] -> sampledIntegralError())/weightFactor;
			counter++;
			}
		}
	
	return vec;
	}

double* integralVectorProduct(double* matrix, double* vec, unsigned int numberOfSlots, unsigned int splineOrder)
	{
	gsl_vector_view b = gsl_vector_view_array (vec, numberOfSlots);
	gsl_matrix_view A = gsl_matrix_view_array (matrix, numberOfSlots, splineOrder);
	gsl_vector * y = gsl_vector_calloc (splineOrder);
	
	gsl_blas_dgemv (CblasTrans, 1, &A.matrix, &b.vector, 0, y);
	
	double* result = new double[splineOrder];
	for(unsigned int i=0;i<splineOrder;i++) result[i]=gsl_vector_get(y,i);
	
	gsl_vector_free(y);
	return result;
	}

	
bool isDataConsistentWithZero(vector< vector< basisSlot* > > slotArray)
	{
	bool consistentWithZero=true; unsigned int slotCounter; double zeroChiSqVec; int consistentWithZeroCounter=0;
	for(unsigned int j=0;j<slotArray.size();j++)
		{
		slotCounter=0;
		zeroChiSqVec=0;
		for(unsigned int i=0;i<slotArray[j].size();i++)
			{
			basisSlot* currentSlot = slotArray[j][i];
			double currentError=currentSlot -> sampledIntegralError();
			if( (currentSlot -> enoughSampled()) && (currentError > VERY_SMALL_NUMBER))
				{
				zeroChiSqVec+=(currentSlot -> sampledIntegral())*(currentSlot -> sampledIntegral())/currentError/currentError;
				slotCounter++;
				}
			}
		if(2*slotCounter>slotArray[j].size())
			{
			zeroChiSqVec=zeroChiSqVec/double(slotCounter);
			double deviation=(zeroChiSqVec-1)/sqrt(2./double(slotCounter));
			if(deviation>4) consistentWithZero=false;
			else if(deviation>3) consistentWithZeroCounter+=4;
			else if(deviation>2) consistentWithZeroCounter+=2;
			}
		}
	if(consistentWithZeroCounter>=8) consistentWithZero=false;
	return consistentWithZero;
	}




