#include "basic.hpp"
#include "slot.hpp"
#include "matrix.hpp"
#include "histogram.hpp"

using namespace std; 

vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps, int totalNumOfBasisFn)
{
	if(numberOverlaps<1) numberOverlaps = 10;
	
	if(totalNumOfBasisFn<0) totalNumOfBasisFn=0;
	
	if(maxVar<minVar) {double newMin=maxVar; maxVar=minVar; minVar=newMin;}
	
	double range = maxVar-minVar;
	if(range<VERY_SMALL_NUMBER) {cout << "ERROR: sampling interval too small; min = " << minVar << ", max = " << maxVar << ", range = " << range << endl; exit(EXIT_FAILURE);}
	
	if((slotWidth>range)||(slotWidth<0)) slotWidth=range/double(numberOverlaps);
	
	double offset = slotWidth/double(numberOverlaps);
	vector<basisSlot*> basisSlotVector;
	
	int counter=0;
	for(int i=1;i<numberOverlaps;i++)
		{
		slotBounds currentBounds(minVar, minVar+i*offset);
		basisSlotVector.push_back(new taylorSlot(currentBounds, totalNumOfBasisFn));
		}
	while(minVar+(counter*offset)+slotWidth<maxVar)
		{
		slotBounds currentBounds(minVar+(counter*offset), minVar+(counter*offset)+slotWidth);
		basisSlotVector.push_back(new taylorSlot(currentBounds, totalNumOfBasisFn));
		counter++;
		if((minVar+(counter*offset)+slotWidth>=maxVar)&&(minVar+(counter*offset)<maxVar)) 
			{
			if(numberOverlaps<2)
				{
				slotBounds finalBounds(minVar+counter*offset, maxVar);
				basisSlotVector.push_back(new taylorSlot(finalBounds, totalNumOfBasisFn));
				}
			else
				{
				for(int i=0;i<numberOverlaps-1;i++)
					{
					if(minVar+(counter+i)*offset<maxVar)
						{
						slotBounds finalBounds(minVar+(counter+i)*offset, maxVar);
						basisSlotVector.push_back(new taylorSlot(finalBounds, totalNumOfBasisFn));
						}
					}
				}
			}
		}
	
	return basisSlotVector;
}


histogramBasis::histogramBasis(vector<basisSlot*> theBasisSlots)
{
	
	valuesOutsideBounds = new excessBin();
	lowerBound=(theBasisSlots[0] -> getBounds()).getLowerBound(); upperBound=0; noUpperBound=false;
	for(unsigned int i=0;i<theBasisSlots.size();i++)
		{
		basisSlots.push_back(theBasisSlots[i]);
		slotBounds currentBounds=(theBasisSlots[i] -> getBounds());
		if(currentBounds.getLowerBound()<lowerBound) lowerBound=currentBounds.getLowerBound();
		if(currentBounds.getIsInfinite()) noUpperBound=true;
		else if(currentBounds.getUpperBound()>upperBound) upperBound=currentBounds.getUpperBound();
		}

}

void histogramBasis::appendSlot(basisSlot* theSlot)
{
	basisSlots.push_back(theSlot);
	if((theSlot -> getBounds()).getLowerBound()<lowerBound) lowerBound=(theSlot -> getBounds()).getLowerBound();
	if((theSlot -> getBounds()).getIsInfinite()) noUpperBound=true;
	else if((theSlot -> getBounds()).getUpperBound()>upperBound) upperBound=(theSlot -> getBounds()).getUpperBound();
}

basisSlot* histogramBasis::combinedSlot(unsigned int startPoint, unsigned int endPoint) const
{
	if((startPoint>endPoint) || (endPoint>=basisSlots.size())) {cout << "ERROR in combinedSlot" << endl; exit(EXIT_FAILURE);}
	else if(startPoint>endPoint) {unsigned int save=endPoint; endPoint=startPoint; startPoint=save;}
	
	for(unsigned int i=startPoint;i<=endPoint;i++)
		{
		for(unsigned int j=i+1;j<=endPoint;j++)
			{
			if((basisSlots[i] -> getBounds()).overlapping((basisSlots[j] -> getBounds()))) 
				{
				cout << "ERROR: trying to combine overlapping slots " << i << " and " << j << endl; 
				basisSlots[i] -> printSlotInfo(); basisSlots[j] -> printSlotInfo();
				exit(EXIT_FAILURE);
				}
			}
		}
	
	slotBounds theBounds( basisSlots[startPoint] -> getBounds().getLowerBound(), basisSlots[endPoint] -> getBounds().getUpperBound());
	basisSlot * combined = new basisSlot(theBounds);
	for(unsigned int i=startPoint;i<=endPoint;i++)
		{
		combined -> combineWithSlot(basisSlots[i]);
		}

	return combined;
}

basisSlot* histogramBasis::getSlot(unsigned int whichSlot) const
{
if(whichSlot>=basisSlots.size()) {cout << "ERROR in getSlot, " << whichSlot << " does not exist" << endl; exit(EXIT_FAILURE);}
return basisSlots[whichSlot];
}


void histogramBasis::sample(double variable, double valueToSample)
{
	bool isOutsideBounds=true;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) ) {basisSlots[i] -> sample(variable,valueToSample); isOutsideBounds=false;}
		}
	if(isOutsideBounds) valuesOutsideBounds -> sample(variable, valueToSample);
}

pair<double,double> histogramBasis::sampledFunctionValueAverage(double variable) const
{
	vector<double> sampledValue;
	int numberSlots=0;
	double average=0;
	double error=0;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) ) 
			{
			sampledValue.push_back(basisSlots[i] -> sampledFunctionValue(variable));
			average+=sampledValue[numberSlots];
			numberSlots++;
			}
		}
	
	if(numberSlots==0) average=valuesOutsideBounds -> getExcessValues();
	else if(numberSlots>1)
		{
		average=average/double(numberSlots);
		double roundofferror=0;
		double delta;
		for(int i=0;i<numberSlots;i++)
			{
			delta=sampledValue[i]-average;
			error+=delta*delta;
			roundofferror+=delta;
			}
		roundofferror=roundofferror*roundofferror/double(numberSlots);
		error=error-roundofferror;
		if(error<0) error=0;
		error=sqrt(error/double(numberSlots*(numberSlots-1)));
		}
		
	pair<double,double> thisPair(average,error);
	return thisPair;
}

pair<double,double> histogramBasis::sampledFunctionValueWeightedAverage(double variable) const
{
	double sampledValue=0; double value;
	double sumVariance=0; double variance;
	double error;
	int numberSlots=0;
	for(unsigned int i=0;i<basisSlots.size();i++)
	{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) )
			{
			value = basisSlots[i] -> sampledFunctionValue(variable);
			variance = basisSlots[i] -> sampledFunctionError(variable);
			if(variance>0)
				{
				variance=variance*variance;
				sumVariance+=1./variance;
				sampledValue += value/variance;
				}
			numberSlots++;
			}
	}
	
	if(numberSlots==0) sampledValue=valuesOutsideBounds -> getExcessValues();
	else if(sumVariance==0) {error=0;}
	else {sampledValue=sampledValue/sumVariance; error=sqrt(1./sumVariance);}
	
	pair<double,double> result(sampledValue,error);
	
	return result;
}

void histogramBasis::scale(double norm)
{
	for(unsigned int i=0;i<basisSlots.size();i++)
	{
		basisSlots[i] -> scale(norm);
	}
}


spline histogramBasis::oneSplineFit(unsigned int splineOrder)
{
	unsigned int dataSize = basisSlots.size();
	double* U = designMatrix(basisSlots, splineOrder);
	double* b = integralVector(basisSlots);
	gsl_matrix_view gslU = gsl_matrix_view_array (U, dataSize, splineOrder);
	gsl_vector_view gslb = gsl_vector_view_array (b, dataSize);
	
	gsl_matrix * V = gsl_matrix_alloc (splineOrder, splineOrder);
	gsl_vector * diagonal = gsl_vector_alloc (splineOrder);
	gsl_vector *x = gsl_vector_alloc (splineOrder);
	
	double chisq=solveSVD(&gslU.matrix, V, &gslb.vector, diagonal, x);
	
	vector<double> theCoefficients;
	for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,i));
	
	double diagScaling[splineOrder];
	for(unsigned int j=0;j < splineOrder;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./(gsl_vector_get(diagonal,j)*gsl_vector_get(diagonal,j));
		}
		
	vector<double> theErrorCoefficients;
	for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0);
		
	for(unsigned int j=0;j < splineOrder;j++)
		for(unsigned int k=j;k < splineOrder;k++)
			{
			double sum=0;
			for(unsigned int i=0;i<splineOrder;i++) sum+=diagScaling[i]*gsl_matrix_get(V,j,i)*gsl_matrix_get(V,k,i);
			theErrorCoefficients[j+k]+=sum;
			if(j!=k) theErrorCoefficients[j+k]+=sum;
			}
			
	int dof = dataSize-splineOrder;
	
	slotBounds theBounds(lowerBound,upperBound);
	spline solution(theBounds, splineOrder); solution.setSpline(theCoefficients,theErrorCoefficients,chisq,dof);
	return solution;
}

splineArray histogramBasis::splineFit(vector<unsigned int> intervalBoundaries, double gluingFactor, unsigned int splineOrder){
	
	unsigned int numberIntervals=intervalBoundaries.size();
	unsigned int dataSize = basisSlots.size();
	int matrixRows=dataSize;
	if(gluingFactor>0) matrixRows+=(splineOrder-1)*numberIntervals;
	int matrixCols=(numberIntervals+1)*splineOrder;
	
	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows, matrixCols);
	double* b = integralVector(basisSlots);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows);
	for(unsigned int i=0;i<dataSize;i++) gsl_vector_set(gslb,i,b[i]);
	
	//regular part of block design matrix containing data integrals
	unsigned int start=0; unsigned int end;
	int rowOffset, colOffset;
	for(unsigned int i=0;i<=numberIntervals;i++)
		{			
		if(i!=numberIntervals) end = intervalBoundaries[i];
		else end = dataSize;
		rowOffset=start;
		colOffset=i*splineOrder;
		
		vector<basisSlot*> intervalSlots(basisSlots.begin() + start, basisSlots.begin() + end);
		double* U = designMatrix(intervalSlots, splineOrder);
	
		for(unsigned int j=0;j < end-start;j++)
			for(unsigned int k=0; k < splineOrder;k++)
				gsl_matrix_set(fullDesignMatrix,rowOffset+j,colOffset+k,U[j*splineOrder+k]);
	
		start=end;
		}
		
	//extra part of design matrix responsible for matching splines at boundaries up to second derivative
	if(gluingFactor>0)
		{
		rowOffset=dataSize;
		colOffset=0;
		double* binMat=binomialMatrix(splineOrder);
		double theBound;
		for(unsigned int i=0;i<numberIntervals;i++)
			{
			theBound=(basisSlots[intervalBoundaries[i]] -> getBounds()).getLowerBound();
			for(unsigned int j=0;j < splineOrder-1;j++)
				for(unsigned int k=j;k < splineOrder;k++)
					{
					double variable=pow(theBound,k-j);
					gsl_matrix_set(fullDesignMatrix, rowOffset+j, colOffset+k, gluingFactor*binMat[j*splineOrder+k]*variable);
				gsl_matrix_set(fullDesignMatrix, rowOffset+j, colOffset+splineOrder+k, -gluingFactor*binMat[j*splineOrder+k]*variable);
					}
			rowOffset+=splineOrder-1;
			colOffset+=splineOrder;
			}
		}
		
	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
	
	double totalChiSquared=solveSVD(fullDesignMatrix, V, gslb, diagonal, x);
	int totalDegreesOfFreedom=dataSize-4*(numberIntervals+1);
	vector<spline*> splines;
	
	/*for (int i = 0; i < 16; i++)
	{
		for (int j = 0; j < 16; j++)
			cout << fixed << setprecision(4) << gsl_matrix_get (V,i, j) << '\t';
		cout << endl;
	}*/
	
	double diagScaling[(numberIntervals+1)*splineOrder];
	for(unsigned int j=0;j < (numberIntervals+1)*splineOrder;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./gsl_vector_get(diagonal,j);
		}
		
	gsl_matrix * Vscaled = gsl_matrix_alloc (matrixCols, matrixCols);
	for (int i = 0; i < matrixCols; i++)
		for (int j = 0; j < matrixCols; j++) gsl_matrix_set(Vscaled,i,j,gsl_matrix_get(V,i,j)*diagScaling[j]);
	gsl_matrix * fullCovarianceMatrix = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1, Vscaled, 0, fullCovarianceMatrix);
	
	start=0; double currentLowerBound; double currentUpperBound; vector<basisSlot*>:: iterator it;
	for(unsigned int n=0;n<=numberIntervals;n++)
		{
		if(n!=numberIntervals) end = intervalBoundaries[n];
		else end = dataSize;
		it=basisSlots.begin() + start;
		currentLowerBound = (it[0] -> getBounds()).getLowerBound();
		it=basisSlots.begin() + end - 1;
		currentUpperBound = (it[0] -> getBounds()).getUpperBound();
			
		vector<double> theCoefficients;
		for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,n*splineOrder+i));
		
		vector<double> theErrorCoefficients;
		for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0); 
		
		for(unsigned int j=0;j < splineOrder;j++)
			for(unsigned int k=j;k < splineOrder;k++)
				{
				theErrorCoefficients[j+k]+=gsl_matrix_get(fullCovarianceMatrix,n*splineOrder+j,n*splineOrder+k);
				if(j!=k) theErrorCoefficients[j+k]+=gsl_matrix_get(fullCovarianceMatrix,n*splineOrder+j,n*splineOrder+k);
				}

		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, start, end-start);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		int dof = end-start-splineOrder;
		
		slotBounds theBounds(currentLowerBound,currentUpperBound);
		spline * solution = new spline(theBounds, splineOrder); solution -> setSpline(theCoefficients,theErrorCoefficients,chisq,dof);
		splines.push_back(solution);
		
		start=end;
		}
		
	splineArray result(splines);
	result.updateProperties(totalChiSquared,totalDegreesOfFreedom);
	return result;
}






