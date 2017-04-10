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

histogramBasis histogramBasis::coarseGrainedHistogram(int minNumberTimesSampled)
{
	vector<basisSlot*> coarseGrainedSlots;
	unsigned int counter=0; bool enough;
	for(unsigned int i=0;i<basisSlots.size();i++) basisSlots[i] -> updateEnoughSampled(minNumberTimesSampled);
	while(counter<basisSlots.size())
	{
		enough = basisSlots[counter] -> enoughSampled();
		if(enough)
			{
			basisSlot* currentSlot = basisSlots[counter] -> Clone();
			coarseGrainedSlots.push_back(currentSlot);
			counter++;
			}
		else {
			for(unsigned int more=counter+1;more<basisSlots.size();more++)
				{
				basisSlot* combined = combinedSlot(counter, more);
				combined -> updateEnoughSampled(minNumberTimesSampled);
				enough = combined -> enoughSampled();
				if(enough) {coarseGrainedSlots.push_back(combined); counter=more+1; break;}
				}
			if(!enough && (coarseGrainedSlots.size()>0)) 
				{
				slotBounds theBounds( (coarseGrainedSlots.back() -> getBounds()).getLowerBound(), upperBound);
				basisSlot * combined = new basisSlot(theBounds);
				combined -> combineWithSlot(coarseGrainedSlots.back());
				for(unsigned int i=counter;i<basisSlots.size();i++) combined -> combineWithSlot(basisSlots[i]);
				coarseGrainedSlots.back() = combined;
				counter=basisSlots.size();
				}
			}
	}
	
	histogramBasis result(coarseGrainedSlots);
	return result;
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

histogramBasis histogramBasis::scaledHistogram(long norm)
{
	vector<basisSlot*> scaledSlots;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		//basisSlot* currentSlot = new basisSlot(*basisSlots[i]);
		basisSlot* currentSlot = basisSlots[i] -> Clone();
		currentSlot -> scale(norm);
		scaledSlots.push_back(currentSlot);
		}
	
	histogramBasis result(scaledSlots);
	return result;

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

splineArray histogramBasis::matchedSplineFitSmallBins(std::vector< unsigned int > intervalBoundaries, unsigned int splineOrder)
{
	unsigned int numberIntervals=intervalBoundaries.size();
	unsigned int matrixRows = basisSlots.size();
	unsigned int matrixCols=numberIntervals+splineOrder;
	
	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows, matrixCols);
	double* b = integralVector(basisSlots);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows);
	for(unsigned int i=0;i<matrixRows;i++) gsl_vector_set(gslb,i,b[i]);
	
	double* binVec=binomialVector(splineOrder);
	double intervalMatrix[numberIntervals][splineOrder-1];
	double lastIntervalMatrix[numberIntervals][splineOrder-1];
	
	//design matrix
	double* basicDesignMatrix = designMatrix(basisSlots, splineOrder);
	unsigned int start=0; unsigned int end;
	double currentLowerBound=0; double currentUpperBound; vector<basisSlot*>:: iterator it;
	for(unsigned int i=0;i<=numberIntervals;i++)
		{
		if(i==numberIntervals) end = matrixRows;
		else
			{
			end = intervalBoundaries[i];
			it=basisSlots.begin() + end - 1;
			currentUpperBound = (it[0] -> getBounds()).getUpperBound();
		
			for(unsigned int k=0; k < splineOrder-1;k++)
				{
				intervalMatrix[i][k]=(pow(currentUpperBound,splineOrder-1-k)-pow(currentLowerBound,splineOrder-1-k))*binVec[k]*pow(-1,splineOrder-2-k);
				lastIntervalMatrix[i][k]=-pow(currentUpperBound,splineOrder-1-k)*binVec[k]*pow(-1,splineOrder-2-k);
				}
			}
			
		for(unsigned int j=0;j < end-start;j++)
			{
			for(unsigned int k=0; k < splineOrder-1;k++)
				gsl_matrix_set(fullDesignMatrix,start+j,k,basicDesignMatrix[(start+j)*splineOrder+k]);
			
			gsl_matrix_set(fullDesignMatrix,start+j,i+splineOrder-1,basicDesignMatrix[(start+j)*splineOrder+splineOrder-1]);
		
			if(i>0)
				{
				for(unsigned int kk=0; kk < splineOrder-1;kk++)
				    {
				    for(unsigned int k=0; k < i;k++)
					gsl_matrix_set(fullDesignMatrix, start+j, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,start+j,k+splineOrder-1)+intervalMatrix[k][kk]*basicDesignMatrix[(start+j)*splineOrder+kk]);
				    
				    gsl_matrix_set(fullDesignMatrix, start+j, i+splineOrder-1, gsl_matrix_get(fullDesignMatrix,start+j,i+splineOrder-1)+lastIntervalMatrix[i-1][kk]*basicDesignMatrix[(start+j)*splineOrder+kk]);
				    }
				}
			}
			
		currentLowerBound=currentUpperBound;
		start=end;
		}

	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
	
	double totalChiSquared=solveSVD(fullDesignMatrix, V, gslb, diagonal, x);
	int totalDegreesOfFreedom=matrixRows-(numberIntervals+splineOrder-1);
	vector<spline*> splines;
	
	double diagScaling[matrixCols];
	for(unsigned int j=0;j < matrixCols;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./gsl_vector_get(diagonal,j);
		}
	
	gsl_matrix * Vscaled = gsl_matrix_alloc (matrixCols, matrixCols);
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			gsl_matrix_set(Vscaled,i,j,gsl_matrix_get(V,i,j)*diagScaling[j]);
	gsl_matrix * fullCovarianceMatrix = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1, Vscaled, 0, fullCovarianceMatrix);
	double covarianceMatrix[matrixCols][matrixCols];
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			covarianceMatrix[i][j]= gsl_matrix_get (fullCovarianceMatrix, i,j);
		
		/*double correlationMatrix[matrixCols][matrixCols];
		for (unsigned int i = 0; i < matrixCols; i++)
		{
			for (unsigned int j = 0; j < matrixCols; j++)
				{correlationMatrix[i][j]=covarianceMatrix[i][j]/sqrt(covarianceMatrix[i][i]*covarianceMatrix[j][j]);
					if( (i<12) && (j<12) ) cout << setw(16) << fixed << setprecision(8) << correlationMatrix[i][j] << " ";}
			if( i<14 ) cout << endl;
		}*/
	
	start=0;
	vector<double> theCoefficients;
	for(unsigned int n=0;n<=numberIntervals;n++)
		{
		if(n!=numberIntervals) end = intervalBoundaries[n];
		else end = matrixRows;
		it=basisSlots.begin() + start;
		currentLowerBound = (it[0] -> getBounds()).getLowerBound();
		it=basisSlots.begin() + end - 1;
		currentUpperBound = (it[0] -> getBounds()).getUpperBound();
		
		if(n==0) 
			{
			for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,i));
			}
		else
			{
			for(unsigned int i=0;i<splineOrder-1;i++) theCoefficients[i]+=(gsl_vector_get(x,n+splineOrder-1)-theCoefficients[splineOrder-1])*lastIntervalMatrix[n-1][i];
			theCoefficients[splineOrder-1]=gsl_vector_get(x,n+splineOrder-1);
			}
		
		vector<double> theErrorCoefficients;
		for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0); 
		
		if(n>0)
			{
			for(unsigned int j=0;j < splineOrder-1;j++)
				{
				for(unsigned int k=j;k < splineOrder-1;k++)
					{
					covarianceMatrix[j][k]+=(covarianceMatrix[j][splineOrder+n-1]-covarianceMatrix[j][splineOrder-1])*lastIntervalMatrix[n-1][k];
					covarianceMatrix[j][k]+=(covarianceMatrix[k][splineOrder+n-1]-covarianceMatrix[k][splineOrder-1])*lastIntervalMatrix[n-1][j];
					covarianceMatrix[j][k]+=(covarianceMatrix[splineOrder-1][splineOrder-1] + covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - 2*covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][k]*lastIntervalMatrix[n-1][j];
					}
				covarianceMatrix[j][splineOrder-1]=covarianceMatrix[j][splineOrder+n-1];
				covarianceMatrix[j][splineOrder-1]+=(covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][j];
			
				for(unsigned int k=n;k < matrixCols-splineOrder;k++)
					covarianceMatrix[j][k+splineOrder]+=(covarianceMatrix[splineOrder+n-1][k+splineOrder] - covarianceMatrix[splineOrder+n-2][k+splineOrder])*lastIntervalMatrix[n-1][j];
				}
			covarianceMatrix[splineOrder-1][splineOrder-1]=covarianceMatrix[n+splineOrder-1][n+splineOrder-1];
			}
		
		for(unsigned int j=0;j < splineOrder;j++)
			for(unsigned int k=j;k < splineOrder;k++)
			{
				theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				if(j!=k) theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
			}
		
		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, start, end-start);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		int dof = end-start-splineOrder;
		//cout << theCoefficients[splineOrder-1] << '\t' << sqrt(covarianceMatrix[splineOrder-1][splineOrder-1]) << '\t' << chisq/double(dof) << endl;
		
		slotBounds theBounds(currentLowerBound,currentUpperBound);
		spline * solution = new spline(theBounds, splineOrder); solution -> setSpline(theCoefficients,theErrorCoefficients,chisq,dof);
		splines.push_back(solution);
		
		start=end;
		}
	
	/*------------- testing interval distribution ----------------------*/
// 	start=0; 
// 	//ofstream output("corrchisqratio.dat");
// 	for(unsigned int n=0;n<=numberIntervals;n++)
// 		{
// 		if(n!=numberIntervals) end = intervalBoundaries[n];
// 		else end = matrixRows;
// 		
// 		int cut=7;
// 		
// 		cout << "Interval " << n << endl; 
// 		for(int m=-cut;m<=cut;m++) if( (n+m>=0) && (n+m<=numberIntervals) ) cout << left << setw(20) << n+m << setw(2);
// 		cout << endl;
// 	
// 		gsl_matrix * currentDesignMatrix = gsl_matrix_alloc (end-start, splineOrder);
// 		gsl_vector * currentb= gsl_vector_alloc (end-start);
// 		gsl_vector * currentCoeffs=gsl_vector_alloc (splineOrder);
// 		for(unsigned int i=0;i<end-start;i++) 
// 			for(unsigned int j=0;j<splineOrder;j++)
// 				gsl_matrix_set(currentDesignMatrix,i,j,basicDesignMatrix[(start+i)*splineOrder+j]);
// 		double currentChisq;
// 		int dof = end-start-splineOrder;
// 
// 		for(int m=-cut;m<=cut;m++)
// 			{
// 			if( (n+m>=0) && (n+m<=numberIntervals) )
// 				{
// 				for(unsigned int j=0;j<splineOrder;j++) gsl_vector_set(currentCoeffs,j, (splines[n+m] -> getCoefficients())[j]);
// 				for(unsigned int i=0;i<end-start;i++) gsl_vector_set(currentb,i,b[start+i]);
// 				gsl_blas_dgemv (CblasNoTrans, 1, currentDesignMatrix, currentCoeffs, -1, currentb);
// 				gsl_blas_ddot (currentb, currentb, &currentChisq);
// 				cout << setw(10) << setprecision(2) << currentChisq << setw(10) << currentChisq/double(dof) << setw(2);
// 				/*
// 				if(m<0) {
// 					double variable=(splines[n+m] -> getBounds()).getLowerBound(); double curError=splines[n+m] -> splineError(variable); if(curError<VERY_SMALL_NUMBER) curError=1e-6;
// 					output << currentChisq/double(dof) << '\t' << abs(((splines[n]-> splineValue(variable)) - (splines[n+m]-> splineValue(variable)))/curError) << endl;
// 					}
// 				if(m>0) {
// 					double variable=(splines[n+m] -> getBounds()).getUpperBound();  double curError=splines[n+m] -> splineError(variable); if(curError<VERY_SMALL_NUMBER) curError=1e-6;
// 					output << currentChisq/double(dof) << '\t' << abs(((splines[n]-> splineValue(variable)) - (splines[n+m]-> splineValue(variable)))/curError) << endl;
// 					}
// 				*/
// 				}
// 			}
// 		cout << endl;
// 		
// 		for(int m=-cut;m<=cut;m++) if( (n+m>=0) && (n+m<=numberIntervals) ) cout << setw(20) << setprecision(4) << (splines[n+m] -> getBounds()).getLowerBound() << setw(2);
// 		if( n+cut<=numberIntervals ) cout << setw(20) << setprecision(4) << (splines[n+cut] -> getBounds()).getUpperBound() << setw(2);
// 		cout << endl; 
// 		for(int m=-cut;m<=cut;m++)
// 		    if( (n+m>=0) && (n+m<=numberIntervals) )
// 			{
// 				double variable=(splines[n+m] -> getBounds()).getLowerBound(); double curError=splines[n+m] -> splineError(variable); if(curError<VERY_SMALL_NUMBER) curError=1e-6;
// 		    cout << setw(20) << setprecision(8) << ((splines[n]-> splineValue(variable)) - (splines[n+m]-> splineValue(variable)))/curError << setw(2);
// 			}
// 		if( n+cut<=numberIntervals )
// 			{
// 			double variable=(splines[n+cut] -> getBounds()).getUpperBound();  double curError=splines[n+cut] -> splineError(variable); if(curError<VERY_SMALL_NUMBER) curError=1e-6;
// 			cout << setw(20) << setprecision(8) << ((splines[n]-> splineValue(variable)) - (splines[n+cut]-> splineValue(variable)))/curError << setw(2);
// 			}
// 		cout << endl; cout << endl;
// 		
// 		start=end;
// 		}
	
	/*------------------------------------------------------------------*/
	

	splineArray result(splines);
	result.updateProperties(totalChiSquared,totalDegreesOfFreedom);
	return result;
	
}


splineArray histogramBasis::splineProcedureSmallBins(unsigned int splineOrder, unsigned int minNumberBins, double fitAcceptanceThreshold)
{
	spline spline1 = oneSplineFit(splineOrder);
	
	if(spline1.isSplineGood(fitAcceptanceThreshold))
		{
		vector<spline*> splines;
		spline * solution = new spline(spline1.getBounds(), splineOrder);
		solution -> setSpline(spline1.getCoefficients(),spline1.getErrorCoefficients(),(spline1.getChisquared())*(spline1.getDOF()),spline1.getDOF());
		splines.push_back(solution);
		splineArray result(splines);
		result.updateProperties( (spline1.getChisquared())*(spline1.getDOF()),spline1.getDOF());

		return result;
		}
		
	else 	{
		if(minNumberBins<2*splineOrder) minNumberBins=2*splineOrder;
		bool finishedUpdating=false;
		unsigned int matrixRows = basisSlots.size();
		vector<unsigned int> intervalBoundaries;
		intervalBoundaries.push_back(matrixRows/2);
	
		while(finishedUpdating==false)
			{				
			splineArray result = matchedSplineFitSmallBins(intervalBoundaries, splineOrder);	//we don't need errors, coefficients etc until the final spline -> refactor this!!!!!!!!!!!!
			finishedUpdating=true;
		
			unsigned int originalSize=intervalBoundaries.size();
			unsigned int adjustedBoundaries=0;
		
			for(unsigned int i=0;i<=originalSize;i++)
				{
				bool currentSplineGood = result.getSpline(i) -> isSplineGood(fitAcceptanceThreshold);
				if(!currentSplineGood)
					{
					int start=0; if(i>0) start = intervalBoundaries[i+adjustedBoundaries-1];
					int end=matrixRows; if(i+adjustedBoundaries<intervalBoundaries.size()) end = intervalBoundaries[i+adjustedBoundaries];
					int delta = (end-start)/2;
					if(delta>=minNumberBins)
						{
						intervalBoundaries.insert(intervalBoundaries.begin()+i+adjustedBoundaries, start+delta);
						adjustedBoundaries++;
						finishedUpdating=false;
						} 
					}
				}
				if(finishedUpdating) return result;
			}
		}
}


splineArray histogramBasis::splineProcedure(unsigned int splineOrder, unsigned int minLevel, long norm, double fitAcceptanceThreshold)
{
	if(minLevel<2) minLevel=2;
	unsigned int numberElementaryBins = basisSlots.size();
	unsigned int maxLevel=rounding(log(double(numberElementaryBins))/log(2));
	unsigned int currentLevel=minLevel;
	
	//only works for 2^n elementary bins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//make a bin hierarchy in advance by combining bins; position in vector denotes bin level: 0 is largest bin, 1 are second level bins etc
	//intervals can be marked in the same way: as a sequence 0 or 11 or 221 or 2331 etc
	vector< vector<basisSlot*> > analysisBins;
	vector<basisSlot*> currentLevelBins;
	for(unsigned int i=0;i<numberElementaryBins;i++) {basisSlots[i] -> updateEnoughSampled(); currentLevelBins.push_back(basisSlots[i] -> Clone());}
	for(unsigned int j=0;j<=maxLevel;j++)
		{
		analysisBins.insert(analysisBins.begin(),currentLevelBins);
		unsigned int currentSize = currentLevelBins.size();
		if(currentSize==1) break;
		currentLevelBins.resize(0);
		for(unsigned int i=0;i<currentSize;i+=2)
			{
			slotBounds theBounds((analysisBins[0][i] -> getBounds()).getLowerBound(), (analysisBins[0][i+1] -> getBounds()).getUpperBound());
			basisSlot * combined = new basisSlot(theBounds);
			combined -> combineWithSlot(analysisBins[0][i]);
			combined -> combineWithSlot(analysisBins[0][i+1]);
			combined -> updateEnoughSampled();
			currentLevelBins.push_back(combined);
			}
		}
		
	for(unsigned int j=0;j<=maxLevel;j++)
		for(unsigned int i=0;i<analysisBins[j].size();i++)
			{
			analysisBins[j][i] -> scale(norm);
			}
		
	vector<slotBounds> intervalBounds;
	vector<unsigned int> intervalOrders;
	vector<unsigned int> intervalNumbers;
	slotBounds histogramBounds(lowerBound,upperBound);
	intervalBounds.push_back(histogramBounds);
	unsigned int currentNumberIntervals;
	intervalOrders.push_back(0); intervalNumbers.push_back(0);
	vector< vector<basisSlot*> > currentAnalysisBins;	//the bins on each given level don't have to be in order, only analysisBins do
	for(unsigned int j=0;j<=minLevel;j++) currentAnalysisBins.push_back(analysisBins[j]);
	
/*------------ Testing alternative setup when all bins are used ------------*/
for(unsigned int j=minLevel+1;j<=maxLevel;j++) currentAnalysisBins.push_back(analysisBins[j]);
/*------------ Testing alternative setup when all bins are used ------------*/
	
	//divide bins by a level and check if intervals need to be split
	//oneSplineFit should be special case of matchedSplineFit
	while(currentLevel<maxLevel)
		{
		//check matched spline
		splineArray result = matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder);	//we don't need errors, coefficients etc until the final spline -> refactor this!!!!!!!!!!!!

		//check if all fits good: now have multiple chi^2 contributions, should redefine isSplineGood function for this
		currentNumberIntervals=intervalBounds.size();
		bool isIntervalGood[currentNumberIntervals];
		bool allSplinesGood=true;
		for(unsigned int i=0;i<currentNumberIntervals;i++)
			{
			cout << "Checking interval " << i << " (order: " << intervalOrders[i] << ", number: " << intervalNumbers[i] << ")" << endl;
			spline * currentSpline = result.getSpline(i);
			bool currentSplineGood=true;
			for(unsigned int j=0;j<=maxLevel-intervalOrders[i];j++)
				{
				//test all bins fully within spline order by order
				double currentChisq=0;
				unsigned int numberSlotsAtCurrentLevel = pow(2,j);
				for(unsigned int k=0;k<pow(2,j);k++)
					{
					basisSlot* currentSlot = analysisBins[intervalOrders[i]+j][intervalNumbers[i]*numberSlotsAtCurrentLevel+k];
					if(currentSlot -> enoughSampled())
						{
						double currentSplineIntegral = currentSpline -> splineIntegral(currentSlot -> getBounds());
						currentSplineIntegral-=currentSlot -> sampledIntegral();
						currentSplineIntegral*=1./(currentSlot -> sampledIntegralError());
						currentChisq+=currentSplineIntegral*currentSplineIntegral;
						}
					else numberSlotsAtCurrentLevel--;
					}
				if(numberSlotsAtCurrentLevel>pow(2,j)/2.)
					{
					currentChisq*=1./double(numberSlotsAtCurrentLevel);
					cout << intervalOrders[i]+j << '\t' << numberSlotsAtCurrentLevel << '\t' << currentChisq << '\t' << 1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel)) << endl;
					if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {currentSplineGood=false; allSplinesGood=false; break;}
					}
				else {cout << "Only " << numberSlotsAtCurrentLevel << " good bins, not enough for evaluation" << endl; break;}	//already not enough data in slots at this level, so can stop checking
				}
				
			isIntervalGood[i]=currentSplineGood;
			cout << "This interval fit is "; if(!currentSplineGood) cout << "not "; cout << "good" << endl;
			}
						
		unsigned int intervalCounter=0;
/*------------ Testing alternative setup when all bins are used ------------
		for(unsigned int i=0;i<currentNumberIntervals;i++)
			if(!isIntervalGood[i])
				{
				unsigned int nextOrder=intervalOrders[i]+minLevel+1;
				if(nextOrder<=maxLevel)
					{
					if(nextOrder==currentAnalysisBins.size()) {vector<basisSlot*> v; currentAnalysisBins.push_back(v);}
					for(unsigned int k=0;k<pow(2,minLevel+1);k++) currentAnalysisBins[nextOrder].push_back(analysisBins[nextOrder][pow(2,minLevel+1)*intervalNumbers[i]+k]);
					}
				}
------------ Testing alternative setup when all bins are used ------------*/
				
		for(unsigned int i=0;i<currentNumberIntervals;i++)
			{
			if(!isIntervalGood[i])
				{
				intervalOrders[intervalCounter]++;
				intervalOrders.insert(intervalOrders.begin()+intervalCounter+1,intervalOrders[intervalCounter]);
				intervalNumbers[intervalCounter]*=2;
				intervalNumbers.insert(intervalNumbers.begin()+intervalCounter+1,intervalNumbers[intervalCounter]+1);
				intervalCounter++;
				}
			intervalCounter++;
			}
			
		intervalBounds.resize(0);
		for(unsigned int i=0;i<intervalOrders.size();i++)
			{
			slotBounds currentBounds = analysisBins[intervalOrders[i]][intervalNumbers[i]] -> getBounds();
			intervalBounds.push_back(currentBounds);
			}

		if(allSplinesGood) {cout << "Good spline found!" << endl; cout << endl; break;}
		currentLevel++;
		cout << endl;
		}
		
	splineArray result = matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder);
	return result;
}


splineArray matchedSplineFit(vector< vector<basisSlot*> > currentAnalysisBins, vector< slotBounds > intervalBounds, unsigned int splineOrder)
{
	//test weighting bins at different orders differently----------------------------
	double orderWeighting[currentAnalysisBins.size()];
	for(unsigned int i=0;i<currentAnalysisBins.size();i++) orderWeighting[i]=1.;
	//-------------------------------------------------------------------------------
	
	unsigned int numberIntervals=intervalBounds.size()-1;
	unsigned int matrixRows = 0;
	for(unsigned int i=0;i<currentAnalysisBins.size();i++) matrixRows+=currentAnalysisBins[i].size();
	unsigned int matrixCols=numberIntervals+splineOrder;
	
	//this vector (size = number intervals) stores for each interval the bins that are fully inside the interval -> for chi^2 on interval evaluation
	vector< vector<unsigned int> > binsFullyInsideInterval(numberIntervals+1);
	
	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows, matrixCols);
	double* b = integralVector(currentAnalysisBins);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows);
	for(unsigned int i=0;i<matrixRows;i++) gsl_vector_set(gslb,i,b[i]);

	double* binVec=binomialVector(splineOrder);
	double intervalMatrix[numberIntervals][splineOrder-1];
	double lastIntervalMatrix[numberIntervals][splineOrder-1];
	
	double* basicDesignMatrix = designMatrix(currentAnalysisBins, splineOrder);
	double currentLowerBound=0; double currentUpperBound; vector<basisSlot*>:: iterator it;
	for(unsigned int i=0;i<numberIntervals;i++)
		{
		currentUpperBound = intervalBounds[i].getUpperBound();
		
		for(unsigned int k=0; k < splineOrder-1;k++)
			{
			intervalMatrix[i][k]=(pow(currentUpperBound,splineOrder-1-k)-pow(currentLowerBound,splineOrder-1-k))*binVec[k]*pow(-1,splineOrder-2-k);
			lastIntervalMatrix[i][k]=-pow(currentUpperBound,splineOrder-1-k)*binVec[k]*pow(-1,splineOrder-2-k);
			}
		
		currentLowerBound=currentUpperBound;
		}

	//go through all bins after intervals
	//when bins are bigger than intervals their boundaries always coincide with some interval boundaries, because splitting in the same manner
	//level of slot denotes how many slots we actually use of the given level (for now) -> might not need this anymore
	unsigned int currentMatrixRow=0;
	for(unsigned int j1=0;j1 < currentAnalysisBins.size();j1++)
		for(unsigned int j2=0;j2 < currentAnalysisBins[j1].size();j2++)
			{
			slotBounds currentBounds=(currentAnalysisBins[j1][j2] -> getBounds());
			double currentError=currentAnalysisBins[j1][j2] -> sampledIntegralError();
			unsigned int start; unsigned int end;
			for(unsigned int i=0;i<=numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {start=i; break;}
				}
			for(unsigned int i=start;i<=numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {end=i;}
				}
		
			//regular elements for a0, a1, a2
			for(unsigned int k=0; k < splineOrder-1;k++)
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,k,basicDesignMatrix[currentMatrixRow*splineOrder+k]/orderWeighting[j1]);
		
			//very last element (x_i+1^4-x_i^4)/4
			if(start==end)
				{
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,start+splineOrder-1,basicDesignMatrix[currentMatrixRow*splineOrder+splineOrder-1]/orderWeighting[j1]);
				binsFullyInsideInterval[start].push_back(currentMatrixRow);
				}
			else
				{
				for(unsigned int i=start;i<=end;i++)
					{
					if( currentError < VERY_SMALL_NUMBER) gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1,0);
					else gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1, splineBasisFunction(intervalBounds[i],splineOrder-1)/currentError/sqrt(double(currentAnalysisBins[j1].size()))/orderWeighting[j1]);
					}
				}
		
			//other a3 contributions
			if((start==end)&&(start>0))
				{
				for(unsigned int kk=0; kk < splineOrder-1;kk++)
					{
					for(unsigned int k=0; k < start; k++)
						gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]/orderWeighting[j1]);
				
					gsl_matrix_set(fullDesignMatrix, currentMatrixRow, start+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,start+splineOrder-1)+lastIntervalMatrix[start-1][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]/orderWeighting[j1]);
					}
				}
			else if(start!=end)
				{
				if(start==0) start++;
				for(unsigned int i=start;i<=end;i++)
					{
					if( currentError > VERY_SMALL_NUMBER)
						{
						for(unsigned int kk=0; kk < splineOrder-1;kk++)
							{
							for(unsigned int k=0; k < i; k++)
								gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/sqrt(double(currentAnalysisBins[j1].size()))/orderWeighting[j1]);
						
							gsl_matrix_set(fullDesignMatrix, currentMatrixRow, i+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,i+splineOrder-1)+lastIntervalMatrix[i-1][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/sqrt(double(currentAnalysisBins[j1].size()))/orderWeighting[j1]);
							}
						}
					}
				}
				
			gsl_vector_set(gslb,currentMatrixRow,gsl_vector_get(gslb,currentMatrixRow)/orderWeighting[j1]);
			currentMatrixRow++;
			}
	
	//cout << "design matrix" << endl;
	//for(unsigned int i=0;i<matrixRows;i++) {for(unsigned int j=0;j<matrixCols;j++) cout << gsl_matrix_get(fullDesignMatrix,i,j) << '\t'; cout << endl;}
	//cout << endl;
	
	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
	
	double totalChiSquared=solveSVD(fullDesignMatrix, V, gslb, diagonal, x);	//these need to be adjusted
	int totalDegreesOfFreedom=matrixRows-(numberIntervals+splineOrder-1);
	vector<spline*> splines;

	//this can actually be done in splineProcedure and probably should
	currentMatrixRow=0;
	cout << "Checking separate chi_n^2/n in spline fit" << endl;
	cout << "bin order" << '\t' << "n" << '\t' << "chi_n^2/n" << '\t' << "sqrt(2/n)" << endl;
	for(unsigned int j=0;j < currentAnalysisBins.size();j++)
		{
		int dof = currentAnalysisBins[j].size();
		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, currentMatrixRow, dof);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		cout << j << '\t' << dof << '\t' << chisq << '\t' << sqrt(2./double(dof)) << endl;
		currentMatrixRow+=dof;
		}
	
	double diagScaling[matrixCols];
	for(unsigned int j=0;j < matrixCols;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./gsl_vector_get(diagonal,j);
		}
	
	gsl_matrix * Vscaled = gsl_matrix_alloc (matrixCols, matrixCols);
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			gsl_matrix_set(Vscaled,i,j,gsl_matrix_get(V,i,j)*diagScaling[j]);
		
	gsl_matrix * fullCovarianceMatrix = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1, Vscaled, 0, fullCovarianceMatrix);
	double covarianceMatrix[matrixCols][matrixCols];
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			covarianceMatrix[i][j]= gsl_matrix_get (fullCovarianceMatrix, i,j);
		
	vector<double> theCoefficients;
	for(unsigned int n=0;n<=numberIntervals;n++)
		{
		if(n==0) 
			{
			for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,i));
			}
		else
			{
			for(unsigned int i=0;i<splineOrder-1;i++) theCoefficients[i]+=(gsl_vector_get(x,n+splineOrder-1)-theCoefficients[splineOrder-1])*lastIntervalMatrix[n-1][i];
			theCoefficients[splineOrder-1]=gsl_vector_get(x,n+splineOrder-1);
			}
		
		vector<double> theErrorCoefficients;
		for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0); 
		
		if(n>0)
			{
			for(unsigned int j=0;j < splineOrder-1;j++)
				{
				for(unsigned int k=j;k < splineOrder-1;k++)
					{
					covarianceMatrix[j][k]+=(covarianceMatrix[j][splineOrder+n-1]-covarianceMatrix[j][splineOrder-1])*lastIntervalMatrix[n-1][k];
					covarianceMatrix[j][k]+=(covarianceMatrix[k][splineOrder+n-1]-covarianceMatrix[k][splineOrder-1])*lastIntervalMatrix[n-1][j];
					covarianceMatrix[j][k]+=(covarianceMatrix[splineOrder-1][splineOrder-1] + covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - 2*covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][k]*lastIntervalMatrix[n-1][j];
					}
				covarianceMatrix[j][splineOrder-1]=covarianceMatrix[j][splineOrder+n-1];
				covarianceMatrix[j][splineOrder-1]+=(covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][j];
				
				for(unsigned int k=n;k < matrixCols-splineOrder;k++)
					covarianceMatrix[j][k+splineOrder]+=(covarianceMatrix[splineOrder+n-1][k+splineOrder] - covarianceMatrix[splineOrder+n-2][k+splineOrder])*lastIntervalMatrix[n-1][j];
				}
			covarianceMatrix[splineOrder-1][splineOrder-1]=covarianceMatrix[n+splineOrder-1][n+splineOrder-1];
			}
		
		for(unsigned int j=0;j < splineOrder;j++)
			for(unsigned int k=j;k < splineOrder;k++)
				{
				theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				if(j!=k) theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				}
			
		gsl_vector * currentgslb = gsl_vector_alloc (binsFullyInsideInterval[n].size());
		for(unsigned int j=0;j < binsFullyInsideInterval[n].size();j++) gsl_vector_set(currentgslb,j, gsl_vector_get(gslb,binsFullyInsideInterval[n][j]));
		double chisq; gsl_blas_ddot (currentgslb, currentgslb, &chisq);
		int dof = binsFullyInsideInterval[n].size();
		
		spline * solution = new spline(intervalBounds[n], splineOrder); solution -> setSpline(theCoefficients,theErrorCoefficients,chisq,dof);
		splines.push_back(solution);
		}
	
	splineArray result(splines);
	result.updateProperties(totalChiSquared,totalDegreesOfFreedom);
	return result;
	
}

splineArray histogramBasis::splineProcedureOverlap(unsigned int splineOrder, long norm, double fitAcceptanceThreshold)
{
	unsigned int numberElementaryBins = basisSlots.size();
	unsigned int maxLevel=rounding(log(double(numberElementaryBins))/log(2));
	unsigned int currentLevel=2;
	
	//only works for 2^n elementary bins!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//make a bin hierarchy in advance by combining bins; position in vector denotes bin level: 0 is largest bin, 1 are second level bins etc
	//intervals can be marked in the same way: as a sequence 0 or 11 or 221 or 2331 etc
	vector< vector<basisSlot*> > analysisBins;
	vector<basisSlot*> currentLevelBins;
	for(unsigned int i=0;i<numberElementaryBins;i++) {basisSlots[i] -> updateEnoughSampled(); currentLevelBins.push_back(basisSlots[i] -> Clone());}
	for(unsigned int j=0;j<=maxLevel;j++)
		{
		analysisBins.insert(analysisBins.begin(),currentLevelBins);
		unsigned int currentSize = currentLevelBins.size();
		if(currentSize==1) break;
		currentLevelBins.resize(0);
		unsigned int next=2;
		if(j==0) next=1;
		for(unsigned int i=0;i<currentSize-next;i+=next)
			{
			slotBounds theBounds((analysisBins[0][i] -> getBounds()).getLowerBound(), (analysisBins[0][i+next] -> getBounds()).getUpperBound());
			basisSlot * combined = new basisSlot(theBounds);
			combined -> combineWithSlot(analysisBins[0][i]);
			combined -> combineWithSlot(analysisBins[0][i+next]);
			combined -> updateEnoughSampled();
			currentLevelBins.push_back(combined);
			}
		}
	
	for(unsigned int j=0;j<=maxLevel;j++)
		for(unsigned int i=0;i<analysisBins[j].size();i++)
			{
			analysisBins[j][i] -> scale(norm);
			}
		
	vector<slotBounds> intervalBounds;
	vector<unsigned int> intervalOrders;
	vector<unsigned int> intervalNumbers;
	slotBounds histogramBounds(lowerBound,upperBound);
	intervalBounds.push_back(histogramBounds);
	unsigned int currentNumberIntervals;
	intervalOrders.push_back(0); intervalNumbers.push_back(0);
	vector< vector<basisSlot*> > currentAnalysisBins;	//the bins on each given level don't have to be in order, only analysisBins do
	for(unsigned int j=0;j<=maxLevel;j++) currentAnalysisBins.push_back(analysisBins[j]);
	
	//divide bins by a level and check if intervals need to be split
	while(currentLevel<maxLevel)
		{
		//check matched spline
		splineArray result = matchedSplineFitOverlap(currentAnalysisBins, intervalBounds, splineOrder);	//we don't need errors, coefficients etc until the final spline -> refactor this!!!!!!!!!!!!
		
		//check if all fits good: now have multiple chi^2 contributions, should redefine isSplineGood function for this
		currentNumberIntervals=intervalBounds.size();
		bool isIntervalGood[currentNumberIntervals];
		bool allSplinesGood=true;
		for(unsigned int i=0;i<currentNumberIntervals;i++)
			{
			cout << "Checking interval " << i << " (order: " << intervalOrders[i] << ", number: " << intervalNumbers[i] << ")" << endl;
			spline * currentSpline = result.getSpline(i);
			bool currentSplineGood=true;
			for(unsigned int j=0;j<=maxLevel-intervalOrders[i];j++)
				{
				//test all bins fully within spline order by order
				double currentChisq=0;
				unsigned int numberSlotsAtCurrentLevel=analysisBins[intervalOrders[i]+j].size();
				for(unsigned int k=0;k<analysisBins[intervalOrders[i]+j].size();k++)
					{
					basisSlot* currentSlot = analysisBins[intervalOrders[i]+j][k];
					if( (currentSlot -> enoughSampled()) && (currentSlot->getBounds().getLowerBound()>=intervalBounds[i].getLowerBound()) && (currentSlot->getBounds().getUpperBound()<=intervalBounds[i].getUpperBound()) )
						{
						double currentSplineIntegral = currentSpline -> splineIntegral(currentSlot -> getBounds());
						currentSplineIntegral-=currentSlot -> sampledIntegral();
						currentSplineIntegral*=1./(currentSlot -> sampledIntegralError());
						currentChisq+=currentSplineIntegral*currentSplineIntegral;
						}
					else numberSlotsAtCurrentLevel--;
					}
				if(numberSlotsAtCurrentLevel>pow(2,j)/2.)
					{
					currentChisq*=1./double(numberSlotsAtCurrentLevel);
					cout << intervalOrders[i]+j << '\t' << numberSlotsAtCurrentLevel << '\t' << currentChisq << '\t' << 1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel)) << endl;
					if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {currentSplineGood=false; allSplinesGood=false; break;}
					}
				else {cout << "Only " << numberSlotsAtCurrentLevel << " good bins, not enough for evaluation" << endl; break;}	//already not enough data in slots at this level, so can stop checking
				}
			
			isIntervalGood[i]=currentSplineGood;
			cout << "This interval fit is "; if(!currentSplineGood) cout << "not "; cout << "good" << endl;
			}
		
		unsigned int intervalCounter=0;
		for(unsigned int i=0;i<currentNumberIntervals;i++)
			{
			if(!isIntervalGood[i])
				{
				intervalOrders[intervalCounter]++;
				intervalOrders.insert(intervalOrders.begin()+intervalCounter+1,intervalOrders[intervalCounter]);
				intervalNumbers[intervalCounter]*=2;
				intervalNumbers.insert(intervalNumbers.begin()+intervalCounter+1,intervalNumbers[intervalCounter]+1);
				intervalCounter++;
				}
			intervalCounter++;
			}
		
		intervalBounds.resize(0);
		for(unsigned int i=0;i<intervalOrders.size();i++)
			{
			slotBounds currentBounds = analysisBins[intervalOrders[i]][2*intervalNumbers[i]] -> getBounds();
			intervalBounds.push_back(currentBounds);
			}
		
		if(allSplinesGood) {cout << "Good spline found!" << endl; cout << endl; break;}
		currentLevel++;
		cout << endl;
		}
	
	splineArray result = matchedSplineFitOverlap(currentAnalysisBins, intervalBounds, splineOrder);
	return result;
}

splineArray matchedSplineFitOverlap(vector< vector<basisSlot*> > currentAnalysisBins, vector< slotBounds > intervalBounds, unsigned int splineOrder)
{
	unsigned int numberIntervals=intervalBounds.size()-1;
	unsigned int matrixRows = 0;
	for(unsigned int i=0;i<currentAnalysisBins.size();i++) matrixRows+=currentAnalysisBins[i].size();
	unsigned int matrixCols=numberIntervals+splineOrder;
	
	//this vector (size = number intervals) stores for each interval the bins that are fully inside the interval -> for chi^2 on interval evaluation
	vector< vector<unsigned int> > binsFullyInsideInterval(numberIntervals+1);
	
	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows, matrixCols);
	double* b = integralVector(currentAnalysisBins);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows);
	for(unsigned int i=0;i<matrixRows;i++) gsl_vector_set(gslb,i,b[i]);
	
	double* binVec=binomialVector(splineOrder);
	double intervalMatrix[numberIntervals][splineOrder-1];
	double lastIntervalMatrix[numberIntervals][splineOrder-1];
	
	double* basicDesignMatrix = designMatrix(currentAnalysisBins, splineOrder);
	double currentLowerBound=0; double currentUpperBound; vector<basisSlot*>:: iterator it;
	for(unsigned int i=0;i<numberIntervals;i++)
		{
		currentUpperBound = intervalBounds[i].getUpperBound();
		
		for(unsigned int k=0; k < splineOrder-1;k++)
			{
			intervalMatrix[i][k]=(pow(currentUpperBound,splineOrder-1-k)-pow(currentLowerBound,splineOrder-1-k))*binVec[k]*pow(-1,splineOrder-2-k);
			lastIntervalMatrix[i][k]=-pow(currentUpperBound,splineOrder-1-k)*binVec[k]*pow(-1,splineOrder-2-k);
			}
		
		currentLowerBound=currentUpperBound;
		}
	
	//go through all bins after intervals
	unsigned int currentMatrixRow=0;
	for(unsigned int j1=0;j1 < currentAnalysisBins.size();j1++)
		for(unsigned int j2=0;j2 < currentAnalysisBins[j1].size();j2++)
			{
			slotBounds currentBounds=(currentAnalysisBins[j1][j2] -> getBounds());
			double currentError=currentAnalysisBins[j1][j2] -> sampledIntegralError();
			unsigned int start; unsigned int end;
			for(unsigned int i=0;i<=numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {start=i; break;}
				}
			for(unsigned int i=start;i<=numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {end=i;}
				}
			
			slotBounds startBounds(max(intervalBounds[start].getLowerBound(), currentBounds.getLowerBound()), intervalBounds[start].getUpperBound());
			slotBounds endBounds(intervalBounds[end].getLowerBound(), min(intervalBounds[end].getUpperBound(), currentBounds.getUpperBound()));
			
			//regular elements for a0, a1, a2
			for(unsigned int k=0; k < splineOrder-1;k++)
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,k,basicDesignMatrix[currentMatrixRow*splineOrder+k]);
			
			//very last element (x_i+1^4-x_i^4)/4
			if(start==end)
				{
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,start+splineOrder-1,basicDesignMatrix[currentMatrixRow*splineOrder+splineOrder-1]);
				binsFullyInsideInterval[start].push_back(currentMatrixRow);
				}
			else
				{
				for(unsigned int i=start;i<=end;i++)
					{
					slotBounds adjustedBounds=intervalBounds[i];
					if(i==start) adjustedBounds=startBounds;
					else if(i==end) adjustedBounds=endBounds;
						
					if( currentError < VERY_SMALL_NUMBER) gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1,0);
					else gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1, splineBasisFunction(adjustedBounds,splineOrder-1)/currentError/sqrt(double(currentAnalysisBins[j1].size())));
					}
				}
			
			//other a3 contributions
			if((start==end)&&(start>0))
				{
				for(unsigned int kk=0; kk < splineOrder-1;kk++)
					{
					for(unsigned int k=0; k < start; k++)
						gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]);
					
					gsl_matrix_set(fullDesignMatrix, currentMatrixRow, start+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,start+splineOrder-1)+lastIntervalMatrix[start-1][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]);
					}
				}
			else if(start!=end)
				{
				bool adjust=true;
				if(start==0) {start++; adjust=false;}
				for(unsigned int i=start;i<=end;i++)
					{
					if( currentError > VERY_SMALL_NUMBER)
						{
						slotBounds adjustedBounds=intervalBounds[i];
						if( (i==start)  && adjust ) adjustedBounds=startBounds;
						else if(i==end) adjustedBounds=endBounds;

						for(unsigned int kk=0; kk < splineOrder-1;kk++)
							{
							for(unsigned int k=0; k < i; k++)
								gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*splineBasisFunction(adjustedBounds,kk)/currentError/sqrt(double(currentAnalysisBins[j1].size())));
							
							gsl_matrix_set(fullDesignMatrix, currentMatrixRow, i+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,i+splineOrder-1)+lastIntervalMatrix[i-1][kk]*splineBasisFunction(adjustedBounds,kk)/currentError/sqrt(double(currentAnalysisBins[j1].size())));
							}
						}
					}
				}
			
			gsl_vector_set(gslb,currentMatrixRow,gsl_vector_get(gslb,currentMatrixRow));
			currentMatrixRow++;
			}
		
		//cout << "design matrix" << endl;
		//for(unsigned int i=0;i<matrixRows;i++) {for(unsigned int j=0;j<matrixCols;j++) cout << gsl_matrix_get(fullDesignMatrix,i,j) << '\t'; cout << endl;}
		//cout << endl;
		
	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
		
	double totalChiSquared=solveSVD(fullDesignMatrix, V, gslb, diagonal, x);	//these need to be adjusted
	int totalDegreesOfFreedom=matrixRows-(numberIntervals+splineOrder-1);
	vector<spline*> splines;
		
	//this can actually be done in splineProcedure and probably should
	currentMatrixRow=0;
	cout << "Checking separate chi_n^2/n in spline fit" << endl;
	cout << "bin order" << '\t' << "n" << '\t' << "chi_n^2/n" << '\t' << "sqrt(2/n)" << endl;
	for(unsigned int j=0;j < currentAnalysisBins.size();j++)
		{
		int dof = currentAnalysisBins[j].size();
		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, currentMatrixRow, dof);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		cout << j << '\t' << dof << '\t' << chisq << '\t' << sqrt(2./double(dof)) << endl;
		currentMatrixRow+=dof;
		}
		
	double diagScaling[matrixCols];
	for(unsigned int j=0;j < matrixCols;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./gsl_vector_get(diagonal,j);
		}
		
	gsl_matrix * Vscaled = gsl_matrix_alloc (matrixCols, matrixCols);
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			gsl_matrix_set(Vscaled,i,j,gsl_matrix_get(V,i,j)*diagScaling[j]);
			
	gsl_matrix * fullCovarianceMatrix = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1, Vscaled, 0, fullCovarianceMatrix);
	double covarianceMatrix[matrixCols][matrixCols];
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			covarianceMatrix[i][j]= gsl_matrix_get (fullCovarianceMatrix, i,j);
			
	vector<double> theCoefficients;
	for(unsigned int n=0;n<=numberIntervals;n++)
		{
		if(n==0) 
			{
			for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,i));
			}
		else
			{
			for(unsigned int i=0;i<splineOrder-1;i++) theCoefficients[i]+=(gsl_vector_get(x,n+splineOrder-1)-theCoefficients[splineOrder-1])*lastIntervalMatrix[n-1][i];
			theCoefficients[splineOrder-1]=gsl_vector_get(x,n+splineOrder-1);
			}
			
		vector<double> theErrorCoefficients;
		for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0); 
			
		if(n>0)
			{
			for(unsigned int j=0;j < splineOrder-1;j++)
				{
				for(unsigned int k=j;k < splineOrder-1;k++)
					{
					covarianceMatrix[j][k]+=(covarianceMatrix[j][splineOrder+n-1]-covarianceMatrix[j][splineOrder-1])*lastIntervalMatrix[n-1][k];
					covarianceMatrix[j][k]+=(covarianceMatrix[k][splineOrder+n-1]-covarianceMatrix[k][splineOrder-1])*lastIntervalMatrix[n-1][j];
					covarianceMatrix[j][k]+=(covarianceMatrix[splineOrder-1][splineOrder-1] + covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - 2*covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][k]*lastIntervalMatrix[n-1][j];
					}
				covarianceMatrix[j][splineOrder-1]=covarianceMatrix[j][splineOrder+n-1];
				covarianceMatrix[j][splineOrder-1]+=(covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][j];
					
				for(unsigned int k=n;k < matrixCols-splineOrder;k++)
					covarianceMatrix[j][k+splineOrder]+=(covarianceMatrix[splineOrder+n-1][k+splineOrder] - covarianceMatrix[splineOrder+n-2][k+splineOrder])*lastIntervalMatrix[n-1][j];
				}
			covarianceMatrix[splineOrder-1][splineOrder-1]=covarianceMatrix[n+splineOrder-1][n+splineOrder-1];
			}
			
		for(unsigned int j=0;j < splineOrder;j++)
			for(unsigned int k=j;k < splineOrder;k++)
				{
				theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				if(j!=k) theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				}
				
		gsl_vector * currentgslb = gsl_vector_alloc (binsFullyInsideInterval[n].size());
		for(unsigned int j=0;j < binsFullyInsideInterval[n].size();j++) gsl_vector_set(currentgslb,j, gsl_vector_get(gslb,binsFullyInsideInterval[n][j]));
		double chisq; gsl_blas_ddot (currentgslb, currentgslb, &chisq);
		int dof = binsFullyInsideInterval[n].size();
			
		spline * solution = new spline(intervalBounds[n], splineOrder); solution -> setSpline(theCoefficients,theErrorCoefficients,chisq,dof);
		splines.push_back(solution);
		}
		
	splineArray result(splines);
	result.updateProperties(totalChiSquared,totalDegreesOfFreedom);
	return result;
		
}



