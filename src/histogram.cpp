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
		basisSlot* currentSlot = basisSlots[i] -> Clone();
		currentSlot -> scale(norm);
		scaledSlots.push_back(currentSlot);
		}
	
	histogramBasis result(scaledSlots);
	return result;

}

bool histogramBasis::addAnotherHistogram(histogramBasis anotherHistogram)
{
	bool addable=false;
	if(basisSlots.size()==anotherHistogram.getSize())
		{
		addable=true;
		for(unsigned int i=0;i<basisSlots.size();i++)
			{
			addable = basisSlots[i] -> addAnotherSlot(anotherHistogram.getSlot(i));
			}
		}
	return addable;
}


vector< vector< basisSlot* > > histogramBasis::binHierarchy(long norm)
{
	unsigned int numberElementaryBins = basisSlots.size();
	unsigned int maxLevel=rounding(log(double(numberElementaryBins))/log(2));
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
	
	for(unsigned int j=0;j<analysisBins.size();j++)
		{
		int initialAnalysisBinsLevelSize=analysisBins[j].size()-1;
		for(int i=initialAnalysisBinsLevelSize;i>=0;i--)
			{
			if(analysisBins[j][i] -> enoughSampled()) analysisBins[j][i] -> scale(norm);
			else analysisBins[j].erase(analysisBins[j].begin()+i);
			}
		}
	
	return analysisBins;
}


splineArray histogramBasis::splineProcedure(unsigned int splineOrder, unsigned int minLevel, long norm, double fitAcceptanceThreshold, double gluingFactor)
{
	vector<double> a3array; vector<double> chisqArray;
	
	if(minLevel<2) minLevel=2;
	unsigned int numberElementaryBins = basisSlots.size();
	unsigned int maxLevel=rounding(log(double(numberElementaryBins))/log(2));
	unsigned int currentLevel=minLevel;
	
	//only works for 2^n elementary bins!
	//make a bin hierarchy in advance by combining bins; position in vector denotes bin level: 0 is largest bin, 1 are second level bins etc
	//intervals can be marked in the same way: as a sequence 0 or 11 or 221 or 2331 etc
	vector< vector<basisSlot*> > analysisBins=binHierarchy(norm);

	if(isDataConsistentWithZero(analysisBins)==true) cout << "Data is consistent with zero on the interval" << endl;
		
	vector<slotBounds> intervalBounds;
	vector<unsigned int> intervalOrders;
	vector<unsigned int> intervalNumbers;
	slotBounds histogramBounds(lowerBound,upperBound);
	intervalBounds.push_back(histogramBounds);
	unsigned int currentNumberIntervals;
	intervalOrders.push_back(0); intervalNumbers.push_back(0);
	vector< vector<basisSlot*> > currentAnalysisBins;					//the bins on each given level don't have to be in order, only analysisBins do
	for(unsigned int j=0;j<=maxLevel;j++) currentAnalysisBins.push_back(analysisBins[j]);	//all bins are used

	//divide bins by a level and check if intervals need to be split
	bool allSplinesGood;
	while(currentLevel<maxLevel)
		{
		//check matched spline
		splineArray result = matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder, 0, a3array, chisqArray);	//don't need errors, coefficients etc until the final spline -> refactor this!!!!!!!!

		if(result.checkOverallAcceptance(fitAcceptanceThreshold)==true) allSplinesGood=true;
		else {
			//check if all fits good: now have multiple chi^2 contributions, should redefine isSplineGood function for this
			currentNumberIntervals=intervalBounds.size();
			bool isIntervalGood[currentNumberIntervals];
			allSplinesGood=true;
			chisqArray.resize(0);
			for(unsigned int i=0;i<currentNumberIntervals;i++)
				{
				cout << "Checking interval " << i << " (order: " << intervalOrders[i] << ", number: " << intervalNumbers[i] << ")" << endl;
				spline * currentSpline = result.getSpline(i);
				bool currentSplineGood=true;
				chisqArray.push_back(1+fitAcceptanceThreshold*sqrt(2.));
				for(unsigned int j=0;j<=maxLevel-intervalOrders[i];j++)
					{
					//test all bins fully within spline order by order
					double currentChisq=0;
					unsigned int numberSlotsAtCurrentLevel=0;
					for(unsigned int k=0;k<analysisBins[j+intervalOrders[i]].size();k++)
						{
						basisSlot* currentSlot = analysisBins[j+intervalOrders[i]][k];
						if((currentSlot -> getBounds()).overlapping(currentSpline -> getBounds()))
							{
							numberSlotsAtCurrentLevel++;
							double currentSplineIntegral = currentSpline -> splineIntegral(currentSlot -> getBounds());
							currentSplineIntegral-=currentSlot -> sampledIntegral();
							currentSplineIntegral*=1./(currentSlot -> sampledIntegralError());
							currentChisq+=currentSplineIntegral*currentSplineIntegral;
							}
						}
					
					if(numberSlotsAtCurrentLevel>pow(2,j)/2.)
						{
						currentChisq*=1./double(numberSlotsAtCurrentLevel);
						double delta=1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))-currentChisq; if(delta<0) delta=0;
						if(delta<chisqArray[i]) chisqArray[i]=delta;
						cout << intervalOrders[i]+j << '\t' << numberSlotsAtCurrentLevel << '\t' << currentChisq << '\t' << 1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel)) << endl;
						if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {currentSplineGood=false; allSplinesGood=false; break;}
						}
					else {cout << "Only " << numberSlotsAtCurrentLevel << " good bins, not enough for evaluation" << endl; break;}
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
				slotBounds currentBounds( (basisSlots[pow(2,maxLevel-intervalOrders[i])*intervalNumbers[i]] -> getBounds()).getLowerBound(), (basisSlots[pow(2,maxLevel-intervalOrders[i])*(intervalNumbers[i]+1)-1] -> getBounds()).getUpperBound());
			
				intervalBounds.push_back(currentBounds);
				}
		}
		
		if(allSplinesGood) {cout << "Good spline found!" << endl; cout << endl; break;}
		currentLevel++;
		cout << endl;
		}
	
	if(!allSplinesGood) gluingFactor=0;
	/**/if(gluingFactor>0)
	{
		a3array.resize(0);
		splineArray geta3=matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder, 0, a3array, chisqArray);
		for(unsigned int i=0;i<intervalOrders.size();i++) a3array.push_back((geta3.getSpline(i) -> getCoefficients())[splineOrder-1]);
	}
	
	splineArray result = matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder, gluingFactor, a3array, chisqArray);
	
	/**/if(gluingFactor>0)
	{
		cout << endl; cout << "Iterative consistent constraints procedure" << endl; cout << endl;
		for(int repeat=0;repeat<1;repeat++)
		{
			bool CCprocedureConverged=false; int iterations=0; double maximalGoodGluingFactor=0; double minimalBadGluingFactor=1/VERY_SMALL_NUMBER; int maxIterations=20;
			while(CCprocedureConverged==false)
			{
				iterations++;
				allSplinesGood=true;
				for(unsigned int i=0;i<intervalBounds.size();i++)
				{
					//cout << "Final check after CC; interval " << i << " (order: " << intervalOrders[i] << ", number: " << intervalNumbers[i] << ")" << endl;
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
							//cout << intervalOrders[i]+j << '\t' << numberSlotsAtCurrentLevel << '\t' << currentChisq << '\t' << 1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel)) << endl;
							if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {currentSplineGood=false; allSplinesGood=false; break;}
						}
						else {/*cout << "Only " << numberSlotsAtCurrentLevel << " good bins, not enough for evaluation" << endl;*/ break;}	//already not enough data in slots at this level, so can stop checking
					}
					
					//cout << "This interval fit is "; if(!currentSplineGood) cout << "not "; cout << "good" << endl;
				}
				cout << "Gluing factor " << gluingFactor << "; the spline is "; if(!allSplinesGood) cout << "not "; cout << "good" << endl;
				
				if(!allSplinesGood)
				{
					if(gluingFactor<minimalBadGluingFactor) minimalBadGluingFactor=gluingFactor;
					
					if(iterations>maxIterations) {gluingFactor=maximalGoodGluingFactor; CCprocedureConverged=true; allSplinesGood=true;}
					else if(maximalGoodGluingFactor>0) {gluingFactor=(maximalGoodGluingFactor+minimalBadGluingFactor)/2.;}
					else gluingFactor*=0.1;
				}
				else	{
					if(gluingFactor>maximalGoodGluingFactor) maximalGoodGluingFactor=gluingFactor;
					
					if(iterations>maxIterations) {CCprocedureConverged=true; gluingFactor=maximalGoodGluingFactor;}
					else if(minimalBadGluingFactor<1/VERY_SMALL_NUMBER) {gluingFactor=(maximalGoodGluingFactor+minimalBadGluingFactor)/2.;}
					else gluingFactor*=2;
				}
				
				result = matchedSplineFit(currentAnalysisBins, intervalBounds, splineOrder, gluingFactor, a3array, chisqArray);
			}
			/**/a3array.resize(0); chisqArray.resize(0);
			for(unsigned int i=0;i<intervalOrders.size();i++) a3array.push_back((result.getSpline(i) -> getCoefficients())[splineOrder-1]);
			for(unsigned int i=0;i<intervalBounds.size();i++)
			{
				spline * currentSpline = result.getSpline(i);
				chisqArray.push_back(1+fitAcceptanceThreshold*sqrt(2.));
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
						double delta=1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))-currentChisq; if(delta<0) delta=0;
						if(delta<chisqArray[i]) chisqArray[i]=delta;
						if(currentChisq>1+fitAcceptanceThreshold*sqrt(2./double(numberSlotsAtCurrentLevel))) {break;}
					}
					else break;
				}
			}/**/
			gluingFactor=1;
		}
	}
	/**/
	
	result.updateGoodness(allSplinesGood);
	return result;
}




splineArray matchedSplineFit(vector< vector<basisSlot*> > currentAnalysisBins, vector< slotBounds > intervalBounds, unsigned int splineOrder, double gluingFactor, vector<double> a3array, vector<double> chisqArray)
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
	
/**/ unsigned int matrixRows2=matrixRows; if(gluingFactor>0) {matrixRows2+=numberIntervals+1; if(numberIntervals>0) matrixRows2++;}
	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows2, matrixCols);
	double* b = integralVector(currentAnalysisBins);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows2);
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
			//double weightFactor=sqrt(double(currentAnalysisBins[j1].size()));
			double weightFactor=sqrt(double(pow(2,j1)));
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
					else gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1, splineBasisFunction(intervalBounds[i],splineOrder-1)/currentError/weightFactor/orderWeighting[j1]);
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
								gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/weightFactor/orderWeighting[j1]);
						
							gsl_matrix_set(fullDesignMatrix, currentMatrixRow, i+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,i+splineOrder-1)+lastIntervalMatrix[i-1][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/weightFactor/orderWeighting[j1]);
							}
						}
					}
				}
				
			gsl_vector_set(gslb,currentMatrixRow,gsl_vector_get(gslb,currentMatrixRow)/orderWeighting[j1]);
			currentMatrixRow++;
			}
			
/**/ if(gluingFactor>0)
	{
	for(unsigned int i=0;i<numberIntervals;i++) 
		{
		double minChisq=min(chisqArray[i+1],chisqArray[i]); if(i>0) minChisq=min(minChisq,chisqArray[i-1]); if(i<numberIntervals-1) minChisq=min(minChisq,chisqArray[i+2]);
		double setToWhat=sqrt(gluingFactor/double(numberIntervals))*minChisq/(a3array[i+1]-a3array[i]);
		gsl_matrix_set(fullDesignMatrix, matrixRows+i, splineOrder-1+i, -setToWhat);
		gsl_matrix_set(fullDesignMatrix, matrixRows+i, splineOrder+i, setToWhat);
		}
	gsl_matrix_set(fullDesignMatrix, matrixRows+numberIntervals, splineOrder-1, sqrt(gluingFactor/double(numberIntervals))*chisqArray[0]/(a3array[0]));
	if(numberIntervals>0) gsl_matrix_set(fullDesignMatrix, matrixRows+numberIntervals+1, splineOrder-1+numberIntervals, sqrt(gluingFactor/double(numberIntervals))*chisqArray[numberIntervals]/(a3array[numberIntervals]));
	}
/**/
	
	//cout << "design matrix" << endl;
	//for(unsigned int i=0;i<matrixRows2;i++) {for(unsigned int j=0;j<matrixCols;j++) cout << gsl_matrix_get(fullDesignMatrix,i,j) << '\t'; cout << endl;}
	//cout << endl;
	
	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
	
	double totalChiSquared=solveSVD(fullDesignMatrix, V, gslb, diagonal, x);	//these need to be adjusted
	int totalDegreesOfFreedom=matrixRows-(numberIntervals+splineOrder-1);
	vector<spline*> splines;

	//this can actually be done in splineProcedure and probably should
	currentMatrixRow=0;
	vector<double> theChisq; vector<int> theDOF;
	cout << "Checking separate chi_n^2/n in spline fit" << endl;
	cout << "level" << '\t' << "n" << '\t' << "chi_n^2/n" << '\t' << "sqrt(2/n)" << endl;
	for(unsigned int j=0;j < currentAnalysisBins.size();j++)
		{
		int dof = currentAnalysisBins[j].size();
		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, currentMatrixRow, dof);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		theChisq.push_back(chisq); theDOF.push_back(dof);
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

		spline * solution = new spline(intervalBounds[n], splineOrder); solution -> setSpline(theCoefficients,theErrorCoefficients);
		splines.push_back(solution);
		}

	splineArray result(splines);
	result.updateProperties(totalChiSquared,totalDegreesOfFreedom);
	result.updateLevelProperties(theChisq,theDOF);
	return result;
	
}






