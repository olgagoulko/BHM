#include "basic.hpp"
#include "slot.hpp"
#include "histogram.hpp"

using namespace std; 

vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps, int totalNumOfBasisFn)
{
	if(numberOverlaps<1) numberOverlaps = 10;
	
	if(totalNumOfBasisFn<1) totalNumOfBasisFn=1;
	
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
	for(unsigned int i=0;i<theBasisSlots.size();i++) basisSlots.push_back(theBasisSlots[i]);
	
}

void histogramBasis::appendSlot(basisSlot* theSlot)
{
	basisSlots.push_back(theSlot);
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

pair<double,double> histogramBasis::sampledFunctionValueAverage(double variable)
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

pair<double,double> histogramBasis::sampledFunctionValueWeightedAverage(double variable)
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









