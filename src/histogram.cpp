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
	while(minVar+(counter*offset)+slotWidth<maxVar)
		{
		slotBounds currentBounds(minVar+(counter*offset), minVar+(counter*offset)+slotWidth);
		basisSlotVector.push_back(new taylorSlot(currentBounds, totalNumOfBasisFn));
		counter++;
		if((minVar+(counter*offset)+slotWidth>=maxVar)&&(minVar+(counter*offset)<maxVar)) 
			{
			slotBounds finalBounds(minVar+(counter*offset), maxVar);
			basisSlotVector.push_back(new taylorSlot(finalBounds, totalNumOfBasisFn));
			}
		}
	
	return basisSlotVector;
}


histogramBasis::histogramBasis(vector<basisSlot*> theBasisSlots)
{
	
	valuesOutsideBounds = new excessBin();
	for(int i=0;i<theBasisSlots.size();i++) basisSlots.push_back(theBasisSlots[i]);
	
}

void histogramBasis::sample(double variable, double valueToSample)
{
	bool isOutsideBounds=true;
	for(int i=0;i<basisSlots.size();i++)
		{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) ) {basisSlots[i] -> sample(variable,valueToSample); isOutsideBounds=false;}
		}
	if(isOutsideBounds) valuesOutsideBounds -> sample(variable, valueToSample);
}

double histogramBasis::sampledFunctionValueAverage(double variable)
{
	double sampledValue=0;
	int numberSlots=0;
	for(int i=0;i<basisSlots.size();i++)
		{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) ) {sampledValue += basisSlots[i] -> sampledFunctionValue(variable); numberSlots++;}
		}
	
	if(numberSlots==0) sampledValue=valuesOutsideBounds -> getExcessValues();
	else sampledValue=sampledValue/double(numberSlots);
	
	return sampledValue;
}










