#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "basic.hpp" 
#include "slot.hpp" 

std::vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps = 10, int totalNumOfBasisFn = 5);

class histogramBasis {
	
private:
	
	excessBin* valuesOutsideBounds;
	std::vector<basisSlot*> basisSlots;
	
public:
	
	histogramBasis(std::vector<basisSlot*> theBasisSlots);
	
	void appendSlot(basisSlot* theSlot);
	
	void sample(double variable, double valueToSample);
	std::pair<double,double> sampledFunctionValueAverage(double variable);
	std::pair<double,double> sampledFunctionValueWeightedAverage(double variable);
	
	//void addAnotherHistogram(histogramBasis anotherHistogram);
	//void subtractAnotherHistogram(histogramBasis anotherHistogram);
	void scale(double norm);
};


#endif