#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "basic.hpp" 
#include "slot.hpp" 
#include "matrix.hpp"

std::vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps = 10, int totalNumOfBasisFn = 5);

class histogramBasis {
	
private:
	
	excessBin* valuesOutsideBounds;
	double lowerBound;
	double upperBound;
	bool noUpperBound;
	std::vector<basisSlot*> basisSlots;
	
public:
	
	histogramBasis(std::vector<basisSlot*> theBasisSlots);
	
	void appendSlot(basisSlot* theSlot);
	basisSlot* combinedSlot(unsigned int startPoint, unsigned int endPoint) const;
	basisSlot* getSlot(unsigned int whichSlot) const;
	unsigned int getSize() const {return basisSlots.size();}
	histogramBasis coarseGrainedHistogram(int minNumberTimesSampled = 10);
	
	void sample(double variable, double valueToSample);
	std::pair<double,double> sampledFunctionValueAverage(double variable) const;
	std::pair<double,double> sampledFunctionValueWeightedAverage(double variable) const;
	
	histogramBasis scaledHistogram(long norm);
	bool addAnotherHistogram(histogramBasis anotherHistogram);
	
	std::vector< std::vector< basisSlot* > > binHierarchy(long int norm);
	splineArray splineProcedure(unsigned int splineOrder, unsigned int minLevel, long norm, double fitAcceptanceThreshold, double gluingFactor);
};

splineArray matchedSplineFit(std::vector< std::vector< basisSlot* > > currentAnalysisBins, std::vector< slotBounds > intervalBounds, unsigned int splineOrder, double gluingFactor, std::vector< double > a3array, std::vector< double > chisqArray);

#endif