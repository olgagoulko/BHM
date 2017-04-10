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
	
	spline oneSplineFit(unsigned int splineOrder = 4);
	splineArray splineFit(std::vector<unsigned int> intervalBoundaries, double gluingFactor, unsigned int splineOrder = 4);
	splineArray matchedSplineFitSmallBins(std::vector<unsigned int> intervalBoundaries, unsigned int splineOrder = 4);
	splineArray splineProcedureSmallBins(unsigned int splineOrder, unsigned int minNumberBins, double fitAcceptanceThreshold);
	splineArray splineProcedure(unsigned int splineOrder, unsigned int minLevel, long norm, double fitAcceptanceThreshold);
	splineArray splineProcedureOverlap(unsigned int splineOrder, long int norm, double fitAcceptanceThreshold);
};

splineArray matchedSplineFit(std::vector< std::vector<basisSlot*> > currentAnalysisBins, std::vector<slotBounds> intervalBounds, unsigned int splineOrder = 4);
splineArray matchedSplineFitOverlap(std::vector< std::vector<basisSlot*> > currentAnalysisBins, std::vector<slotBounds> intervalBounds, unsigned int splineOrder = 4);

#endif