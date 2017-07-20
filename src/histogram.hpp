#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "basic.hpp" 
#include "slot.hpp" 
#include "spline.hpp"

std::vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps = 1, int totalNumOfBasisFn = 0);

class histogramBasis {
	
private:
	
	excessBin* valuesOutsideBounds;
	double lowerBound;
	double upperBound;
	bool noUpperBound;
	std::vector<basisSlot*> basisSlots;
	
public:
	
	histogramBasis(std::vector<basisSlot*> theBasisSlots);
	~histogramBasis();
	histogramBasis& operator= (const histogramBasis& toBeAssigned);
	histogramBasis (const histogramBasis& toBeCopied);
	
	void appendSlot(basisSlot* theSlot);
	basisSlot* combinedSlot(unsigned int startPoint, unsigned int endPoint) const;
	
	basisSlot* getSlot(unsigned int whichSlot) const;
	unsigned int getSize() const {return basisSlots.size();}
	long getExcessCounter() const {return valuesOutsideBounds->getExcessCounter();}
	double getExcessValues(double norm) const {return valuesOutsideBounds->getExcessValues(norm);}
	
	histogramBasis coarseGrainedHistogram(int minNumberTimesSampled = defaultMinNumberTimesSampled);
	
	void sample(double variable, double valueToSample);
	void sampleUniform(double variable, double valueToSample);
	std::pair<double,double> sampledFunctionValueAverage(double variable) const;
	std::pair<double,double> sampledFunctionValueWeightedAverage(double variable) const;
	
	histogramBasis scaledHistogram(long norm);
	histogramBasis normalizedHistogram(double norm);
	bool addAnotherHistogram(histogramBasis anotherHistogram);
	
	std::vector< std::vector< basisSlot* > > binHierarchy(long int norm);
	splineArray BHMfit(unsigned int splineOrder, unsigned int minLevel, long norm, double fitAcceptanceThreshold, double jumpSuppression);
	
};

splineArray matchedSplineFit(std::vector< std::vector< basisSlot* > > currentAnalysisBins, std::vector< slotBounds > intervalBounds, unsigned int splineOrder, double jumpSuppression, std::vector< double > aMaxVector, std::vector< double > chisqArray);

#endif