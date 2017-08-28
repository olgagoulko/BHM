#ifndef HISTOGRAM_HPP
#define HISTOGRAM_HPP

#include "basic.hpp" 
#include "slot.hpp" 
#include "spline.hpp"

#include <iosfwd>

struct fitAcceptanceThreshold {
	fitAcceptanceThreshold() : max(0), steps(0) {}
	double min;
	double max;
	int steps;
} ;

/// Read basis slots from a formatted stream (a file)
std::vector<basisSlot*> readBasisSlots(std::istream&);

std::vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps = 1, int totalNumOfBasisFn = 0);

class histogramBasis {
	
private:
	
	excessBin* valuesOutsideBounds;
	double lowerBound;
	double upperBound;
	bool noUpperBound;
	std::vector<basisSlot*> basisSlots;
        unsigned long numberOfSamples;
	
public:

        class ConsistentWithZero_Error : public std::runtime_error {
          public:
            ConsistentWithZero_Error() : std::runtime_error("Data consistent with zero")
            {}
        };

        class NotEnoughData_Error : public std::runtime_error {
          public:
            NotEnoughData_Error() : std::runtime_error("Not enough data for analysis")
            {}
        };

        class InvalidFileFormatError : public std::runtime_error {
          public:
            explicit InvalidFileFormatError(const std::string& reason) :
                std::runtime_error("Invalid format of the input file: "+reason)
            {}
        };

        class OverlappingSlot : public std::runtime_error {
          public:
            explicit OverlappingSlot(const std::string& reason) :
                std::runtime_error("Overlapping slot when reading from file: "+reason)
            {}
        };

    
        /// Create the histogram from the vector of slots (e.g., obtained via `generateBasisSlots()`
        /** @warning
            The histogram object takes the ownership of the pointers to slots. Caller must not deallocate the pointers.
            @bug
            The pointers must be heap-allocated via `new`, otherwise the behavior is undefined.
        */
	histogramBasis(std::vector<basisSlot*> theBasisSlots);

        /// Construct from a formatted stream (e.g., a file)
        histogramBasis(std::istream&);

        ~histogramBasis();

        
    
	histogramBasis& operator= (const histogramBasis& toBeAssigned);
	histogramBasis (const histogramBasis& toBeCopied);

        /// Add the slot to the histogram.
        /** @warning
            The histogram object takes the ownership of the pointer. Caller must not deallocate the pointer.
            @bug
            The pointer must be heap-allocated via `new`, otherwise the behavior is undefined.
        */
	void appendSlot(basisSlot* theSlot);
	basisSlot* combinedSlot(unsigned int startPoint, unsigned int endPoint) const;
	
	basisSlot* getSlot(unsigned int whichSlot) const;
	unsigned int getSize() const {return basisSlots.size();}
        double getLowerBound() const { return lowerBound; }
        double getUpperBound() const { return upperBound; } // FIXME: what should it return for "no upper bound"? +INF?
        bool hasUpperBound() const { return !noUpperBound; }
        unsigned long getNumberOfSamples() const { return numberOfSamples; } // FIXME: does NOT count out-of-bounds samples 
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
        splineArray BHMfit(unsigned int splineOrder, unsigned int minLevel, long norm, fitAcceptanceThreshold theThreshold, double jumpSuppression, bool verbose, bool fail_if_zero);

    // DEBUG!
    friend std::ostream& operator<<(std::ostream& , const histogramBasis&);
	
};
std::ostream& operator<<(std::ostream& , const histogramBasis&);


splineArray matchedSplineFit(std::vector< std::vector< basisSlot* > > currentAnalysisBins, std::vector< slotBounds > intervalBounds, unsigned int splineOrder, double jumpSuppression, std::vector< double > aMaxVector, std::vector< double > chisqArray);

#endif
