#ifndef SLOT_HPP
#define SLOT_HPP

#include "basic.hpp"

class excessBin {

private:

	long excessCounter;
	double sumExcessValues;

public:
	
	excessBin();
	excessBin(long theExcessCounter, double theSumExcessValues);
	
	int getExcessCounter() const {return excessCounter;}
	double getExcessValues(double norm = 1.) const {return sumExcessValues/norm;}
	
	void sample(double variable, double valueToSample);
	
};

class slotBounds {
	
protected:
	
	double lowerBound;
	double upperBound;
	bool noUpperBound;
	
public:
	
	slotBounds(double theLowerBound = 0.);
	slotBounds(double theLowerBound, double theUpperBound);
	
	bool operator==(const slotBounds& other) const { return ((lowerBound==other.lowerBound) && (noUpperBound==other.noUpperBound) && ( (upperBound==other.upperBound)||!noUpperBound) ); }
	
	double getLowerBound() const {return lowerBound;}
	double getUpperBound() const {return upperBound;}
	bool getIsInfinite() const {return noUpperBound;}
	
	bool checkIfInBounds(double valueToCheck) const;
	bool overlapping(const slotBounds& compareBounds) const;
	double slotWidth() const;
	
	void printBoundsInfo() const;
	
};

class basisSlot {
	
protected:
	
	slotBounds bounds;
	int totalNumOfBasisFn;
	
	std::vector< std::vector<double> >GramSchmidtCoeffs;
	std::vector<double>sampledCoeffsValues;
	std::vector<double>sampledCoeffsVariance;
	double integral;
	double variance;
	long numberTimesSampled;
	bool enoughData;
	
public:
	
	basisSlot(slotBounds theBounds, int theTotalNumOfBasisFn = 0);
	virtual basisSlot* Clone() { return new basisSlot(*this);}
	virtual ~basisSlot() {}
	
	virtual void initializeGramSchmidt();
	
	virtual double GramSchmidtBasisFn(int numOfBasisFn, double variable) const {return 1.;}
	virtual double bareBasisFn(int numOfBasisFn, double variable) const {return 1.;};
	virtual double pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const {return 0.;}
	virtual double weight(double variable) const {return 1.;}
	virtual bool checkIfInBasisSlot(double variable) const {return bounds.checkIfInBounds(variable);}
	
	void sample(double variable, double valueToSample);
	bool addAnotherSlot(basisSlot* anotherSlot, int sign = 1);
	void combineWithSlot(basisSlot* anotherSlot);
	void scale(long norm);
	void normalize(double norm);
	
	long getNumberTimesSampled() const {return numberTimesSampled;}
	slotBounds getBounds() const {return bounds;}
	
	void updateEnoughSampled(int minNumberTimesSampled = defaultMinNumberTimesSampled);
	bool enoughSampled() const {return enoughData;}
	
	double sampledIntegral() const;
	double sampledIntegralVariance() const;
	double sampledIntegralError() const;
	double sampledFunctionValue(double variable) const;
	double sampledFunctionVariance(double variable) const;
	double sampledFunctionError(double variable) const;
	std::vector<double> bareBasisSampledCoeffs() const;
	
	void printSlotInfo() const;
	void printGramSchmidtCoeffs() const;
	void printSampledCoeffs() const;
	
}; 

class taylorSlot: public basisSlot {
	
private:
	
public:
	
	taylorSlot(slotBounds theBounds, int theTotalNumOfBasisFn = 0);
	taylorSlot* Clone() { return new taylorSlot(*this);}
	~taylorSlot() {}
	
	void initializeGramSchmidt();
	
	double GramSchmidtBasisFn(int numOfBasisFn, double variable) const;
	double bareBasisFn(int numOfBasisFn, double variable) const;
	double pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const;
	
};

class sqrtSlot: public basisSlot {
	
private:
	
public:
	
	sqrtSlot(slotBounds theBounds, int theTotalNumOfBasisFn = 0);
	sqrtSlot* Clone() { return new sqrtSlot(*this);}
	~sqrtSlot() {}
	
	double GramSchmidtBasisFn(int numOfBasisFn, double variable) const;
	double bareBasisFn(int numOfBasisFn, double variable) const;
	double weight(double variable) const;
	double pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const;
	
};


#endif
