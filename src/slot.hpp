/*** LICENCE: ***
Bin histogram method for restoration of smooth functions from noisy integrals. Copyright (C) 2017 Olga Goulko

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.

*** END OF LICENCE ***/

/*** 

This is an edited version of original BHM code so that it may function with 2-dimensional data. 
Edits done in the Fall 2022 by Juan Pablo Varela.

***/

#ifndef SLOT_HPP
#define SLOT_HPP

#include "basic.hpp"
#include <vector>

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
	
	double lowerLeft; /// Lower-left corner of bin
	double upperLeft; /// Upper-left corner of bin
	double lowerRight;  /// Lower-right edge of bin
	double upperRight; ///  Upper-right corner of bin 
	///bool noUpperBound; /// not needed for 2D, or 1D really - more for basis projection 
	
public:
	
	slotBounds(double theLowerLeft = 0.);
	slotBounds(double theLowerLeft, double theUpperLeft, double theLowerRight, double theUpperRight);
	
	bool operator==(const slotBounds& other) const { return ((lowerLeft==other.lowerLeft) && 
		(upperLeft==other.upperLeft) && (lowerRight == other.lowerRight) && (upperRight == other.upperRight)
		); }
	
	double getLowerLeft() const {return lowerLeft;}
	double getLowerRight() const {return lowerLeft;}
	double getUpperLeft() const {return UpperLeft;}
	double getUpperRight() const {return UpperRight;}

	//bool getIsInfinite() const {return noUpperBound;}
	
	bool checkIfInBounds(double valueToCheck) const;
	bool overlapping(const slotBounds& compareBounds) const;
	double slotWidth() const;
	
        std::ostream& printBoundsInfo(std::ostream& strm, verbosity_level_type vlevel) const;
	
};

/// For basis projection method - don't need to modify for thesis
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

        /// Initialize empty slot
	basisSlot(slotBounds theBounds, int theTotalNumOfBasisFn = 0);
        /// Initialize and pre-fill a slot without basis functions
    basisSlot(slotBounds theBounds, long nhits, double theIntegral, double theVariance);

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
	
	void updateEnoughSampled(unsigned int minNumberTimesSampled = defaultMinNumberTimesSampled);
	bool enoughSampled() const {return enoughData;}
	
	double sampledIntegral() const;
	double getVariance() const {return variance;}
	double sampledIntegralVariance() const;
	double sampledIntegralError() const;
	double sampledFunctionValue(double variable) const;
	double sampledFunctionVariance(double variable) const;
	double sampledFunctionError(double variable) const;
	std::vector<double> bareBasisSampledCoeffs() const;
	
	std::ostream& printSlotInfo(std::ostream&) const;
        void printGramSchmidtCoeffs(std::ostream& =std::cout) const;
    void printSampledCoeffs(std::ostream& =std::cout) const;
	
}; 


/// Both taylorSlot and sqrtSlot are child slots of basisSlot related to
/// basis projection - don't need to care about
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
