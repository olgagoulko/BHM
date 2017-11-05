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
#ifndef SPLINE_HPP
#define SPLINE_HPP

#include "basic.hpp"  
#include "slot.hpp"

class splinePiece {
	
private:
	
	slotBounds bounds;
	unsigned int splineOrder;			//this is total number of spline parameters (m+1 in the paper), e.g. splineOrder=4 for cubic polynomial
	std::vector<double>splineCoefficients;
	std::vector<double>splineErrorCoefficients;
	
public:
	
	splinePiece(slotBounds theBounds, unsigned int theSplineOrder = 4);
	splinePiece* Clone() { return new splinePiece(*this);}
	
	void setSplinePiece(std::vector< double > theCoefficients, std::vector< double > theCovariance);
	
	std::vector<double> getCoefficients() const {return splineCoefficients;}
	std::vector<double> getErrorCoefficients() const {return splineErrorCoefficients;}
	unsigned int getSplineOrder() const {return splineOrder;}
	slotBounds getBounds() const {return bounds;}
	
	bool checkIntervalAcceptance(std::vector< std::vector< basisSlot* > > currentAnalysisBins, double currentFitAcceptanceThreshold, double chisqArrayElement, unsigned int intervalOrder, bool checkIntervals, double usableBinFraction) const;
	
	double splineValue(double variable) const;
	double splineError(double variable) const;
	double splineIntegral(slotBounds theBounds) const;
	double splineDerivative(double variable, unsigned int derivativeOrder) const;
	
	std::ostream& printSplinePiece(std::ostream&, verbosity_level_type verb=VERBOSE) const;
	
};


class splineArray {
	
private:
	
	double lowerBound;
	double upperBound;
	unsigned int splineOrder;
	std::vector<splinePiece*> splines;
	std::vector<double> intervalBoundaries;
	std::vector<double> levelsChiSquared;
	std::vector<int> levelsDegreesOfFreedom;
	double currentThreshold;
	bool acceptableSpline;
	
public:
	
	splineArray();
	splineArray(std::vector< splinePiece* > theSplines);
	~splineArray();
	splineArray& operator= (const splineArray& toBeAssigned);
	splineArray (const splineArray& toBeCopied);
	
	void updateLevelProperties(std::vector<double> theChisq, std::vector<int> theDOF);
	void updateGoodness(bool acceptable, double threshold);
	
	splinePiece * getSplinePiece(unsigned int whichPiece) const;
	bool getAcceptance() const {return acceptableSpline;}
	double getThreshold() const {return currentThreshold;}
	bool checkOverallAcceptance(double currentFitAcceptanceThreshold) const;
	std::vector<slotBounds> getBounds() const;
	int numberKnots() const {return intervalBoundaries.size();}
	
	double splineValue(double variable) const;
	double splineError(double variable) const;
	double splineDerivative(double variable, unsigned int derivativeOrder) const;

        double getLowerBound() const { return lowerBound; }
        double getUpperBound() const { return upperBound; }
    
        std::ostream& printSplineArrayInfo(std::ostream&) const;
        std::ostream& printSplines(std::ostream&) const;
	
};


double splineBasisFunction(slotBounds theBounds, unsigned int functionOrder);

double* binomialMatrix(int vectorSize);
double* binomialVector(int vectorSize);

void solveLinearEquation(double* matrix, double* vec, unsigned int vectorSize);
void solveSVD(gsl_matrix* A, gsl_matrix* V, gsl_vector* b, gsl_vector* diagonal, gsl_vector* x);

double* designMatrix(std::vector< basisSlot* > slotArray, unsigned int splineOrder);
double* designMatrix(std::vector< std::vector< basisSlot*> > slotArray, unsigned int splineOrder);
double* designMatrixProduct(double* matrix, unsigned int numberOfSlots, unsigned int splineOrder);
double* integralVector(std::vector< basisSlot* > slotArray);
double* integralVector(std::vector< std::vector< basisSlot*> > slotArray);
double* integralVectorProduct(double* matrix, double* vec, unsigned int numberOfSlots, unsigned int splineOrder);

bool isDataConsistentWithZero(std::vector< std::vector<basisSlot*> > slotArray);


#endif
