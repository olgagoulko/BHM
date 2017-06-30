#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "basic.hpp"  
#include "slot.hpp"

class spline {
	
private:
	
	slotBounds bounds;
	unsigned int splineOrder;
	std::vector<double>splineCoefficients;
	std::vector<double>splineErrorCoefficients;
	double chiSquared;
	int degreesOfFreedom;
	
public:
	
	spline(slotBounds theBounds, unsigned int theSplineOrder = 4);
	
	void setSpline(std::vector< double > theCoefficients, std::vector< double > theCovariance, double theChiSquared = 0, int theDOF = 1);
	
	std::vector<double> getCoefficients() const {return splineCoefficients;}
	std::vector<double> getErrorCoefficients() const {return splineErrorCoefficients;}
	unsigned int getSplineOrder() const {return splineOrder;}
	slotBounds getBounds() const {return bounds;}
	double getChisquared() const {return chiSquared/double(degreesOfFreedom);}
	int getDOF() const {return degreesOfFreedom;}
	double goodnessOfFitQ() const;
	
	bool isSplineGood(double fitAcceptanceThreshold) const {return goodnessOfFitQ()>fitAcceptanceThreshold;}
	
	double splineValue(double variable) const;
	double splineError(double variable) const;
	double splineIntegral(slotBounds theBounds) const;
	double splineDerivative(double variable, unsigned int derivativeOrder) const;
	
	void printSpline() const;
	
};


class splineArray {
	
private:
	
	double lowerBound;
	double upperBound;
	unsigned int splineOrder;
	std::vector<spline*> splines;
	std::vector<double> intervalBoundaries;
	double totalChiSquared;
	std::vector<double> levelsChiSquared;
	int totalDegreesOfFreedom;
	std::vector<int> levelsDegreesOfFreedom;
	bool accetableSpline;
	
public:
	
	splineArray();
	splineArray(std::vector< spline* > theSplines);
	void updateProperties(double theTotalChisq, int theTotalDOF);
	void updateLevelProperties(std::vector<double> theChisq, std::vector<int> theDOF);
	void updateGoodness(bool acceptable);
	
	spline * getSpline(unsigned int whichSpline) const;
	double getChisquared() const {return totalChiSquared/double(totalDegreesOfFreedom);}
	double goodnessOfFitQ() const;
	bool getAcceptance() const {return accetableSpline;}
	bool checkOverallAcceptance(double fitAcceptanceThreshold) const;
	std::vector<slotBounds> getBounds() const;
	int numberKnots() const {return intervalBoundaries.size();}
	
	double splineValue(double variable) const;
	double splineError(double variable) const;
	double splineDerivative(double variable, unsigned int derivativeOrder) const;
	
	void printSplineArrayInfo() const;
	void printSplines() const;
	
};


double splineBasisFunction(slotBounds theBounds, unsigned int functionOrder);

double* binomialMatrix(int vectorSize);
double* binomialVector(int vectorSize);

void solveLinearEquation(double* matrix, double* vec, unsigned int vectorSize);
double solveSVD(gsl_matrix* A, gsl_matrix* V, gsl_vector* b, gsl_vector* diagonal, gsl_vector* x);

double* secondDerivativeMatrix(unsigned int vectorSize);

double* designMatrix(std::vector< basisSlot* > slotArray, unsigned int splineOrder);
double* designMatrix(std::vector< std::vector< basisSlot*> > slotArray, unsigned int splineOrder);
double* designMatrixProduct(double* matrix, unsigned int numberOfSlots, unsigned int splineOrder);
double* integralVector(std::vector< basisSlot* > slotArray);
double* integralVector(std::vector< std::vector< basisSlot*> > slotArray);
double* integralVectorProduct(double* matrix, double* vec, unsigned int numberOfSlots, unsigned int splineOrder);

bool isDataConsistentWithZero(std::vector< std::vector<basisSlot*> > slotArray);


#endif