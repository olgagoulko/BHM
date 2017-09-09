#include "basic.hpp"
#include "slot.hpp"

#include <stdexcept>

using namespace std;

//----------------------------------------------------------------------------------------------------------------------------------

excessBin::excessBin()
	{
	excessCounter=0;
	sumExcessValues=0;
	}
	
excessBin::excessBin(long theExcessCounter, double theSumExcessValues)
	{
	excessCounter=theExcessCounter;
	sumExcessValues=theSumExcessValues;
	}

void excessBin::sample(double variable, double valueToSample)
	{
	excessCounter++;
	sumExcessValues+=valueToSample;
	}

//----------------------------------------------------------------------------------------------------------------------------------

slotBounds::slotBounds(double theLowerBound)
	{
	lowerBound=theLowerBound;
	upperBound=lowerBound;
	noUpperBound=true;
	}

slotBounds::slotBounds(double theLowerBound, double theUpperBound)
	{
	lowerBound=theLowerBound;
	upperBound=theUpperBound;
	noUpperBound=false;
	}

bool slotBounds::checkIfInBounds(double valueToCheck) const
	{
	bool isInBounds=true;
	if(valueToCheck<lowerBound-VERY_SMALL_NUMBER) isInBounds=false;
	else if(noUpperBound==false) {
		if(valueToCheck>=upperBound) isInBounds=false;
		}
	return isInBounds;
	}

bool slotBounds::overlapping(const slotBounds& compareBounds) const
	{
	bool overlapping = ( checkIfInBounds(compareBounds.getLowerBound()) && (noUpperBound || !isAround(upperBound,compareBounds.getLowerBound())) );
	if(overlapping == false) overlapping = ( compareBounds.checkIfInBounds(lowerBound) && (compareBounds.getIsInfinite() || !isAround(lowerBound,compareBounds.getUpperBound())) );
	return overlapping;
	}

double slotBounds::slotWidth() const
	{
	double width=0;
	if(noUpperBound==false) {width= upperBound-lowerBound;}
	else {cerr << "ERROR: infinite slot width" << endl; throw std::runtime_error("ERROR: infinite slot width");}
	return width;
	}

std::ostream& slotBounds::printBoundsInfo(std::ostream& strm, verbosity_level_type vlevel) const
	{
        if (vlevel==VERBOSE) {
            strm << "lowerBound = " << fixed << setprecision(16) << lowerBound << endl;
            if (noUpperBound) {
                strm << "no upper bound" << endl;
            } else {
                strm << "upperBound = " << upperBound << endl;
            }
        } else /* not VERBOSE */ {
            strm << scientific << setprecision(6)
                 << lowerBound
                 << " ";
            if (noUpperBound) {
                strm << "INF" << endl;
            } else {
                strm << upperBound << endl;
            }
        }
        return strm;
        }


//----------------------------------------------------------------------------------------------------------------------------------

basisSlot::basisSlot(slotBounds theBounds, long nhits, double theIntegral, double theVariance)
    : bounds(theBounds), totalNumOfBasisFn(0),
      integral(theIntegral), variance(theVariance), numberTimesSampled(nhits),
      enoughData(false) // FIXME: is this correct?
	{	}

basisSlot::basisSlot(slotBounds theBounds, int theTotalNumOfBasisFn)
    : bounds(theBounds), totalNumOfBasisFn(theTotalNumOfBasisFn>0? theTotalNumOfBasisFn : 0),
      integral(0), variance(0), numberTimesSampled(0), enoughData(false)
	{
	GramSchmidtCoeffs.resize(totalNumOfBasisFn);
	for(int i=0; i<totalNumOfBasisFn; i++)
		{
		sampledCoeffsValues.push_back(0);
		sampledCoeffsVariance.push_back(0);
		for(int j=0; j<totalNumOfBasisFn; j++)
			{
			if (i == j) GramSchmidtCoeffs[i].push_back(1.);
			else GramSchmidtCoeffs[i].push_back(0);
			}
		}
	}

void basisSlot::initializeGramSchmidt()
	{
	//modified Gram-Schmidt (numerically more stable, but not as good as Householder method)
	for(int i=0; i<totalNumOfBasisFn; i++)
		{
		double newnorm=0;
		for(int j1=0; j1<=i; j1++)
			for(int j2=0; j2<=i; j2++)
				{
				newnorm+=GramSchmidtCoeffs[i][j1]*GramSchmidtCoeffs[i][j2]*pairwiseIntegral(j1,j2);
				}
		if(newnorm<=0) {
                    cerr << "ERROR: negative norm in Gram-Schmidt initialization in function number " << i << endl;
                    printSlotInfo(std::cerr);
                    throw std::runtime_error("ERROR: negative norm in Gram-Schmidt initialization in function number");
                }
		newnorm=sqrt(newnorm);
		for(int j1=0; j1<=i; j1++) GramSchmidtCoeffs[i][j1]=GramSchmidtCoeffs[i][j1]/newnorm;
		for(int j1=i+1; j1<totalNumOfBasisFn; j1++)
			{
			newnorm=0;
			for(int k=0; k<=i; k++)
				for(int m=0; m<=j1; m++)
					{
					newnorm+=GramSchmidtCoeffs[i][k]*GramSchmidtCoeffs[j1][m]*pairwiseIntegral(k,m);
					}
			for(int j2=0; j2<=i; j2++)
				{
				GramSchmidtCoeffs[j1][j2]-=GramSchmidtCoeffs[i][j2]*newnorm;
				}
			}
		}
	}


void basisSlot::sample(double variable, double valueToSample)
	{
	if(this -> checkIfInBasisSlot(variable))
		{
		numberTimesSampled++;
		double samplingTotal; double delta;    
		delta=valueToSample-integral;
		integral+=delta/double(numberTimesSampled);
		variance+=delta*(valueToSample-integral);
		for(int j=0; j<totalNumOfBasisFn; j++)
			{
			samplingTotal=weight(variable)*valueToSample*GramSchmidtBasisFn(j, variable);
			delta=samplingTotal-sampledCoeffsValues[j];
			sampledCoeffsValues[j]+=delta/double(numberTimesSampled);
			sampledCoeffsVariance[j]+=delta*(samplingTotal-sampledCoeffsValues[j]);
			}
		}
	else {cerr << "ERROR: trying to sample in the wrong slot" << endl; throw std::runtime_error("ERROR: trying to sample in the wrong slot"); }
	}

	
bool basisSlot::addAnotherSlot(basisSlot* anotherSlot, int sign)
	{
	bool addable=false;
	if( (totalNumOfBasisFn==(anotherSlot -> totalNumOfBasisFn)) && (bounds==(anotherSlot -> bounds)) && (numberTimesSampled+sign*(anotherSlot -> numberTimesSampled)>=0))
		{
		addable=true;
		long newNumberTimesSampled=numberTimesSampled+sign*(anotherSlot -> numberTimesSampled);
		if(newNumberTimesSampled>0)
			{
			double auxiliary=integral-(anotherSlot -> integral); auxiliary=auxiliary*auxiliary;
			auxiliary*=numberTimesSampled*(anotherSlot -> numberTimesSampled);
			integral=(integral*numberTimesSampled+sign*(anotherSlot -> integral)*(anotherSlot -> numberTimesSampled))/double(newNumberTimesSampled);
			variance+=sign*(anotherSlot -> variance)+sign*auxiliary/double(newNumberTimesSampled);
			for(int j=0; j<totalNumOfBasisFn; j++)
				{
				sampledCoeffsValues[j]*=numberTimesSampled;
				sampledCoeffsValues[j]+=sign*(anotherSlot -> sampledCoeffsValues)[j]*(anotherSlot -> numberTimesSampled);
				sampledCoeffsValues[j]*=1./double(newNumberTimesSampled);
			
				auxiliary=sampledCoeffsValues[j]-(anotherSlot -> sampledCoeffsValues)[j]; auxiliary=auxiliary*auxiliary;
				auxiliary*=numberTimesSampled*(anotherSlot -> numberTimesSampled);
				sampledCoeffsVariance[j]+=sign*(anotherSlot -> sampledCoeffsVariance)[j]+sign*auxiliary/double(newNumberTimesSampled);
				}
			}
		numberTimesSampled=newNumberTimesSampled;
		}
	return addable;
	}

void basisSlot::combineWithSlot(basisSlot* anotherSlot)
	{
	//only combines integrals, not basis functions for any type of slot
	double auxiliary=integral-(anotherSlot -> integral); auxiliary=auxiliary*auxiliary;
	auxiliary*=numberTimesSampled*(anotherSlot -> numberTimesSampled);
	integral=integral*numberTimesSampled+(anotherSlot -> integral)*(anotherSlot -> numberTimesSampled);
	variance+=(anotherSlot -> variance);
	numberTimesSampled+=(anotherSlot -> numberTimesSampled);
	if(numberTimesSampled>0) integral=integral/double(numberTimesSampled);
	if(numberTimesSampled>0) variance+=auxiliary/double(numberTimesSampled);
	}


void basisSlot::scale(long norm)
	{
	//adding norm-numberTimesSampled zeros to the sample retrospectively
	//this way we don't have to sample 0 in every slot at every step
	if(norm>numberTimesSampled)
		{
		double normScaling=numberTimesSampled/double(norm);
		double varScaling=normScaling*(norm-numberTimesSampled);
		
		variance+=integral*integral*varScaling;
		integral*=normScaling;
		
		for(int j=0; j<totalNumOfBasisFn; j++)
			{
			sampledCoeffsVariance[j]+=sampledCoeffsValues[j]*sampledCoeffsValues[j]*varScaling;
			sampledCoeffsValues[j]*=normScaling;
			}
		numberTimesSampled=norm;
		}
	}

void basisSlot::normalize(double norm)
	{
	if(norm>0)
		{			
		variance*=1./norm/norm;
		integral*=1./norm;
			
		for(int j=0; j<totalNumOfBasisFn; j++)
			{
			sampledCoeffsVariance[j]*=1./norm/norm;
			sampledCoeffsValues[j]*=1./norm;
			}
		}
	else {cerr << "ERROR: norm not > 0" << endl; throw std::runtime_error("ERROR: norm not > 0");}
	}


void basisSlot::updateEnoughSampled(int minNumberTimesSampled)
	{
	if(numberTimesSampled>=minNumberTimesSampled) enoughData=true;
	else enoughData=false;
	}

double basisSlot::sampledIntegral() const
	{
	return integral;
	}


double basisSlot::sampledIntegralVariance() const
	{
	double result=0;
	if(numberTimesSampled>1) result=variance/double(numberTimesSampled-1);
	return result;
	}


double basisSlot::sampledIntegralError() const
	{
	double result=0;
	if(numberTimesSampled>1) result= sqrt(sampledIntegralVariance()/double(numberTimesSampled));
	return result;
	}


double basisSlot::sampledFunctionValue(double variable) const
	{
	double result=0;
	for(int i=0; i<totalNumOfBasisFn; i++) result+=sampledCoeffsValues[i]*GramSchmidtBasisFn(i,variable);
    
	if(totalNumOfBasisFn==0) result=integral/bounds.slotWidth();
    
	return result;
	}

double basisSlot::sampledFunctionVariance(double variable) const
	{
	double result=0;
	if(numberTimesSampled>1)
		{
		for(int i=0; i<totalNumOfBasisFn; i++) {result+=sampledCoeffsVariance[i]*GramSchmidtBasisFn(i,variable)*GramSchmidtBasisFn(i,variable)/double(numberTimesSampled-1);}
		if(totalNumOfBasisFn==0) result=variance/double(numberTimesSampled-1)/bounds.slotWidth()/bounds.slotWidth();
		}
	return result;
	}

double basisSlot::sampledFunctionError(double variable) const
	{
	double result=0;
	if(numberTimesSampled>1) result = sqrt(sampledFunctionVariance(variable)/double(numberTimesSampled));
	return result;
	}

vector<double> basisSlot::bareBasisSampledCoeffs() const
	{
	vector<double> result(totalNumOfBasisFn);

	for(int i=0; i<totalNumOfBasisFn; i++)
		{
		result[i]=0;
		for(int j=i; j<totalNumOfBasisFn; j++)
			result[i]+=sampledCoeffsValues[j]*GramSchmidtCoeffs[j][i];
		}
	return result;
	}

std::ostream& basisSlot::printSlotInfo(std::ostream& strm) const
	{
        bounds.printBoundsInfo(strm, VERBOSE);
	strm << "totalNumOfBasisFn = " << totalNumOfBasisFn << endl;
        return strm;
	}


void basisSlot::printGramSchmidtCoeffs(std::ostream& strm) const
	{
	for(int i=0; i<totalNumOfBasisFn; i++)
		{
		for(int j=0; j<totalNumOfBasisFn; j++)
			{
			strm << fixed << setprecision(12) << GramSchmidtCoeffs[i][j] << '\t';
			}
		strm << endl;
		}
	}

void basisSlot::printSampledCoeffs(std::ostream& strm) const
	{
	strm << integral << '\t' << variance << '\t' << numberTimesSampled << endl;
	for(int i=0; i<totalNumOfBasisFn; i++) strm << sampledCoeffsValues[i] << '\t' << sampledCoeffsVariance[i] << endl;
	}

//----------------------------------------------------------------------------------------------------------------------------------

taylorSlot::taylorSlot(slotBounds theBounds, int theTotalNumOfBasisFn) : basisSlot(theBounds, theTotalNumOfBasisFn)
	{
	if(theBounds.getIsInfinite()==true) {cerr << "ERROR: trying to initialize Taylor expansion slot with open bounds" << endl; throw std::runtime_error("ERROR: trying to initialize Taylor expansion slot with open bounds");}
	initializeGramSchmidt();
	}

void taylorSlot::initializeGramSchmidt()
	{
	//we don't need initializeGramSchmidt() here, because we are calculating the Legendre polynomials explicitly
	}


double taylorSlot::bareBasisFn(int numOfBasisFn, double variable) const
	{      
	/* 1 x x^2 ... x^n*/
	return pow(variable,numOfBasisFn);
	}


double taylorSlot::GramSchmidtBasisFn(int numOfBasisFn, double variable) const
	{
	double shiftedvariable=2*variable-(bounds.getUpperBound()+bounds.getLowerBound());
	shiftedvariable=shiftedvariable/bounds.slotWidth();
	
	if(abs(shiftedvariable)>1)
		{
		if(isAround(abs(shiftedvariable),1)) {if(shiftedvariable>1) shiftedvariable=1; else shiftedvariable=-1;}
		else {
                    std::cerr << "ERROR: wrong variable in taylorSlot; variable = " << variable
                              << ", number basis function = " << numOfBasisFn
                              << endl;
                    printSlotInfo(std::cerr);
                    throw std::runtime_error("ERROR: wrong variable in taylorSlot");
                }
		}
	
	return sqrt(2*numOfBasisFn+1)*gsl_sf_legendre_Pl (numOfBasisFn, shiftedvariable)/sqrt(bounds.slotWidth());
	}


double taylorSlot::pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const
	{
	int sumExp=numOfBasisFn1+numOfBasisFn2+1;
	return (pow(bounds.getUpperBound(),sumExp)-pow(bounds.getLowerBound(),sumExp))/double(sumExp);
	}


//----------------------------------------------------------------------------------------------------------------------------------

sqrtSlot::sqrtSlot(slotBounds theBounds, int theTotalNumOfBasisFn) : basisSlot(theBounds, theTotalNumOfBasisFn)
	{
	if(theBounds.getIsInfinite()==true) {cerr << "ERROR: trying to initialize square root expansion slot with open bounds" << endl; throw std::runtime_error("ERROR: trying to initialize square root expansion slot with open bounds");}
	initializeGramSchmidt();
	}


double sqrtSlot::bareBasisFn(int numOfBasisFn, double variable) const
	{       
	/* x^-1/2
	 1                                   *
	 x^1/2
	 x
	 x^3/2... */
	//double result=1./sqrt(variable);
	//result*=pow(sqrt(variable),numOfBasisFn);
	
	/* x^-1/2
	 x^1/2
	 x^3/2... */
	double result=1./sqrt(variable);
	result*=pow(variable,numOfBasisFn);
	
	return result;
	}


double sqrtSlot::GramSchmidtBasisFn(int numOfBasisFn, double variable) const
	{
	double result=0;
	for(int j=0; j<=numOfBasisFn; j++)
		{
		result+=GramSchmidtCoeffs[numOfBasisFn][j]*bareBasisFn(j,variable);
		}
	return result;
	}

double sqrtSlot::weight(double variable) const
	{
	return variable;
	}

double sqrtSlot::pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const
	{
	//double sumExp=double(numOfBasisFn1+numOfBasisFn2)/2.+1;
	double sumExp=numOfBasisFn1+numOfBasisFn2+1;
	return 1./sumExp*(pow(bounds.getUpperBound(),sumExp)-pow(bounds.getLowerBound(),sumExp));
	}











