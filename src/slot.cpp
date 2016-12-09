#include "basic.hpp"
#include "slot.hpp"

using namespace std;


excessBin::excessBin()
{
excessCounter=0;
sumExcessValues=0;
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

double slotBounds::slotWidth() const
{
	double width=0;
	if(noUpperBound==false) {width= upperBound-lowerBound;}
	else {cout << "ERROR: infinite slot width " << endl; exit(EXIT_FAILURE);}
	return width;
}


void slotBounds::printBoundsInfo() const
{
cout << fixed << setprecision(16) << "lowerBound = " << lowerBound << endl;
if(noUpperBound) cout << "no upper bound" << endl;
else cout << fixed << setprecision(16) << "upperBound = " << upperBound << endl;
}


//----------------------------------------------------------------------------------------------------------------------------------

basisSlot::basisSlot(slotBounds theBounds, int theTotalNumOfBasisFn)
{
    bounds=theBounds;
    totalNumOfBasisFn=1;
    if(theTotalNumOfBasisFn>1) totalNumOfBasisFn=theTotalNumOfBasisFn;
    GramSchmidtCoeffs.resize(totalNumOfBasisFn);
    numberTimesSampled=0;
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
	if(newnorm<=0) {cout << "ERROR: negative norm in Gram-Schmidt initialization in function number " << i << endl; printSlotInfo(); exit(EXIT_FAILURE);}
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
    for(int j=0; j<totalNumOfBasisFn; j++)
        {
		samplingTotal=weight(variable)*valueToSample*GramSchmidtBasisFn(j, variable);
		delta=samplingTotal-sampledCoeffsValues[j];
		sampledCoeffsValues[j]+=delta/double(numberTimesSampled);
		sampledCoeffsVariance[j]+=delta*(samplingTotal-sampledCoeffsValues[j]);
        }
    }
    else {cout << "ERROR: trying to sample in the wrong slot" << endl; exit(EXIT_FAILURE);}
}

bool basisSlot::addAnotherSlot(basisSlot* anotherSlot)
{
    bool addable=false;
    if( (totalNumOfBasisFn==(anotherSlot -> totalNumOfBasisFn)) && (bounds==(anotherSlot -> bounds)))
    {
	    addable=true; long newNumberTimesSampled=numberTimesSampled+(anotherSlot -> numberTimesSampled);
	for(int j=0; j<totalNumOfBasisFn; j++)
		{
		sampledCoeffsValues[j]*=numberTimesSampled;
		sampledCoeffsValues[j]+=(anotherSlot -> sampledCoeffsValues)[j]*(anotherSlot -> numberTimesSampled);
		sampledCoeffsValues[j]*=1./double(newNumberTimesSampled);
		sampledCoeffsVariance[j]+=(anotherSlot -> sampledCoeffsVariance)[j];
		}
	numberTimesSampled=newNumberTimesSampled;
    }
    return addable;
}

void basisSlot::scale(double norm)
{
	if(norm<VERY_SMALL_NUMBER) {cout << "ERROR: norm in basis slot scaling too small: " << norm << endl; exit(EXIT_FAILURE);}
	for(int j=0; j<totalNumOfBasisFn; j++)
		{
			sampledCoeffsValues[j]*=1./norm;
			sampledCoeffsVariance[j]*=1./norm/norm;
		}
}

double basisSlot::sampledFunctionValue(double variable) const
{
    double result=0;
    for(int i=0; i<totalNumOfBasisFn; i++) result+=sampledCoeffsValues[i]*GramSchmidtBasisFn(i,variable);
    result=result*bounds.slotWidth();
    return result;
}

double basisSlot::sampledFunctionVariance(double variable) const
{
	double result=0;
	if(numberTimesSampled>1)
		{
		for(int i=0; i<totalNumOfBasisFn; i++) {result+=sampledCoeffsVariance[i]*GramSchmidtBasisFn(i,variable)*GramSchmidtBasisFn(i,variable)/double(numberTimesSampled-1);}
		}
	result=result*bounds.slotWidth()*bounds.slotWidth();
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

void basisSlot::printSlotInfo() const
{
bounds.printBoundsInfo();
cout << "totalNumOfBasisFn = " << totalNumOfBasisFn << endl;
}


void basisSlot::printGramSchmidtCoeffs() const
{
    for(int i=0; i<totalNumOfBasisFn; i++)
    {
        for(int j=0; j<totalNumOfBasisFn; j++)
        {
            cout << fixed << setprecision(12) << GramSchmidtCoeffs[i][j] << '\t';
        }
        cout << endl;
    }
}

void basisSlot::printSampledCoeffs() const
{
	for(int i=0; i<totalNumOfBasisFn; i++) cout << sampledCoeffsValues[i] << '\t' << sampledCoeffsVariance[i] << endl;
}

//----------------------------------------------------------------------------------------------------------------------------------

taylorSlot::taylorSlot(slotBounds theBounds, int theTotalNumOfBasisFn) : basisSlot(theBounds, theTotalNumOfBasisFn)
{
	if(theBounds.getIsInfinite()==true) {cout << "ERROR: trying to initialize Taylor expansion slot with open bounds" << endl; exit(EXIT_FAILURE);}
	initializeGramSchmidt();
}

void taylorSlot::initializeGramSchmidt()
{
    
}


double taylorSlot::bareBasisFn(int numOfBasisFn, double variable) const
{       /* 1 x x^2 ... x^n*/
	return pow(variable,numOfBasisFn);
}


double taylorSlot::GramSchmidtBasisFn(int numOfBasisFn, double variable) const
{

	double shiftedvariable=2*variable-(bounds.getUpperBound()+bounds.getLowerBound());
	shiftedvariable=shiftedvariable/bounds.slotWidth();
	
	if(abs(shiftedvariable)>1)
		{
		if(isAround(abs(shiftedvariable),1)) {if(shiftedvariable>1) shiftedvariable=1; else shiftedvariable=-1;}
		else {cout << "ERROR: wrong variable in taylorSlot; variable = " << variable << ", number basis function = " << numOfBasisFn << endl; printSlotInfo(); exit(EXIT_FAILURE);}
		}
	
	return sqrt(2*numOfBasisFn+1)*gsl_sf_legendre_Pl (numOfBasisFn, shiftedvariable)/sqrt(bounds.slotWidth());

}


double taylorSlot::pairwiseIntegral(int numOfBasisFn1, int numOfBasisFn2) const
{
	int sumExp=numOfBasisFn1+numOfBasisFn2+1;
	return (pow(bounds.getUpperBound(),sumExp)-pow(bounds.getLowerBound(),sumExp))/double(sumExp);
}



