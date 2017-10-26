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
#include "basic.hpp"
#include "slot.hpp"
#include "spline.hpp"
#include "histogram.hpp"

#include <cstdio>
#include <cassert>
#include <complex>

using namespace std; 

vector<basisSlot*> generateBasisSlots(double minVar, double maxVar, double slotWidth, int numberOverlaps, int totalNumOfBasisFn)
	{
	if(numberOverlaps<1) numberOverlaps = 1;
	
	if(totalNumOfBasisFn<0) totalNumOfBasisFn=0;
	
	if(maxVar<minVar) {double newMin=maxVar; maxVar=minVar; minVar=newMin;}
	
	double range = maxVar-minVar;
	if(range<VERY_SMALL_NUMBER) {std::cerr << "ERROR: sampling interval too small; min = " << minVar << ", max = " << maxVar << ", range = " << range << endl; throw std::runtime_error("samling interval too small"); }
	
	if((slotWidth>range)||(slotWidth<0)) slotWidth=range/double(numberOverlaps);
	
	double offset = slotWidth/double(numberOverlaps);
	vector<basisSlot*> basisSlotVector;
	
	int counter=0;
	for(int i=1;i<numberOverlaps;i++)
		{
		slotBounds currentBounds(minVar, minVar+i*offset);
		basisSlotVector.push_back(new taylorSlot(currentBounds, totalNumOfBasisFn));
		}
	while(minVar+(counter*offset)+slotWidth<maxVar)
		{
		slotBounds currentBounds(minVar+(counter*offset), minVar+(counter*offset)+slotWidth);
		basisSlotVector.push_back(new taylorSlot(currentBounds, totalNumOfBasisFn));
		counter++;
		if((minVar+(counter*offset)+slotWidth>=maxVar)&&(minVar+(counter*offset)<maxVar)) 
			{
			if(numberOverlaps<2)
				{
				slotBounds finalBounds(minVar+counter*offset, maxVar);
				basisSlotVector.push_back(new taylorSlot(finalBounds, totalNumOfBasisFn));
				}
			else
				{
				for(int i=0;i<numberOverlaps-1;i++)
					{
					if(minVar+(counter+i)*offset<maxVar)
						{
						slotBounds finalBounds(minVar+(counter+i)*offset, maxVar);
						basisSlotVector.push_back(new taylorSlot(finalBounds, totalNumOfBasisFn));
						}
					}
				}
			}
		}
	
	return basisSlotVector;
	}


histogramBasis::histogramBasis(vector<basisSlot*> theBasisSlots) : numberOfInboundsSamples(0)
	{
	valuesOutsideBounds = new excessBin();
        //histogram with zero basis slots is permitted
	if(theBasisSlots.size()==0) {lowerBound=0;upperBound=0;}
        else {
		lowerBound=(theBasisSlots[0] -> getBounds()).getLowerBound();
		upperBound=(theBasisSlots[0] -> getBounds()).getUpperBound();
		}
        noUpperBound=false;
	for(unsigned int i=0;i<theBasisSlots.size();i++)
		{
		basisSlots.push_back(theBasisSlots[i]);
		slotBounds currentBounds=(theBasisSlots[i] -> getBounds());
		if(currentBounds.getLowerBound()<lowerBound) lowerBound=currentBounds.getLowerBound();
		if(currentBounds.getIsInfinite()) noUpperBound=true;
		else if(currentBounds.getUpperBound()>upperBound) upperBound=currentBounds.getUpperBound();
		}
	normalizationFactor=0;	//corresponds to omitting normalization, or norm=1 in this context
	}

histogramBasis::~histogramBasis()
	{
	delete valuesOutsideBounds;
	for (vector<basisSlot*>::iterator i = basisSlots.begin(); i != basisSlots.end(); ++i ) delete *i;
	basisSlots.clear();
	}
	
histogramBasis& histogramBasis::operator=(const histogramBasis& toBeAssigned)
	{
	if (this != &toBeAssigned)
		{
		delete valuesOutsideBounds;
		valuesOutsideBounds = new excessBin(toBeAssigned.getExcessCounter(),toBeAssigned.getExcessValues(1.));

		for (vector<basisSlot*>::iterator i = basisSlots.begin(); i != basisSlots.end(); i++ ) delete *i;
		basisSlots.clear();
		
		if(toBeAssigned.getSize()==0) lowerBound=0;
		else lowerBound=(toBeAssigned.getSlot(0) -> getBounds()).getLowerBound(); 
		upperBound=0; noUpperBound=false;
                numberOfInboundsSamples=toBeAssigned.numberOfInboundsSamples;
		normalizationFactor=toBeAssigned.normalizationFactor;
			
		for (unsigned int i=0;i<toBeAssigned.getSize();i++)
			{
			basisSlot* theSlot=toBeAssigned.getSlot(i) -> Clone();
			basisSlots.push_back(theSlot);
		
			slotBounds currentBounds=(basisSlots[i] -> getBounds());
			if(currentBounds.getLowerBound()<lowerBound) lowerBound=currentBounds.getLowerBound();
			if(currentBounds.getIsInfinite()) noUpperBound=true;
			else if(currentBounds.getUpperBound()>upperBound) upperBound=currentBounds.getUpperBound();
			}
		}
	return *this;
	}
	
histogramBasis::histogramBasis(const histogramBasis& toBeCopied): numberOfInboundsSamples(toBeCopied.numberOfInboundsSamples)
	{
	valuesOutsideBounds = new excessBin(toBeCopied.getExcessCounter(),toBeCopied.getExcessValues(1.));
	
	if(toBeCopied.getSize()==0) lowerBound=0;
	else lowerBound=(toBeCopied.getSlot(0) -> getBounds()).getLowerBound(); 
	upperBound=0; noUpperBound=false;
	
	normalizationFactor=toBeCopied.normalizationFactor;
	
	for (unsigned int i=0;i<toBeCopied.getSize();i++)
		{
		basisSlot* theSlot=toBeCopied.getSlot(i) -> Clone();
		basisSlots.push_back(theSlot);
			
		slotBounds currentBounds=(basisSlots[i] -> getBounds());
		if(currentBounds.getLowerBound()<lowerBound) lowerBound=currentBounds.getLowerBound();
		if(currentBounds.getIsInfinite()) noUpperBound=true;
		else if(currentBounds.getUpperBound()>upperBound) upperBound=currentBounds.getUpperBound();
		}
	}


namespace {
    /// convenience class to hold histogram line
    struct hist_line_s {
        unsigned long nhits;
        double xmin, avg, m2;
        
        hist_line_s(): nhits(), xmin(), avg(1.0), m2(0.0) {}

        /// reads a line with 1..4 numbers from the stream; returns number of values read, or -1 on read error
        int read(std::istream& istrm, std::string& linebuf)
        {
            if (!std::getline(istrm, linebuf)) return -1;
            int count=std::sscanf(linebuf.c_str(), "%lf %lu %lf %lf",
                                  &xmin, &nhits, &avg, &m2);
            return count;
        }
    };
}

histogramBasis::histogramBasis(std::istream& istrm) : valuesOutsideBounds(0),
                                                      lowerBound(0), upperBound(0),
                                                      noUpperBound(false),
                                                      basisSlots(), numberOfInboundsSamples(0), normalizationFactor(0)
{
    try {
        std::string linebuf; // mostly for diagnostics
        hist_line_s histline;

        // first line must contain 2 numbers
        if (histline.read(istrm, linebuf)!=2) {
            throw InvalidFileFormatError("Invalid or empty first line:\n>"+linebuf);
        }
        valuesOutsideBounds=new excessBin(histline.nhits, histline.xmin);	//histline.xmin is a dummy
	normalizationFactor=histline.xmin;
        
        // second line must contain 4 or 2 numbers
	int ndata2=histline.read(istrm, linebuf);
	if ( (ndata2!=4) && (ndata2!=2) ) {
            throw InvalidFileFormatError("Invalid or empty second line:\n>"+linebuf);
        }
        lowerBound=histline.xmin;
        upperBound=histline.xmin;
    
        for ( ; ; ) {
            hist_line_s histline_next; // this is the "next" (or the "last" (truncated)) histogram line
            int ndata=histline_next.read(istrm, linebuf);
            if (ndata==-1) throw InvalidFileFormatError("Unexpected EOF");
            if (ndata==0 || ndata==3 || ndata>4) throw InvalidFileFormatError("Invalid input line:\n>"+linebuf);

            if (histline_next.xmin <= histline.xmin) throw OverlappingSlot("Invalid input line:\n>"+linebuf);

            slotBounds bounds=slotBounds(histline.xmin, histline_next.xmin);
            basisSlot* slot_ptr=new basisSlot(bounds, histline.nhits, histline.avg, histline.m2);
            appendSlot(slot_ptr);
        
            if (ndata==1) {
                break;
            }
            histline=histline_next;

        }
        if (!(numberOfInboundsSamples>0)) throw InvalidFileFormatError("No data points");
    } catch (...) {
        // FIXME: This is a very bad design, but better than leaking memory;
        // FIXME: `valuesOutsideBounds` should not be a raw pointer (or a pointer at all, for that matter!)
        // FIXME: and a proper RAII approach should be used instead.
        delete valuesOutsideBounds;
        throw;
    }
}


void histogramBasis::appendSlot(basisSlot* theSlot)
	{
	basisSlots.push_back(theSlot);
        numberOfInboundsSamples += theSlot->getNumberTimesSampled();
	if((theSlot -> getBounds()).getLowerBound()<lowerBound) lowerBound=(theSlot -> getBounds()).getLowerBound();
	if((theSlot -> getBounds()).getIsInfinite()) noUpperBound=true;
	else if((theSlot -> getBounds()).getUpperBound()>upperBound) upperBound=(theSlot -> getBounds()).getUpperBound();
	}

basisSlot* histogramBasis::combinedSlot(unsigned int startPoint, unsigned int endPoint) const
	{
            if((startPoint>endPoint) || (endPoint>=basisSlots.size())) {std::cerr << "ERROR in combinedSlot" << endl; throw std::runtime_error("ERROR in combinedSlot"); }
	else if(startPoint>endPoint) {unsigned int save=endPoint; endPoint=startPoint; startPoint=save;}
	
	for(unsigned int i=startPoint;i<=endPoint;i++)
		{
		for(unsigned int j=i+1;j<=endPoint;j++)
			{
			if((basisSlots[i] -> getBounds()).overlapping((basisSlots[j] -> getBounds()))) 
				{
				cerr << "ERROR: trying to combine overlapping slots " << i << " and " << j << endl; 
				basisSlots[i] -> printSlotInfo(cerr); basisSlots[j] -> printSlotInfo(cerr);
				throw std::logic_error("ERROR: trying to combine overlapping slots");
				}
			}
		}
	
	slotBounds theBounds( basisSlots[startPoint] -> getBounds().getLowerBound(), basisSlots[endPoint] -> getBounds().getUpperBound());
	basisSlot * combined = new basisSlot(theBounds);
	for(unsigned int i=startPoint;i<=endPoint;i++)
		{
		combined -> combineWithSlot(basisSlots[i]);
		}

	return combined;
	}

basisSlot* histogramBasis::getSlot(unsigned int whichSlot) const
	{
            if(whichSlot>=basisSlots.size()) {cerr << "ERROR in getSlot, " << whichSlot << " does not exist" << endl; throw std::runtime_error("Error in getSlot");}
	return basisSlots[whichSlot];
	}


void histogramBasis::sample(double variable, double valueToSample)
	{
	bool isOutsideBounds=true;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		    if( basisSlots[i] -> checkIfInBasisSlot(variable) )
                        {
                            basisSlots[i] -> sample(variable,valueToSample); isOutsideBounds=false;
                        }
		}
	if(isOutsideBounds)
            {
                valuesOutsideBounds -> sample(variable, valueToSample);
            } else {
                ++numberOfInboundsSamples;
            }
	}

void histogramBasis::sampleUniform(double variable, double valueToSample)
	{
	//faster sampling, assuming all slots are of equal width, without overlaps
	if(variable > upperBound || variable < lowerBound || basisSlots.size()==0) valuesOutsideBounds -> sample(variable, valueToSample);
	else 
		{
		basisSlots[int((variable-lowerBound)/((basisSlots[0] -> getBounds()).slotWidth()))] -> sample(variable,valueToSample);
                ++numberOfInboundsSamples;
		}
	}


histogramBasis histogramBasis::coarseGrainedHistogram(unsigned int minNumberTimesSampled)
	{
	vector<basisSlot*> coarseGrainedSlots;
	unsigned int counter=0; bool enough;
	for(unsigned int i=0;i<basisSlots.size();i++) basisSlots[i] -> updateEnoughSampled(minNumberTimesSampled);
	while(counter<basisSlots.size())
		{
		enough = basisSlots[counter] -> enoughSampled();
		if(enough)
			{
			basisSlot* currentSlot = basisSlots[counter] -> Clone();
			coarseGrainedSlots.push_back(currentSlot);
			counter++;
			}
		else {
			for(unsigned int more=counter+1;more<basisSlots.size();more++)
				{
				basisSlot* combined = combinedSlot(counter, more);
				combined -> updateEnoughSampled(minNumberTimesSampled);
				enough = combined -> enoughSampled();
				if(enough) {coarseGrainedSlots.push_back(combined); counter=more+1; break;}
				}
			if(!enough && (coarseGrainedSlots.size()>0)) 
				{
				slotBounds theBounds( (coarseGrainedSlots.back() -> getBounds()).getLowerBound(), upperBound);
				basisSlot * combined = new basisSlot(theBounds);
				combined -> combineWithSlot(coarseGrainedSlots.back());
				for(unsigned int i=counter;i<basisSlots.size();i++) combined -> combineWithSlot(basisSlots[i]);
				coarseGrainedSlots.back() = combined;
				counter=basisSlots.size();
				}
			}
		}
	
	histogramBasis result(coarseGrainedSlots);
	return result;
	}


pair<double,double> histogramBasis::sampledFunctionValueAverage(double variable) const
	{
	vector<double> sampledValue;
	int numberSlots=0;
	double average=0;
	double error=0;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) ) 
			{
			sampledValue.push_back(basisSlots[i] -> sampledFunctionValue(variable));
			average+=sampledValue[numberSlots];
			numberSlots++;
			}
		}
	
	if(numberSlots==0) average=valuesOutsideBounds -> getExcessValues();
	else if(numberSlots>1)
		{
		average=average/double(numberSlots);
		double roundofferror=0;
		double delta;
		for(int i=0;i<numberSlots;i++)
			{
			delta=sampledValue[i]-average;
			error+=delta*delta;
			roundofferror+=delta;
			}
		roundofferror=roundofferror*roundofferror/double(numberSlots);
		error=error-roundofferror;
		if(error<0) error=0;
		error=sqrt(error/double(numberSlots*(numberSlots-1)));
		}
		
	pair<double,double> thisPair(average,error);
	return thisPair;
	}

pair<double,double> histogramBasis::sampledFunctionValueWeightedAverage(double variable) const
	{
	double sampledValue=0; double value;
	double sumVariance=0; double variance;
	double error;
	int numberSlots=0;
	for(unsigned int i=0;i<basisSlots.size();i++)
	{
		if( basisSlots[i] -> checkIfInBasisSlot(variable) )
			{
			value = basisSlots[i] -> sampledFunctionValue(variable);
			variance = basisSlots[i] -> sampledFunctionError(variable);
			if(variance>0)
				{
				variance=variance*variance;
				sumVariance+=1./variance;
				sampledValue += value/variance;
				}
			numberSlots++;
			}
	}
	
	if(numberSlots==0) sampledValue=valuesOutsideBounds -> getExcessValues();
	else if(sumVariance==0) {error=0;}
	else {sampledValue=sampledValue/sumVariance; error=sqrt(1./sumVariance);}
	
	pair<double,double> result(sampledValue,error);
	
	return result;
	}

histogramBasis histogramBasis::scaledHistogram(long norm)
	{
	vector<basisSlot*> scaledSlots;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		basisSlot* currentSlot = basisSlots[i] -> Clone();
		currentSlot -> scale(norm);
		scaledSlots.push_back(currentSlot);
		}
	
	histogramBasis result(scaledSlots); //FIXME: do these histograms have the right parameters otherwise?
	return result;
	}

histogramBasis histogramBasis::normalizedHistogram(double norm)
	{
	vector<basisSlot*> scaledSlots;
	for(unsigned int i=0;i<basisSlots.size();i++)
		{
		basisSlot* currentSlot = basisSlots[i] -> Clone();
		if( (norm>0) && (norm!=1) ) currentSlot -> normalize(norm);
		scaledSlots.push_back(currentSlot);
		}
	
	histogramBasis result(scaledSlots);//FIXME: do these histograms have the right parameters otherwise?
	return result;
	}

bool histogramBasis::addAnotherHistogram(histogramBasis anotherHistogram)
	{
	bool addable=false;
	if(basisSlots.size()==anotherHistogram.getSize())
		{
		addable=true;
		for(unsigned int i=0;i<basisSlots.size();i++)
			{
			addable = basisSlots[i] -> addAnotherSlot(anotherHistogram.getSlot(i));
			}
		}
	return addable;
	}


vector< vector< basisSlot* > > histogramBasis::binHierarchy(long norm, unsigned int dataPointsMin, double usableBinFraction)
	{
	unsigned int numberElementaryBins = basisSlots.size();
	unsigned int maxLevel=rounding(log(double(numberElementaryBins))/log(2));
	vector< vector<basisSlot*> > analysisBins;
	vector<basisSlot*> currentLevelBins;
	for(unsigned int i=0;i<numberElementaryBins;i++)
		{
		basisSlots[i] -> updateEnoughSampled(dataPointsMin);
		currentLevelBins.push_back(basisSlots[i] -> Clone());
		}
	for(unsigned int j=0;j<=maxLevel;j++)
		{
		analysisBins.insert(analysisBins.begin(),currentLevelBins);
		unsigned int currentSize = currentLevelBins.size();
		if(currentSize==1) break;
		currentLevelBins.resize(0);
		for(unsigned int i=0;i<currentSize;i+=2)
			{
			slotBounds theBounds((analysisBins[0][i] -> getBounds()).getLowerBound(), (analysisBins[0][i+1] -> getBounds()).getUpperBound());
			basisSlot * combined = new basisSlot(theBounds);
			combined -> combineWithSlot(analysisBins[0][i]);
			combined -> combineWithSlot(analysisBins[0][i+1]);
			combined -> updateEnoughSampled(dataPointsMin);
			currentLevelBins.push_back(combined);
			}
		}
	
	//only use slots where enough sampled, disregard levels if too few usable slots
	unsigned int breakPoint=analysisBins.size();
	for(unsigned int j=0;j<analysisBins.size();j++)
		{
		int initialAnalysisBinsLevelSize=analysisBins[j].size()-1;
		for(int i=initialAnalysisBinsLevelSize;i>=0;i--)
			{
			if(analysisBins[j][i] -> enoughSampled()) analysisBins[j][i] -> scale(norm);
			else analysisBins[j].erase(analysisBins[j].begin()+i);
			}
		if( (analysisBins[j].size()<initialAnalysisBinsLevelSize*usableBinFraction) || (analysisBins[j].size()==0))
			{
			breakPoint=j; break;
			}
		}
	unsigned int initialAnalysisBinsSize=analysisBins.size();
	for(unsigned int j=breakPoint;j<initialAnalysisBinsSize;j++) analysisBins.pop_back();
	
	return analysisBins;
	}


splineArray histogramBasis::BHMfit(BHMparameters parameters, long norm, bool fail_if_zero)
	{
	unsigned int splineOrder=parameters.splineOrder;
	unsigned int minLevel=parameters.minLevel;
	fitAcceptanceThreshold theThreshold=parameters.threshold;
	double usableBinFraction=parameters.usableBinFraction;
	double jumpSuppression=parameters.jumpSuppression;
		
	if(splineOrder<1) {
            // rationale: it's caller's job to check parameters and issue warning to an appropriate log
            // here we simply bail out if the input is wrong.
            throw std::invalid_argument("BHMfit(): splineOrder must be at least 1");
        }
		
	vector<double> aMaxVector; vector<double> chisqArray;
	
	if(minLevel<2) minLevel=2;
        int maxLevel=ilog2(basisSlots.size());
        if (maxLevel<0) {
            // rationale: it's caller's responsibility to check this precondition
            throw std::invalid_argument("Number of elementary bins is not a power of 2");
        }
	
	//make bin hierarchy; position in outer vector denotes bin level: 0 is largest bin, 1 are second level bins etc
	//intervals labeled in the same way: as a sequence 0 or 11 or 221 or 2331 etc
	vector< vector<basisSlot*> > analysisBins=binHierarchy(norm, parameters.dataPointsMin, usableBinFraction); //the bins on each given level don't have to be in order, only analysisBins do

	if(analysisBins.size()<minLevel) {
            // rationale: hard to be checked by caller, but no meaningful spline can be created
            throw NotEnoughData_Error();
        }

	if(isDataConsistentWithZero(analysisBins)==true) {
            LOGGER << "WARNING: Data is consistent with zero on the interval";
            // rationale: we have to throw here because the return object is not going to be (meaningfully) constructed
            if (fail_if_zero) throw ConsistentWithZero_Error();
        }
	
	bool allSplinesGood; bool checkIntervals;
	vector<slotBounds> intervalBounds;
	vector<unsigned int> intervalOrders;
	double thresholdIncrease=0;
	if(theThreshold.max<=theThreshold.min) theThreshold.steps=0;
	if(theThreshold.steps>0) thresholdIncrease=(theThreshold.max-theThreshold.min)/double(theThreshold.steps);
	double currentFitAcceptanceThreshold;
	for(int thresholdStep=0; thresholdStep<=theThreshold.steps; thresholdStep++)
		{
		currentFitAcceptanceThreshold=theThreshold.min+thresholdStep*thresholdIncrease;
		LOGGER << "Begin BHM fitting with threshold T = " << currentFitAcceptanceThreshold;
		
		intervalBounds.resize(0);
		intervalOrders.resize(0);
		vector<unsigned int> intervalNumbers;
		slotBounds histogramBounds(lowerBound,upperBound);
		intervalBounds.push_back(histogramBounds);
		unsigned int currentNumberIntervals;
		intervalOrders.push_back(0); intervalNumbers.push_back(0);
		unsigned int currentLevel=minLevel;

		while(currentLevel<analysisBins.size()-1)
			{
			splineArray result = matchedSplineFit(analysisBins, intervalBounds, splineOrder, 0, aMaxVector, chisqArray, currentFitAcceptanceThreshold);

			if(result.checkOverallAcceptance(currentFitAcceptanceThreshold)!=true) checkIntervals=true;
			else checkIntervals=false;
			
			currentNumberIntervals=intervalBounds.size();
			bool isIntervalGood[currentNumberIntervals];
			allSplinesGood=true;
			chisqArray.resize(0);
			for(unsigned int i=0;i<currentNumberIntervals;i++)
				{
				if(checkIntervals) LOGGER << "Checking interval " << i << " (order: " << intervalOrders[i] << ", number: " << intervalNumbers[i] << ")";
				double chisqArrayElement = 1+currentFitAcceptanceThreshold*sqrt(2.);
				bool currentSplineGood=result.getSplinePiece(i) -> checkIntervalAcceptance(analysisBins, currentFitAcceptanceThreshold, chisqArrayElement, intervalOrders[i], checkIntervals);
				chisqArray.push_back(chisqArrayElement);
			
				if(checkIntervals)
					{
					if(currentSplineGood==false) allSplinesGood=false;
					isIntervalGood[i]=currentSplineGood;
					}
				}
					
			if(checkIntervals)
				{
				unsigned int intervalCounter=0;
				for(unsigned int i=0;i<currentNumberIntervals;i++)
					{
					if(!isIntervalGood[i])
						{
						intervalOrders[intervalCounter]++;
						intervalOrders.insert(intervalOrders.begin()+intervalCounter+1,intervalOrders[intervalCounter]);
						intervalNumbers[intervalCounter]*=2;
						intervalNumbers.insert(intervalNumbers.begin()+intervalCounter+1,intervalNumbers[intervalCounter]+1);
						intervalCounter++;
						}
					intervalCounter++;
					}
				
				intervalBounds.resize(0);
				for(unsigned int i=0;i<intervalOrders.size();i++)
					{
					slotBounds currentBounds( (basisSlots[pow(2,maxLevel-intervalOrders[i])*intervalNumbers[i]] -> getBounds()).getLowerBound(), (basisSlots[pow(2,maxLevel-intervalOrders[i])*(intervalNumbers[i]+1)-1] -> getBounds()).getUpperBound());
				
					intervalBounds.push_back(currentBounds);
					}
				}
			
			if(allSplinesGood)
				{
				LOGGER << "Good spline found with threshold T = " << currentFitAcceptanceThreshold;
				break;
				}
			currentLevel++;
			}
		if(allSplinesGood) break;
		}

	if(!allSplinesGood)
                {
                    jumpSuppression=0;
                    LOGGER << "No acceptable fit could be found with the current threshold "
                           << currentFitAcceptanceThreshold;
                }
	else if(intervalBounds.size()==1)
                {
                    jumpSuppression=0; //current setup is not to constrain highest derivate at domain boundaries, only at spline knots
                    LOGGER << "No knots for jump suppression";
                } 
	
	if(jumpSuppression>0)
		{
		LOGGER << "\nIterative consistent constraints procedure\n";
		aMaxVector.resize(0);
		splineArray geta3=matchedSplineFit(analysisBins, intervalBounds, splineOrder, 0, aMaxVector, chisqArray, currentFitAcceptanceThreshold);
		for(unsigned int i=0;i<intervalOrders.size();i++) aMaxVector.push_back((geta3.getSplinePiece(i) -> getCoefficients())[splineOrder-1]);
		}
	
	splineArray result = matchedSplineFit(analysisBins, intervalBounds, splineOrder, jumpSuppression, aMaxVector, chisqArray, currentFitAcceptanceThreshold);
	
	if(jumpSuppression>0)
		{
		for(int repeat=0;repeat<1;repeat++)
			{
			bool CCprocedureConverged=false; int iterations=0; double maximalGoodGluingFactor=0; double minimalBadGluingFactor=1/VERY_SMALL_NUMBER; int maxIterations=20;
			while(CCprocedureConverged==false)
				{
				iterations++;
				allSplinesGood = result.checkOverallAcceptance(currentFitAcceptanceThreshold);
				LOGGER << "Global suppression factor "
                                       << jumpSuppression
                                       << "; the interval fit is "
                                       << (allSplinesGood? "good" : "not good");
				
				if(!allSplinesGood)
					{
					if(jumpSuppression<minimalBadGluingFactor) minimalBadGluingFactor=jumpSuppression;
					
					if(iterations>maxIterations) {jumpSuppression=maximalGoodGluingFactor; CCprocedureConverged=true; allSplinesGood=true;}
					else if(maximalGoodGluingFactor>0) {jumpSuppression=(maximalGoodGluingFactor+minimalBadGluingFactor)/2.;}
					else jumpSuppression*=0.1;
					}
				else	{
					if(jumpSuppression>maximalGoodGluingFactor) maximalGoodGluingFactor=jumpSuppression;
					
					if(iterations>maxIterations) {CCprocedureConverged=true; jumpSuppression=maximalGoodGluingFactor;}
					else if(minimalBadGluingFactor<1/VERY_SMALL_NUMBER) {jumpSuppression=(maximalGoodGluingFactor+minimalBadGluingFactor)/2.;}
					else jumpSuppression*=2;
					}
				
				result = matchedSplineFit(analysisBins, intervalBounds, splineOrder, jumpSuppression, aMaxVector, chisqArray, currentFitAcceptanceThreshold);
				}
			aMaxVector.resize(0); chisqArray.resize(0);
			for(unsigned int i=0;i<intervalOrders.size();i++) aMaxVector.push_back((result.getSplinePiece(i) -> getCoefficients())[splineOrder-1]);
			for(unsigned int i=0;i<intervalBounds.size();i++)
				{
				double chisqArrayElement=1+currentFitAcceptanceThreshold*sqrt(2.);
				result.getSplinePiece(i) -> checkIntervalAcceptance(analysisBins, currentFitAcceptanceThreshold, chisqArrayElement, intervalOrders[i], false);
				chisqArray.push_back(chisqArrayElement);
				}
			jumpSuppression=1;
			}
		}
	
	vector< vector<basisSlot*> >::iterator i;
	vector<basisSlot*>::iterator j;
	for (i = analysisBins.begin(); i != analysisBins.end(); i++) for (j = i->begin(); j != i->end(); j++) delete *j;
	analysisBins.clear();

	result.updateGoodness(allSplinesGood, currentFitAcceptanceThreshold);
	return result;
	}


std::ostream& operator<<(std::ostream& ostrm, const histogramBasis& hist)
{
    ostrm << hist.getExcessCounter() << "\n";
    
    for(unsigned int i=0; i<hist.getSize(); ++i) 
    {
        basisSlot* slot = hist.getSlot(i);
        ostrm << slot->getBounds().getLowerBound() << '\t'
              << slot->getNumberTimesSampled() << '\t'
              << slot->sampledIntegral() << '\t'
              << slot->getVariance()
              << '\n';
    }
    if (hist.getSize()!=0) { 
        ostrm << hist.getSlot(hist.getSize()-1)->getBounds().getUpperBound()
              << endl;
    }
    return ostrm;
}



splineArray matchedSplineFit(vector< vector<basisSlot*> > currentAnalysisBins, vector< slotBounds > intervalBounds, unsigned int splineOrder, double jumpSuppression, vector<double> aMaxVector, vector<double> chisqArray, double currentFitAcceptanceThreshold)
	{
            if(intervalBounds.size()==0) { cerr << "ERROR in matchedSplineFit, intervalBounds size is zero" << endl; throw std::runtime_error("ERROR in matchedSplineFit"); }
	if(splineOrder<1) {LOGGER << "WARNING: splineOrder has to be at least 0, setting splineOrder to 0"; splineOrder=1;} //splineOrder=m+1
	
	unsigned int numberIntervals=intervalBounds.size();
	unsigned int matrixRows = 0;
	for(unsigned int i=0;i<currentAnalysisBins.size();i++) matrixRows+=currentAnalysisBins[i].size();
	unsigned int matrixCols=numberIntervals+splineOrder-1;
	
	vector< vector<unsigned int> > binsFullyInsideInterval(numberIntervals);
	
/**/ 	unsigned int matrixRows2=matrixRows;
	if(jumpSuppression>0)
		{
		matrixRows2+=numberIntervals-1;
		//replace by line below if also want to constrain highest derivative on domain boundaries
		//matrixRows2+=numberIntervals; if(numberIntervals>1) matrixRows2++;
		}
/**/

	gsl_matrix * fullDesignMatrix = gsl_matrix_calloc (matrixRows2, matrixCols);
	double* b = integralVector(currentAnalysisBins);
	gsl_vector * gslb= gsl_vector_calloc (matrixRows2);
	for(unsigned int i=0;i<matrixRows;i++) gsl_vector_set(gslb,i,b[i]);

	double* binVec=binomialVector(splineOrder);
	double intervalMatrix[numberIntervals-1][splineOrder-1];
	double lastIntervalMatrix[numberIntervals-1][splineOrder-1];
	
	double* basicDesignMatrix = designMatrix(currentAnalysisBins, splineOrder);
	double currentLowerBound=0; double currentUpperBound;
	for(unsigned int i=0;i<numberIntervals-1;i++)
		{
		currentUpperBound = intervalBounds[i].getUpperBound();
		
		for(unsigned int k=0; k < splineOrder-1;k++)
			{
			intervalMatrix[i][k]=(pow(currentUpperBound,splineOrder-1-k)-pow(currentLowerBound,splineOrder-1-k))*binVec[k]*pow(-1,splineOrder-k);
			lastIntervalMatrix[i][k]=-pow(currentUpperBound,splineOrder-1-k)*binVec[k]*pow(-1,splineOrder-k);
			}
		
		currentLowerBound=currentUpperBound;
		}

	//construction of full desgin matrix
	unsigned int currentMatrixRow=0;
	for(unsigned int j1=0;j1 < currentAnalysisBins.size();j1++)
		for(unsigned int j2=0;j2 < currentAnalysisBins[j1].size();j2++)
			{
			slotBounds currentBounds=(currentAnalysisBins[j1][j2] -> getBounds());
			double currentError=currentAnalysisBins[j1][j2] -> sampledIntegralError();
			double weightFactor=sqrt(double(pow(2,j1)));
			unsigned int start; unsigned int end;
			for(unsigned int i=0;i<numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {start=i; break;}
				}
			for(unsigned int i=start;i<numberIntervals;i++)
				{
				if(currentBounds.overlapping(intervalBounds[i])) {end=i;}
				}
		
			//regular elements for a0, a1, a2
			for(unsigned int k=0; k < splineOrder-1;k++)
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,k,basicDesignMatrix[currentMatrixRow*splineOrder+k]);
		
			//very last element (x_i+1^4-x_i^4)/4
			if(start==end)
				{
				gsl_matrix_set(fullDesignMatrix,currentMatrixRow,start+splineOrder-1,basicDesignMatrix[currentMatrixRow*splineOrder+splineOrder-1]);
				binsFullyInsideInterval[start].push_back(currentMatrixRow);
				}
			else
				{
				for(unsigned int i=start;i<=end;i++)
					{
					if( currentError < VERY_SMALL_NUMBER) gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1,0);
					else gsl_matrix_set(fullDesignMatrix,currentMatrixRow,i+splineOrder-1, splineBasisFunction(intervalBounds[i],splineOrder-1)/currentError/weightFactor);
					}
				}
		
			//other a3 contributions
			if((start==end)&&(start>0))
				{
				for(unsigned int kk=0; kk < splineOrder-1;kk++)
					{
					for(unsigned int k=0; k < start; k++)
						gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]);
				
					gsl_matrix_set(fullDesignMatrix, currentMatrixRow, start+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,start+splineOrder-1)+lastIntervalMatrix[start-1][kk]*basicDesignMatrix[currentMatrixRow*splineOrder+kk]);
					}
				}
			else if(start!=end)
				{
				if(start==0) start++;
				for(unsigned int i=start;i<=end;i++)
					{
					if( currentError > VERY_SMALL_NUMBER)
						{
						for(unsigned int kk=0; kk < splineOrder-1;kk++)
							{
							for(unsigned int k=0; k < i; k++)
								gsl_matrix_set(fullDesignMatrix, currentMatrixRow, k+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,k+splineOrder-1)+intervalMatrix[k][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/weightFactor);
						
							gsl_matrix_set(fullDesignMatrix, currentMatrixRow, i+splineOrder-1, gsl_matrix_get(fullDesignMatrix,currentMatrixRow,i+splineOrder-1)+lastIntervalMatrix[i-1][kk]*splineBasisFunction(intervalBounds[i],kk)/currentError/weightFactor);
							}
						}
					}
				}
				
			gsl_vector_set(gslb,currentMatrixRow,gsl_vector_get(gslb,currentMatrixRow));
			currentMatrixRow++;
			}
			
	if(jumpSuppression>0)
	{
	double factor=sqrt(jumpSuppression/double(numberIntervals));
	for(unsigned int i=0;i<numberIntervals-1;i++) 
		{
		double minChisq=min(chisqArray[i+1],chisqArray[i]); if(i>0) minChisq=min(minChisq,chisqArray[i-1]); if(i+2<numberIntervals) minChisq=min(minChisq,chisqArray[i+2]);
		double setToWhat=factor*minChisq/(aMaxVector[i+1]-aMaxVector[i]);
		gsl_matrix_set(fullDesignMatrix, matrixRows+i, splineOrder-1+i, -setToWhat);
		gsl_matrix_set(fullDesignMatrix, matrixRows+i, splineOrder+i, setToWhat);
		}
	//uncomment below if also want to constrain highest derivative on domain boundaries
	//gsl_matrix_set(fullDesignMatrix, matrixRows+numberIntervals-1, splineOrder-1, factor*chisqArray[0]/(aMaxVector[0]));
	//if(numberIntervals>1) gsl_matrix_set(fullDesignMatrix, matrixRows+numberIntervals, splineOrder+numberIntervals-2, factor*chisqArray[numberIntervals-1]/(aMaxVector[numberIntervals-1]));
	}
	
	//chi^2 minimization using full design matrix
	gsl_matrix * V = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_vector * diagonal = gsl_vector_alloc (matrixCols);
	gsl_vector *x = gsl_vector_alloc (matrixCols);
	solveSVD(fullDesignMatrix, V, gslb, diagonal, x);
	vector<splinePiece*> splines;

	currentMatrixRow=0;
	vector<double> theChisq; vector<int> theDOF;
	LOGGER << "Checking separate chi_n^2/n in spline fit\n"
		<< left << setw(8) << "level" << setw(8) << "n" 
		<< setw(16) << "chi_n^2/n" << setw(16) << "max chi_n^2/n";
	for(unsigned int j=0;j < currentAnalysisBins.size();j++)
		{
		int dof = currentAnalysisBins[j].size();
		gsl_vector_view currentgslb = gsl_vector_subvector (gslb, currentMatrixRow, dof);
		double chisq; gsl_blas_ddot (&currentgslb.vector, &currentgslb.vector, &chisq);
		chisq*=double(pow(2,j))/double(dof);
		theChisq.push_back(chisq); theDOF.push_back(dof);
	
		LOGGER << left << fixed << setprecision(4) << setw(8) << j << setw(8) << dof 
			<< right << setw(9) << chisq << setw(7) << " "
			<< left << setw(16) << 1+currentFitAcceptanceThreshold*sqrt(2./double(dof));
	
		currentMatrixRow+=dof;
		}

	double diagScaling[matrixCols];
	for(unsigned int j=0;j < matrixCols;j++) 
		{
		diagScaling[j]=0;
		if( gsl_vector_get(diagonal,j)>VERY_SMALL_NUMBER ) diagScaling[j]=1./gsl_vector_get(diagonal,j);
		}
	
	gsl_matrix * Vscaled = gsl_matrix_alloc (matrixCols, matrixCols);
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			gsl_matrix_set(Vscaled,i,j,gsl_matrix_get(V,i,j)*diagScaling[j]);
		
	gsl_matrix * fullCovarianceMatrix = gsl_matrix_alloc (matrixCols, matrixCols);
	gsl_blas_dsyrk (CblasUpper, CblasNoTrans, 1, Vscaled, 0, fullCovarianceMatrix);
	double covarianceMatrix[matrixCols][matrixCols];
	for (unsigned int i = 0; i < matrixCols; i++)
		for (unsigned int j = 0; j < matrixCols; j++)
			covarianceMatrix[i][j]= gsl_matrix_get (fullCovarianceMatrix, i,j);

	//updating coefficients and error coefficients of all spline pieces
	vector<double> theCoefficients;
	for(unsigned int n=0;n<numberIntervals;n++)
		{
		if(n==0) 
			{
			for(unsigned int i=0;i<splineOrder;i++) theCoefficients.push_back(gsl_vector_get(x,i));
			}
		else
			{
			for(unsigned int i=0;i<splineOrder-1;i++) theCoefficients[i]+=(gsl_vector_get(x,n+splineOrder-1)-theCoefficients[splineOrder-1])*lastIntervalMatrix[n-1][i];
			theCoefficients[splineOrder-1]=gsl_vector_get(x,n+splineOrder-1);
			}

		vector<double> theErrorCoefficients;
		for(unsigned int j=0;j < 2*splineOrder-1;j++) theErrorCoefficients.push_back(0); 

		if(n>0)
			{
			for(unsigned int j=0;j < splineOrder-1;j++)
				{
				for(unsigned int k=j;k < splineOrder-1;k++)
					{
					covarianceMatrix[j][k]+=(covarianceMatrix[j][splineOrder+n-1]-covarianceMatrix[j][splineOrder-1])*lastIntervalMatrix[n-1][k];
					covarianceMatrix[j][k]+=(covarianceMatrix[k][splineOrder+n-1]-covarianceMatrix[k][splineOrder-1])*lastIntervalMatrix[n-1][j];
					covarianceMatrix[j][k]+=(covarianceMatrix[splineOrder-1][splineOrder-1] + covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - 2*covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][k]*lastIntervalMatrix[n-1][j];
					}
				covarianceMatrix[j][splineOrder-1]=covarianceMatrix[j][splineOrder+n-1];
				covarianceMatrix[j][splineOrder-1]+=(covarianceMatrix[splineOrder+n-1][splineOrder+n-1] - covarianceMatrix[splineOrder+n-2][splineOrder+n-1])*lastIntervalMatrix[n-1][j];
				
				for(unsigned int k=n;k < matrixCols-splineOrder;k++)
					covarianceMatrix[j][k+splineOrder]+=(covarianceMatrix[splineOrder+n-1][k+splineOrder] - covarianceMatrix[splineOrder+n-2][k+splineOrder])*lastIntervalMatrix[n-1][j];
				}
			covarianceMatrix[splineOrder-1][splineOrder-1]=covarianceMatrix[n+splineOrder-1][n+splineOrder-1];
			}

		for(unsigned int j=0;j < splineOrder;j++)
			for(unsigned int k=j;k < splineOrder;k++)
				{
				theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				if(j!=k) theErrorCoefficients[j+k]+=covarianceMatrix[j][k];
				}

		splinePiece * solution = new splinePiece(intervalBounds[n], splineOrder); solution -> setSplinePiece(theCoefficients,theErrorCoefficients);
		splines.push_back(solution);
		}

	splineArray result(splines);
	result.updateLevelProperties(theChisq,theDOF);
	
	gsl_matrix_free(fullDesignMatrix); gsl_vector_free(gslb);
	gsl_matrix_free(V); gsl_vector_free(diagonal); gsl_vector_free(x); gsl_matrix_free(Vscaled); gsl_matrix_free(fullCovarianceMatrix);
	delete [] binVec; delete [] b; delete [] basicDesignMatrix;
	
	return result;
	
}



