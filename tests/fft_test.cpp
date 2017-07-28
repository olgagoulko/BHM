#include "histogram.hpp"
#include "spline.hpp"
#include "basic.hpp"
#include "slot.hpp"

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_complex_math.h>
#include <complex>

using namespace std;  

complex<double> I (0.0,1.0);

int main(int argc, char **argv) {
	
	cout << "----------------- Test distribution with two Gaussian peaks ----------------" << endl;
	
	long unsigned int seed=345902845;
	gsl_rng * RNG = gsl_rng_alloc (gsl_rng_mt19937);
	gsl_rng_set (RNG, seed);
	
	long samplingSteps=1e6;

	vector<slotBounds> intervalBounds;

	double minVar=0; double maxVar=20.0;
	double slotWidth=(maxVar-minVar)/pow(2.,14);
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth);
	histogramBasis binHistogram(histogramVector);
	
	double currentVar=0;
	for(int i=0; i<samplingSteps;i++)
		{
		currentVar = gsl_ran_exponential(RNG, 1.0);
		binHistogram.sampleUniform(currentVar,1);
		}

	histogramBasis scaledBinHistogram = binHistogram.scaledHistogram(samplingSteps);

	splineArray testBHMfit = binHistogram.BHMfit(4, 2, samplingSteps, 2, 0, true);
	testBHMfit.printSplineArrayInfo(); cout << endl;
	testBHMfit.printSplines();
	
	//FFT of a real valued array (result is "half-complex", symmetry f(omega)=f*(-omega))
	int length=binHistogram.getSize();
	int bhmlength=1*length; double bhmslot=maxVar/double(bhmlength);
	double normalisation=maxVar/double(length);	//for transforms in tau
	double bhmnormalisation=maxVar/double(bhmlength);
	double data[length]; double bhmdata[bhmlength];
	for (int i = 0; i < length; i++)
		{
		data[i] = scaledBinHistogram.sampledFunctionValueWeightedAverage(slotWidth*(i+0.5)).first;
		}
	for (int i = 0; i < bhmlength; i++)
		{
		bhmdata[i] = testBHMfit.splineValue(bhmslot*(i+0.5));
		}
		
	gsl_fft_real_radix2_transform (data, 1, length);
	gsl_fft_real_radix2_transform (bhmdata, 1, bhmlength);
		
	//we will output the same information in the conventional order as a complex blitz array (frequencies from 0 to 1/2dt) -- i do this manually because the gsl function is super intransparent
		
	complex<double> endarray[length/2+1]; complex<double> bhmendarray[bhmlength/2+1];// only positive frequencies
	endarray[0]=normalisation*data[0]; bhmendarray[0]=bhmnormalisation*bhmdata[0];
	for (int i = 1; i < length-i; i++)
		{
		endarray[i]=normalisation*(data[i]-I*data[length-i]);	//minus sign here because of different FT sign convention of GSL
		}
	for (int i = 1; i < bhmlength-i; i++)
		{
		bhmendarray[i]=bhmnormalisation*(bhmdata[i]-I*bhmdata[bhmlength-i]);
		}
	endarray[length/2]=normalisation*data[length/2];
	bhmendarray[bhmlength/2]=bhmnormalisation*bhmdata[bhmlength/2];
		
	double freqstep=2.*PI/maxVar;

	
	ofstream output("histogram_testoutput.dat");
// 	double printStep=0.01; pair<double,double> coarseGrainedResult;
// 	for(int i=0; i<(maxVar-minVar)/printStep;i++) 
// 		{
// 		currentVar=minVar+i*printStep+printStep/2;
// 		coarseGrainedResult=scaledBinHistogram.sampledFunctionValueWeightedAverage(currentVar);
// 		output << currentVar << '\t' << coarseGrainedResult.first << '\t' << coarseGrainedResult.second << '\t' 
// 		<< '\t' << testBHMfit.splineValue(currentVar) << '\t' << testBHMfit.splineError(currentVar)
// 		<< endl;
// 		}
	
	for(int i=0; i<length/2+1;i++) 
		{
		output << i*freqstep << '\t' << real(endarray[i]) << '\t' << imag(endarray[i]) << '\t' << real(bhmendarray[i]) << '\t' << imag(bhmendarray[i]) << endl;
		}
	
	
	cout << "---------------------------------- Testing finished -------------------------------------" << endl;
	
	return 0;
	
}