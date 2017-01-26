#include "histogram.hpp"
#include "basic.hpp"
#include "slot.hpp"
#include "sput.h"

using namespace std; 

static void test_histogram()
{ 
	double minVar=0; double maxVar=5; double slotWidth=0.5; int numberOverlaps=5; int totalNumOfBasisFn=4;
	vector<basisSlot*> histogramVector=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram(histogramVector);
	
	myTestHistogram.sample(0, 2.0);
	myTestHistogram.sample(0.55, 2.0);
	myTestHistogram.sample(4.0, 5.0);
	myTestHistogram.sample(10, 1.5);
	myTestHistogram.sample(8, 3.0);
	
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(-1.).first  == 4.5, "4.5 in excess bin");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(6.).first  == 4.5, "4.5 in excess bin");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(10.).first  == 4.5, "4.5 in excess bin");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0).first, 32.), "value at zero");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(2.5).first == 0, "0 at 2.5 where not sampled");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(4.88).first == 0, "0 at 4.88 where not sampled");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(4.999999999999).first == 0, "0 just below 5.0 where not sampled");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.55).first,8321./1250.), "value at 0.55");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.52).first,1028866./156250.), "value at 0.52");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.2).first,-48641./84375.), "value at 0.2");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.23).first,-87190663./675000000.), "value at 0.23");
	
	//due to numerical precision (last digit) testing at slot boundaries is impossible --> this will never matter in real life
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.900000000000001).first,7168./1250.), "value at 3.9");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.95).first,12352./1250.), "value at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4).first,6816./250.), "value at 4");
	
	myTestHistogram.sample(4,-0.4);
	myTestHistogram.sample(4.23,1.1);

	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.65).first,-4209./31250./2.), "value at 3.65");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.95).first,28932827./3906250./2.), "value at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.95).second,sqrt(4893669436322603./2.)/23437500./2.), "error at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueWeightedAverage(3.95).first,8.504782282449135/2.), "weighted average value at 3.95"); 
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueWeightedAverage(3.95).second,3.310039266235157/2.), "weighted average error at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4).first,70335263./3906250./2.), "value at 4");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.16).first,189704957./244140625./2.), "value at 4.16");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.23).first,3577765089./976562500./2.), "value at 4.23");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.300000000000001).first,11706064./5859375./2.), "value at 4.3");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.44).first,-97984339./1464843750./2.), "value at 4.44");
	
	numberOverlaps=1; totalNumOfBasisFn=2;
	vector<basisSlot*> histogramVector2=generateBasisSlots(minVar, maxVar, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram2(histogramVector2);
	myTestHistogram2.sample(0.61, 3.0);
	sput_fail_unless(myTestHistogram2.sampledFunctionValueAverage(10.).first  == 0, "0 in excess bin");
	sput_fail_unless(myTestHistogram2.sampledFunctionValueAverage(0.4).first == 0, "0 at 0.4 where not sampled");
	sput_fail_unless(myTestHistogram2.sampledFunctionValueAverage(1.1).first == 0, "0 at 1.1 where not sampled");
	sput_fail_unless(myTestHistogram2.sampledFunctionValueAverage(1.1).second == 0, "error 0 where not sampled");
	sput_fail_unless(isAround(myTestHistogram2.sampledFunctionValueAverage(0.5).first,402./50.), "value at 0.5");
	sput_fail_unless(isAround(myTestHistogram2.sampledFunctionValueAverage(0.9).first,-6./250.), "value at 0.9");

}

static void test_combinedSlot()
{ 
	double lowerBound=1.; double upperBound=2.5; double middleBound=1.2;
	
	slotBounds bounds1(lowerBound, upperBound);
	int myTestNumberOfBasisFn=4;
	
	vector<basisSlot*> myBasisSlots; myBasisSlots.push_back(new taylorSlot(bounds1, myTestNumberOfBasisFn));
	histogramBasis myTestHistogram(myBasisSlots);
	
	slotBounds bounds2(lowerBound, middleBound);
	myTestNumberOfBasisFn=2;
	myTestHistogram.appendSlot(new taylorSlot(bounds2, myTestNumberOfBasisFn));
	
	slotBounds bounds3(middleBound, upperBound);
	myTestNumberOfBasisFn=6;
	myTestHistogram.appendSlot(new taylorSlot(bounds3, myTestNumberOfBasisFn));
	
	double stepWidth=0.001; double variable=lowerBound+stepWidth/2.;
	while(variable<upperBound)
	{
		myTestHistogram.sample(variable,cos(variable));
		variable+=stepWidth;
	}
	
	basisSlot* combo=myTestHistogram.combinedSlot(1, 2);
	
	sput_fail_unless(isAround((myTestHistogram.getSlot(0) -> sampledIntegral()),-0.24299885082878656), "sampled integral has the right value");
	sput_fail_unless(isAround(combo -> sampledIntegral(),(myTestHistogram.getSlot(0) -> sampledIntegral())), "combinedSlot integral is the same as when originally sampled");
	sput_fail_unless(isAround((myTestHistogram.getSlot(0) -> sampledIntegralVariance()),0.36561204503695705), "sampled variance has the right value");
	sput_fail_unless(isAround(combo -> sampledIntegralVariance(),(myTestHistogram.getSlot(0) -> sampledIntegralVariance())), "combinedSlot variance is the same as when originally sampled");
	sput_fail_unless(isAround((myTestHistogram.getSlot(0) -> sampledIntegralError()),0.015612218399637232), "sampled error has the right value");
	sput_fail_unless(isAround(combo -> sampledIntegralError(),(myTestHistogram.getSlot(0) -> sampledIntegralError())), "combinedSlot error is the same as when originally sampled");
	
	double slotWidth=0.1; int numberOverlaps=1; int totalNumOfBasisFn=0;
	vector<basisSlot*> histogramVector=generateBasisSlots(lowerBound, upperBound, slotWidth, numberOverlaps, totalNumOfBasisFn);
	histogramBasis myTestHistogram2(histogramVector);
	slotBounds bounds4(1.4, 2.0);
	myTestHistogram2.appendSlot(new taylorSlot(bounds4, myTestNumberOfBasisFn));
	
	variable=lowerBound+stepWidth/2.;
	while(variable<upperBound)
	{
		myTestHistogram2.sample(variable,1 - 3*variable/2. + 2*variable*variable - variable*variable*variable/2.);
		variable+=stepWidth;
	}
	
	combo=myTestHistogram2.combinedSlot(4, 9);
	sput_fail_unless(isAround(combo -> sampledIntegral(),(myTestHistogram2.getSlot(15) -> sampledIntegral())), "big combinedSlot integral is the same as when originally sampled");
	sput_fail_unless(isAround(combo -> sampledIntegralVariance(),(myTestHistogram2.getSlot(15) -> sampledIntegralVariance())), "big combinedSlot variance is the same as when originally sampled");
	sput_fail_unless(isAround(combo -> sampledIntegralError(),(myTestHistogram2.getSlot(15) -> sampledIntegralError())), "big combinedSlot error is the same as when originally sampled");
}


int main(int argc, char **argv) {
	
	sput_start_testing();
	
	sput_enter_suite("test_histogram()");
	sput_run_test(test_histogram);
	sput_enter_suite("test_combinedSlot()");
	sput_run_test(test_combinedSlot);
	
	sput_finish_testing();
	
	return sput_get_return_value();
	
}