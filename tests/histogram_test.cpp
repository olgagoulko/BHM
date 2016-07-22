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
	
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(-1.)  == 4.5, "4.5 in excess bin");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(6.)  == 4.5, "4.5 in excess bin");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(10.)  == 4.5, "4.5 in excess bin");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0), 64.0), "64.0 at zero");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(2.5) == 0, "0 at 2.5 where not sampled");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(4.88) == 0, "0 at 4.88 where not sampled");
	sput_fail_unless(myTestHistogram.sampledFunctionValueAverage(4.999999999999) == 0, "0 just below 5.0 where not sampled");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.55),8321./625.), "value at 0.55");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.52),1028866./78125.), "value at 0.52");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.2),-2816./1875.), "value at 0.2");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(0.23),-636292./234375.), "value at 0.23");
	
	//due to numerical precision (last digit) testing at slot boundaries is impossible --> this will never matter in real life
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.900000000000001),7168./625.), "value at 3.9");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.95),12352./625.), "value at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4),6816./125.), "value at 4");
	
	myTestHistogram.sample(4,-0.4);
	myTestHistogram.sample(4.23,1.1);
	
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.65),-4209./15625.), "value at 3.65");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(3.95),69272481./3906250.), "value at 3.95");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4),190765789./3906250.), "value at 4");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.16),187270944./244140625.), "value at 4.16");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.23),474382003./195312500.), "value at 4.23");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.300000000000001),881504./390625.), "value at 4.3");
	sput_fail_unless(isAround(myTestHistogram.sampledFunctionValueAverage(4.44),459764589./488281250.), "value at 4.44");

}


int main(int argc, char **argv) {
	
	sput_start_testing();
	
	sput_enter_suite("test_histogram()");
	sput_run_test(test_histogram);
	
	sput_finish_testing();
	
	return sput_get_return_value();
	
}