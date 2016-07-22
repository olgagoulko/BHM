#include "basic.hpp"
#include "slot.hpp"
#include "sput.h"

using namespace std; 

static void test_basisSlot()
	{
	double lowerBound=1.; double upperBound=2.;
	slotBounds myTestBounds(lowerBound, upperBound);
	int myTestNumberOfBasisFn=3;
	basisSlot myTestSlot(myTestBounds, myTestNumberOfBasisFn);
	slotBounds myTestOpenBounds(lowerBound);
	basisSlot myTestOpenSlot(myTestOpenBounds, myTestNumberOfBasisFn);
		
	sput_fail_unless(myTestSlot.checkIfInBasisSlot(1.0)  == true, "1.0 is in [1:2) slot");
	sput_fail_unless(myTestSlot.checkIfInBasisSlot(1.5)  == true, "1.5 is in [1:2) slot");
	sput_fail_unless(myTestSlot.checkIfInBasisSlot(2.0)  == false, "2.0 is not in [1:2) slot");
	sput_fail_unless(myTestSlot.checkIfInBasisSlot(2.5)  == false, "2.5 is not in [1:2) slot");
	sput_fail_unless(myTestSlot.checkIfInBasisSlot(0.5)  == false, "0.5 is not in [1:2) slot");
	sput_fail_unless(myTestOpenSlot.checkIfInBasisSlot(0.5)  == false, "0.5 is not in [1:inf) slot");
	sput_fail_unless(myTestOpenSlot.checkIfInBasisSlot(10.5)  == true, "10.5 is in [1:inf) slot");
	}
	
static void test_taylorSlot()
	{
	double lowerBound=1.; double upperBound=2.;
	slotBounds myTestBounds(lowerBound, upperBound);
	int myTestNumberOfBasisFn=4;
	taylorSlot myTestSlot(myTestBounds, myTestNumberOfBasisFn);
		
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(0, 1.5),1), "GramSchmidtFn 1 of 1.5");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(1, 1.5),0), "GramSchmidtFn 2 of 1.5");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(2, 1.5),-sqrt(5)/2.), "GramSchmidtFn 3 of 1.5");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(3, 1.5),0), "GramSchmidtFn 4 of 1.5");
	
	sput_fail_unless(isAround(myTestSlot.pairwiseIntegral(0, 0),1), "pairwiseIntegral(0,0)");
	sput_fail_unless(isAround(myTestSlot.pairwiseIntegral(0, 1),3./2.), "pairwiseIntegral(0,1)");
	sput_fail_unless(isAround(myTestSlot.pairwiseIntegral(1, 3),31./5.), "pairwiseIntegral(1,3)");
	
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(0,1.),1.), "GramSchmidtFn 1 value at 1 is 1");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(1,1.),-sqrt(3)), "GramSchmidtFn 2 value at 1 is -sqrt(3)");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(1,1.5),0), "GramSchmidtFn 2 value at 1.5 is 0");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(2,1.),sqrt(5)), "GramSchmidtFn 3 value at 1 is sqrt(5)");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(2,1.5),-sqrt(5)/2.), "GramSchmidtFn 3 value at 1.5 is -sqrt(5)/2");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(3,1.),-sqrt(7)), "GramSchmidtFn 4 value at 1 is -sqrt(7)");
	sput_fail_unless(isAround(myTestSlot.GramSchmidtBasisFn(3,1.5),0), "GramSchmidtFn 4 value at 1.5 is 0");
	
	myTestSlot.sample(1.,5.);
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.),80.), "sampled function value at 1. (after sampling 5. at 1.) is 80.");
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.5),-7.5), "sampled function value at 1.5 (after sampling 5. at 1.) is -7.5");
	myTestSlot.sample(1.5,2.);
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.),77.), "sampled function value at 1. (after sampling additionally 2. at 1.5) is 77.");
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.5),-3.), "sampled function value at 1.5 (after sampling additionally 2. at 1.5) is -3");
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.2),21./5.), "sampled function value at 1.2 (after sampling additionally 2. at 1.5) is 21/5");
	sput_fail_unless(isAround(myTestSlot.sampledFunctionValue(1.8),57./5.), "sampled function value at 1.8 (after sampling additionally 2. at 1.5) is 57/5");
	
	lowerBound=5.; upperBound=6.;
	slotBounds myTestBounds2(lowerBound, upperBound);
	myTestNumberOfBasisFn=9;
	taylorSlot myTestSlot2(myTestBounds2, myTestNumberOfBasisFn);
	sput_fail_unless(isAround(myTestSlot2.GramSchmidtBasisFn(3,5.5),0.), "GramSchmidtFn 4 value at 5.5 is 0");
	sput_fail_unless(isAround(myTestSlot2.GramSchmidtBasisFn(3,5.8),-9.*sqrt(7)/25.), "GramSchmidtFn 4 value at 5.8 is -9 sqrt(7)/25");
	sput_fail_unless(isAround(myTestSlot2.GramSchmidtBasisFn(6,5.2),2689.*sqrt(13)/15625.), "GramSchmidtFn 7 value at 5.2 is 2689 sqrt(13)/15625");
	sput_fail_unless(isAround(myTestSlot2.GramSchmidtBasisFn(8,5.1),-166553.*sqrt(17)/10000000.), "GramSchmidtFn 9 value at 5.1 is -166553*sqrt(17)/10000000");
	sput_fail_unless(isAround(myTestSlot2.GramSchmidtBasisFn(8,5),sqrt(17)), "GramSchmidtFn 9 value at 5.0 is sqrt(17)");
	}

	
	
int main(int argc, char **argv) {
		
	sput_start_testing();
		
	sput_enter_suite("test_basisSlot()");
	sput_run_test(test_basisSlot);
	sput_enter_suite("test_taylorSlot()");
	sput_run_test(test_taylorSlot);

	sput_finish_testing();
		
	return sput_get_return_value();
		
	}