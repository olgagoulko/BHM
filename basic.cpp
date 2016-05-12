#include "basic.hpp"

using namespace std;

/*---------------------------------------------  Implementing passing options  -------------------------------------------------*/

double get_option(int inputN,char *inputV[],char *was) {
	int n;
	char option[20];
	sprintf(option,"-%s",was);
	for (n=1;n<(inputN-1);n++)
		{
		if (strcmp(inputV[n],option)==0)
			return (double) atof(inputV[n+1]);
		}
	return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------*/

bool acceptreject(double probability, gsl_rng* RNG){

bool accept=false;
double random = gsl_rng_uniform(RNG);
if(random<abs(probability)) accept=true;

return accept;

}
	
/*-------------------------------------------------------------------------------------------------------------------------------*/

bool isAround(double whatWeHave, double whatItShouldBe)
{
bool result=false;
if(abs(whatItShouldBe-whatWeHave)/whatItShouldBe<=VERY_SMALL_NUMBER) result=true;
return result;
}

	
/*---------------------------------------------------------  Sign function  -----------------------------------------------------*/

int whatsign(double a) {

int s=1;
if(a<0) {s=-1;}
return s;

}

/*--------------------------------------------------------  Rounding  -----------------------------------------------------------*/

int rounding(double a) {
return int(a + 0.5);
}

/*------------------------------------------------------  RNG seed  -------------------------------------------------------------*/

long unsigned int rdtsc(){
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((long unsigned int)hi << 32) | lo;
}

/*-------------------------------------------------------------------------------------------------------------------------------*/


