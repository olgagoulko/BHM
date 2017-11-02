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
if(probability<0) {cerr << "ERROR: negative probability in acceptreject" << endl; throw std::runtime_error("ERROR: negative probability in acceptreject"); }

if(random<probability) accept=true;

return accept;

}
	
/*-------------------------------------------------------------------------------------------------------------------------------*/

bool isAround(double whatWeHave, double whatItShouldBe, double accuracy)
{
bool result=false;
if((whatItShouldBe==0)&&(abs(whatWeHave)<accuracy)) result=true;
else if(abs((whatItShouldBe-whatWeHave)/whatItShouldBe)<=accuracy) result=true;
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
int rounded=0;
if(a>=0) rounded=int(a + 0.5);
else rounded=-int(-a + 0.5);
return rounded;
}

/*------------------------------------------------------  RNG seed  -------------------------------------------------------------*/

long unsigned int rdtsc(){
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((long unsigned int)hi << 32) | lo;
}

/*-------------------------------------------------------------------------------------------------------------------------------*/


