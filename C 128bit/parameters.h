#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<stdio.h>

typedef struct {					//define struct PTEsolution
	unsigned short degree;			//degree of a(x) and b(x)
	unsigned short numberRoots;		//can be less than 2*degree if there are double roots
	short *roots;					//all roots of a(x) and b(x) without multiplicity
	short *polyA;					//roots of a(x) with multiplicity
	short *polyB;					//roots of b(x) with multiplicity
	unsigned long C;					//difference C=|a(x)-b(x)|
	unsigned short numberPrimes;	//number of example twin smooth primes
	unsigned __int128 *startPrimes; 	//start values to find those primes within 2^20 steps
}PTEsolution;

//declare available examples
extern PTEsolution n6ex1, n6ex2, n6ex3, n6ex4, n6ex5;	//examples with degree n=6
extern PTEsolution n7ex1;								//example with degree n=7
extern PTEsolution n8ex1;								//example with degree n=8

//set all parameters accordingly to the chosen example PTE solution
//option for startParameter:
//0: 								start stays unchanged
//1, ..., PTEsolution.numerPrimes:	set start to corresponding value
void setParameters (PTEsolution *PTE, unsigned short *degree, unsigned short *numberRoots, short **roots, short **polyA, short **polyB, unsigned long *C, unsigned __int128 *start, unsigned short startParameter, unsigned __int128 totalSize, unsigned __int128 *end);

#endif
