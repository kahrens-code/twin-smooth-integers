#ifndef SMOOTHNESS_H
#define SMOOTHNESS_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//pre-computation for smoothness sieving
void preSmoothness (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes, unsigned long relevantSteps, unsigned short **maxExponents, unsigned short numberRoots, short *roots, unsigned short *maxRoot, unsigned long **logTable, unsigned short **logSmoothPrimes, unsigned long start, unsigned int size, unsigned long **smoothInterval, unsigned short **smoothNumbers, unsigned short maxNumberResults, unsigned long **smoothIntsModC, unsigned short truncation, unsigned short *tolerance, unsigned int *surplusSmooth);

//included in preSmoothness, only needed for testing
void findTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned long *logTable, unsigned short *logSmoothPrimes, unsigned long start, unsigned int size, unsigned long *smoothInterval, unsigned short *smoothNumbers, unsigned long *smoothIntsModC, unsigned short truncation, unsigned short *tolerance, unsigned int *surplusSmooth);

//regular sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned long start, unsigned int size, unsigned long *smoothInterval, unsigned short *smoothNumbers);

//truncated log sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned long *logTable, unsigned long start, unsigned int size, unsigned short *smoothNumbers, unsigned short truncation, unsigned short tolerance);

//truncated log sieve2, still in test phase
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findTruncLogSmoothNumbers2 (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned long *logTable, unsigned long start, unsigned int size, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short tolerance);

#endif
