#ifndef SMOOTHNESS_H
#define SMOOTHNESS_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>

//find primes below smoothness bound (for testing)
void findSmoothPrimes (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes);

//general pre-computation for smoothness
void preSmoothness (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes, mpz_t **logTable, unsigned short maxNumberResults, mpz_t **smoothIntsModC, unsigned short numberRoots, short *roots, unsigned short *maxRoot);

//pre-computation for smoothness sieving
void preSievingCompare (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t maxSizePrimePower, unsigned short **maxExponents, mpz_t *logTable, unsigned short **logSmoothPrimes, mpz_t start, unsigned int size, mpz_t **smoothInterval, unsigned short **smoothNumbers, unsigned short truncPrime, unsigned short truncPower, unsigned short **minExponents, unsigned short *tolerance0, unsigned short *tolerance, unsigned short *tolerancePower, unsigned int *surplusSmooth0, unsigned int *surplusSmooth, unsigned int *surplusSmoothPower);

//pre-computation for smoothness prime truncated log sieving
void preSievingPrimeTruncLog (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t maxSizePrimePower, unsigned short **maxExponents, mpz_t *logTable, unsigned short **logSmoothPrimes, mpz_t start, unsigned int size, mpz_t **smoothInterval, unsigned short **smoothNumbers, unsigned short truncPrime, unsigned short *tolerance, unsigned int *surplusSmooth);

//included in preSmoothness, only needed for comparing
void findTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t *logTable, unsigned short *logSmoothPrimes, mpz_t start, unsigned int size, mpz_t *smoothInterval, unsigned short *smoothNumbers, unsigned short truncation, unsigned short *tolerance, unsigned int *surplusSmooth);

//only needed for comparing with power truncated log sieve
void findPowerTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t *logTable, unsigned short *logSmoothPrimes, mpz_t start, unsigned int size, mpz_t *smoothInterval, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short *tolerance, unsigned int *surplusSmooth);

//regular sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t start, unsigned int size, mpz_t *smoothInterval, unsigned short *smoothNumbers);

//prime truncated log sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findPrimeTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, mpz_t *logTable, mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short truncation, unsigned short tolerance);

//power truncated log sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findPowerTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, mpz_t *logTable, mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short tolerance);

#endif
