#ifndef NEW_H
#define NEW_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>

#include"batchSmoothness.h"

//Chinese remainder theorem
//gives smallest non-negative representative
void CRT (
	unsigned int number,	//number of different moduli
	mpz_t *values, 	//array of residues of each modulus
	mpz_t *n,	//array of moduli
	mpz_t NN,	//product of moduli
	mpz_t *result	//pointer to set to output of CRT
	);
void CRT_ui (
	unsigned int number,	//number of different moduli
	unsigned long *values, 	//array of residues of each modulus
	mpz_t *n,	//array of moduli
	mpz_t NN,	//product of moduli
	mpz_t *result	//pointer to set to output of CRT
	);

//factor input by trial division into powers of smooth primes
void factorIntoSmoothPrimes (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t *logTable, mpz_t input, mpz_t **factors, unsigned int *numberFactors, mpz_t *factorSum);

//find residue classes such that poly(res) = 0 mod prime power factors of C	
void findPrimeResidues (mpz_t *factors, unsigned int numberFactors, unsigned long maxNumberPrimeResidues, unsigned short degree, short *poly, mpz_t *numberResidues, unsigned long **primeResidues, unsigned long **numbersPrimeResidues, unsigned long **indices);

//create ordered list of residue classes such that poly(res) = 0 mod C	 	
void createResidueList (mpz_t *factors, unsigned int numberFactors, unsigned long numberResidues_ui, unsigned long *primeResidues, unsigned long *numbersPrimeResidues, unsigned long *indices, mpz_t C, mpz_t **residues);

//find a residue class such that poly(res) = 0 mod C in case there are too many to save them as a list
int noListResidue (mpz_t *factors, unsigned int numberFactors, unsigned short degree, short *poly, mpz_t C, mpz_t *currentPositions, mpz_t *residue);


//fixed interval:

//for C < size
//check if the residue classes have smooth representatives in a fixed interval [start, start + size]
void checkAllResiduesFixedInterval_SmallC (unsigned short *smoothNumbers, unsigned long numberResidues_ui, mpz_t *residues, mpz_t start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults);

//for C >= size
//find relevant steps to cover the interval
unsigned int findRelevantSteps (unsigned long numberResidues_ui, mpz_t *residues, unsigned int size, mpz_t C);

//for C >= size
//check if the residue classes have smooth representatives in a fixed interval [start, start + size]
void checkAllResiduesFixedInterval (unsigned short *smoothNumbers, unsigned long numberResidues_ui, mpz_t *residues, mpz_t start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, mpz_t C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults);


//fixed residue (for numberResidues > unsigned long):

//find largest power of two such that the product of all elements in an array of that length is still smaller than then product of the smooth primes
void findBatchSize (mpz_t start, mpz_t productPrimes, unsigned int *batchLength);

//check smoothness of representatives of a fixed residue class
void checkFixedResidueFrKlMoWi (mpz_t productPrimes, mpz_t residue, mpz_t start, unsigned int batchLength, mpz_t *batch, mpz_t *productTreeArray, mpz_t *modArray, mpz_t *smoothnessArray, unsigned short numberRoots, mpz_t *roots_gmp, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults);
void checkFixedResidueBernstein (mpz_t productPrimes, mpz_t residue, mpz_t start, unsigned int batchLength, mpz_t *batch, mpz_t *productTreeArray, mpz_t *modArray, mpz_t *smoothnessArray, unsigned short numberRoots, mpz_t *roots_gmp, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults);

#endif
