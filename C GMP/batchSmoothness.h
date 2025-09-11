#ifndef BATCHSMOOTHNESS_H
#define BATCHSMOOTHNESS_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>
	
//product of array as binary tree
void treeProduct (unsigned int arrayLength, mpz_t *array, mpz_t *product);
void treeProduct_ui (unsigned int arrayLength, unsigned int *array, mpz_t *product);

//naive product of array
void naiveProduct (unsigned int arrayLength, mpz_t *array, mpz_t *product);
void naiveProduct_ui (unsigned int arrayLength, unsigned int *array, mpz_t *product);

//product of array as binary tree
//save intermediate steps
//only for arrays of length a power of two
void treeProductSave (unsigned int arrayLength, mpz_t *array, mpz_t *tree);

//input modulo all elements in array
void naiveMod (mpz_t input, unsigned int arrayLength, mpz_t *array, mpz_t *modArray);

//input modulo all leaves in tree as binary tree 
//only for arrays of length a power of two
void treeMod (mpz_t input, unsigned int arrayLength, mpz_t *tree, mpz_t *modArray);

//find all smooth elements in array
//test elements in array for smoothness given a list of all smooth primes 
void smoothBatchFrKlMoWiList (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned int arrayLength, mpz_t *array, unsigned int *numberSmoothElements, mpz_t **smoothnessArray);
void smoothBatchBernsteinList (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned int arrayLength, mpz_t *array, unsigned int *numberSmoothElements, mpz_t **smoothnessArray);

//test elements in array for smoothness given the product of all smooth primes
void smoothBatchFrKlMoWiProduct (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothElements, mpz_t *smoothnessArray);
void smoothBatchBernsteinProduct (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothElements, mpz_t *smoothnessArray);

//test segments of array for smoothness given the product of all smooth primes 
void smoothBatchFrKlMoWiProductSegments (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, unsigned int numberSegments, unsigned int segmentSize, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothSegments, unsigned int *smoothSegments);
void smoothBatchBernsteinProductSegments (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t largestElement, unsigned int numberSegments, unsigned int segmentSize, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothSegments, unsigned int *smoothSegments);

#endif
