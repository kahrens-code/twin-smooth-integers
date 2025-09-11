#include<stdio.h>
#include<math.h>
#include<gmp.h>
#include<time.h>
#include <stdlib.h>

#include"parameters.h"
#include"smoothness.h"
#include"batchSmoothness.h"
#include"naive.h"
#include"new.h"


int main() {

	printf("Test different ways to check smoothness\n");

//Initialize and set parameters and variables
	printf("\nInitalize and set\n");

	//general parameters
	mpz_t start; mpz_init(start); mpz_ui_pow_ui(start, 2, 78);//78);//60);//42);	//start of search interval I
	mpz_t temp; mpz_init(temp);
	mpz_ui_pow_ui(temp, 2, 27);//27);//25);//17);
	unsigned int smoothnessBound = mpz_get_ui(temp);	//smoothness bound
	unsigned short maxNumberResults = 1000;	//fix the size for the array of results
	unsigned int testSize = 1000000;	//number of relevant candidates that are tested for smoothness
	unsigned short togglePrimeResidues = 1;	//toggle computation of residue classes mod prime power factors (1 = on, 0 = off), in case it is too large for RAM
	unsigned short toggleResidues = 0;	//toggle computation of residue classes mod C (1 = on, 0 = off), in case it is too large for RAM
	
	//PTE solution parameters
	unsigned short degree;		//degree of a(x) and b(x)
	unsigned short numberRoots;	//can be less than 2*degree if there are double roots
	short *roots;			//all roots of a(x) and b(x) without multiplicity
	short *polyA;			//roots of a(x) with multiplicity
	short *polyB;			//roots of b(x) with multiplicity
	mpz_t C; mpz_init(C);		//difference C=|a(x)-b(x)|
	setParameters(&n7ex5, &degree, &numberRoots, &roots, &polyA, &polyB, &C);	//set PTE specific parameters

	//parameters for fixed interval
	mpz_t maxSizePrimePower; mpz_init_set(maxSizePrimePower, temp);	//define kind of smoothness for smoothness sieve:
		//powersmooth: maxSizePrimePower = smoothnessBound (restrictive)
		//smooth: maxSizePrimePower = end (inefficient)
	mpz_ui_pow_ui(temp, 2, 20);
	unsigned int size = mpz_get_ui(temp);	//size of chunks
	unsigned short truncPrime = 1;	//number of small primes to be omitted in the primeTruncLogSieve
	unsigned short truncPower = 10;	//power of 2 up to which prime powers are omitted in the powerTruncLogSieve
	
	//parameters for fixed residue
	unsigned short chooseFixedResidues = 3;	//set mode for choosing fixedResidue
		//0: set fixedResidue = 0
		//1: set fixedResidue non-trivial
		//2: set fixedResidue randomly
		//>2: run multiple tests with different residues
	
	//general variables
	mpz_t *logTable;	//array with thresholds for rounded logarithms
	mpz_t currentStart;	mpz_init(currentStart);	//start of current chunk/batch
	unsigned int numberSmoothPrimes;	//number of primes below the smoothness bound
	unsigned int *smoothPrimes;	//pointer to array of such primes
	unsigned int numberFactors;	//number of prime power factors of C
	mpz_t *factors;	//pointer to array of such factors
	mpz_t factorSum; mpz_init(factorSum);	//sum of these factors (if small enough -> maxNumberPrimeResidues)
	unsigned long maxNumberPrimeResidues;	//maximal number of residue classes r mod prime power factors of C s.t. a(r)/C and b(r)/C are integers based on the factorization of C
	unsigned long *numbersPrimeResidues;	//pointer to array of actual numbers of residue classes r per prime power factors of C s.t. a(r)/C and b(r)/C are integers 
	unsigned long *primeResidues;	//pointer to array of such classes
	unsigned long *indices;	//structure of primeResidues (indices[i] is the first residue corresponding to factors[i])
	mpz_t numberResidues; mpz_init(numberResidues);	//number |R| of residue classes r mod C s.t. a(r)/C and b(r)/C are integers
	unsigned long numberResidues_ui;
	mpz_t *residues;	//pointer to array of such classes
	unsigned short numberSmoothIntsModC;	//number of l in I s.t. a(l)/C and b(l)/C are smooth integers
	mpz_t *smoothIntsModC;	//pointer to array of such elements
	unsigned short maxRoot;	//largest root of the polynomials a(x) and b(x)
	unsigned long i;
	unsigned short j;

	//varariables for fixed interval
	mpz_t totalSize; mpz_init(totalSize);	//size of the whole search interval
	mpz_t end; mpz_init(end);	//end of search interval I
	unsigned int relevantSteps;	//number of "consecutive classes" needed to cover a chunk
	unsigned int numberChunks;	//number of chunks for full test
	mpz_t *smoothInterval;	//auxiliary array for regular sieve
	unsigned short *smoothNumbers;	//array of 1s and 0s for smooth and non-smooth integers, respectively
	unsigned short *logSmoothPrimes;	//rounded logarithms of smooth primes
	unsigned short *maxExponents;	//highest power of each prime to be included in the search
	unsigned short *minExponents;	//smallest power of each prime above 2^truncPower
	unsigned short tolerance;	//size of allowed non-smooth factor (2^tolerance) for prime truncated log sieve
	unsigned short tolerance0;	//size of allowed non-smooth factor (2^tolerance) for not truncated log sieve
	unsigned short tolerancePower;	//size of allowed non-smooth factor (2^tolerance) for power truncated log sieve
	unsigned int surplusSmooth;	//approximate number of additional smooth integers per chunk when using primeTruncLogSieve
	unsigned int surplusSmooth0;	//approximate number of additional smooth integers per chunk when using logSieve
	unsigned int surplusSmoothPower;	//approximate number of additional smooth integers per chunk when using powerTruncLogSieve

	//varariables for fixed residue
	unsigned int numberBatches;	//number of batches for full test
	unsigned int batchLength;	//size of batches for batchSmoothness
	mpz_t *batch;	//auxiliary array to reduce number of allocations
	mpz_t *productTreeArray;	//auxiliary array for tree product
	mpz_t *modArray;	//auxiliary array for tree mod
	mpz_t *smoothnessArray;	//auxiliary array for smoothness test
	mpz_t *roots_gmp;	//roots as mpz_t for easier subtraction
	mpz_t productPrimes; mpz_init(productPrimes);	//product of all smooth primes
	unsigned int log2productPrimes;	//log2(productPrimes)
	mpz_t fixedResidue; mpz_init(fixedResidue);	//fixed residue used for this approach
	unsigned long *chosenPrimeResidues; //chosen prime residues for finding fixedResidue using CRT (in case there are too many residues to write them as a list)
	mpz_t *chosenPositions;	//chosen start points for search for next residue using noListResidue
	mpz_t *multipleResidues;	//array with multiple values for fixedResidue (for chooseFixedResidues > 2)
	
	//Timings
	clock_t initialPre, finalPre;
	clock_t initialPreInterval, finalPreInterval;
	clock_t initialPreResidue, finalPreResidue;
	clock_t initialSieve, finalSieve;
	clock_t initialLogSieve, finalLogSieve;
	clock_t initialPrimeTruncLogSieve, finalPrimeTruncLogSieve;
	clock_t initialPowerTruncLogSieve, finalPowerTruncLogSieve;
	clock_t initialNaive, finalNaive;
	clock_t initialNew, finalNew;
	clock_t initialFrKlMoWi, finalFrKlMoWi;
	clock_t initialBernstein, finalBernstein;
	clock_t initialNaiveLong, finalNaiveLong;
	clock_t initialNewLong, finalNewLong;
	clock_t initialFrKlMoWiLong, finalFrKlMoWiLong;
	clock_t initialBernsteinLong, finalBernsteinLong;
	
	//print parameters
	printf("Parameters:\n");
	gmp_printf("smoothnessBound: 2^%.2f, start: %Zd (size in base 2: %lu), test size: %u\n", log2(smoothnessBound), start, mpz_sizeinbase(start, 2), testSize);
	gmp_printf("degree: %u, C: %Zd\n", degree, C);
	printf("a(x) = ");
	for (i = 0; i < degree; i++) {
		printf("(x - %d) ", polyA[i]);
	}printf("\n");
	printf("b(x) = ");
	for (i = 0; i < degree; i++) {
		printf("(x - %d) ", polyB[i]);
	}printf("\n");
	
	
//Pre-computations for smoothness tests
	printf("\nSmoothness pre-computation\n");
	initialPre = time(NULL);
	maxNumberPrimeResidues = 0;	//use maxNumberPrimeResidues as flag
	numberResidues_ui = 0;	//use numberResidues_ui as flag
	preSmoothness(smoothnessBound, &numberSmoothPrimes, &smoothPrimes, &logTable, maxNumberResults, &smoothIntsModC, numberRoots, roots, &maxRoot);
	factorIntoSmoothPrimes(numberSmoothPrimes, smoothPrimes, logTable, C, &factors, &numberFactors, &factorSum);
	printf("number of smooth primes: %u", numberSmoothPrimes); fflush(stdout);
	if (mpz_fits_ulong_p(factorSum) && togglePrimeResidues == 1) {
		//printf("find prime residues ... "); fflush(stdout);
		maxNumberPrimeResidues = mpz_get_ui(factorSum);
		findPrimeResidues(factors, numberFactors, maxNumberPrimeResidues, degree, polyA, &numberResidues, &primeResidues, &numbersPrimeResidues, &indices);
		//printf("done\n");
		gmp_printf(", number of residues: %Zd\n", numberResidues);
		if (mpz_fits_ulong_p(numberResidues) && toggleResidues == 1) {
			printf("find residues ... "); fflush(stdout);
			numberResidues_ui = mpz_get_ui(numberResidues);
			createResidueList(factors, numberFactors, numberResidues_ui, primeResidues, numbersPrimeResidues, indices, C, &residues);
			printf("done\n");
		}else {
			if (toggleResidues != 1) {
				printf("Warning: computation of residues was toggled off. Can not use new fixed interval.\n");
			} else {
				printf("Warning: number of residues mod C does not fit in unsigned long. Can not allocate list of residues. Can not use new fixed interval.\n");
			}
		}
	}else {
		if (togglePrimeResidues != 1) {
			printf("\nWarning: computation of primeResidues was toggled off. Can not use new fixed interval.\n");
		} else {
			printf("\nWarning: number of possible residues mod prime power factors of C does not fit in unsigned long. Can not allocate list of primeResidues. Can not use new fixed interval.\n");
		}
	}
	finalPre = time(NULL);
	printf("general pre-computation took %lu s\n", finalPre - initialPre);

	//for fixed interval
	printf("for fixed interval:\n");
	initialPreInterval = clock();
	preSievingCompare(numberSmoothPrimes, smoothPrimes, maxSizePrimePower, &maxExponents, logTable, &logSmoothPrimes, start, size, &smoothInterval, &smoothNumbers, truncPrime, truncPower, &minExponents, &tolerance0, &tolerance, &tolerancePower, &surplusSmooth0, &surplusSmooth, &surplusSmoothPower);
	printf("logSieve:\n");
	printf("tolerance to include all smooth integers: plain: %u, primeTrunc: %u, powerTrunc: %u\n", tolerance0, tolerance, tolerancePower);
	printf("roughly more smooth integers per chunk: plain: %u, primeTrunc: %u, powerTrunc: %u\n", surplusSmooth0, surplusSmooth, surplusSmoothPower);
	if(maxNumberPrimeResidues == 0) {
		printf("Warning: can not compute relevantSteps and set it to 1\n");
		relevantSteps = 1;
	} else {
		if (numberResidues_ui > 0) {
			relevantSteps = findRelevantSteps(numberResidues_ui, residues, size, C);
		} else {
			mpz_mul_ui(temp, numberResidues, size);
			mpz_cdiv_q(temp, temp, C);
			relevantSteps = mpz_get_ui(temp);	//relevantSteps < size since relevantSteps / C < 1
		}
	}
	printf("relevantSteps: %u\n", relevantSteps);
	mpz_set_ui(totalSize, size);
	numberChunks = (testSize + relevantSteps - 1) / relevantSteps; //ceil(testSize / relevantSteps)
	mpz_mul_ui(totalSize, totalSize, numberChunks);
	mpz_add(end, start, totalSize);
	gmp_printf("size: 2^%.2f, end: %Zd (size in base 2: %lu), number of chunks: %u\n", log2(size), end, mpz_sizeinbase(end, 2), numberChunks);
	finalPreInterval = clock();
	printf("pre-computation for fixed interval took %.3f s\n", (double)(finalPreInterval - initialPreInterval) / CLOCKS_PER_SEC);
	
	//for fixed residue
	printf("for fixed residue:\n");
	initialPreResidue = clock();
	treeProduct_ui(numberSmoothPrimes, smoothPrimes, &productPrimes);	//product of smooth primes (unsigned int) as tree 
	log2productPrimes = (unsigned int) mpz_sizeinbase(productPrimes, 2);
	printf("log2(productPrimes): %u\n", log2productPrimes);
	findBatchSize (start, productPrimes, &batchLength);	//find largest power of two such that the product of all elements in an array of that length is still smaller than then product of the smooth primes
	printf("batchLength: %u\n", batchLength);
	batch = (mpz_t *) malloc(batchLength * sizeof(mpz_t));
	productTreeArray = (mpz_t *) malloc((2 * batchLength - 1) * sizeof(mpz_t));
	modArray = (mpz_t *) malloc((batchLength) * sizeof(mpz_t));
	smoothnessArray = (mpz_t *) malloc((batchLength) * sizeof(mpz_t));
	roots_gmp = (mpz_t *) malloc((numberRoots) * sizeof(mpz_t));
	if (batch == NULL || productTreeArray == NULL || modArray == NULL || smoothnessArray == NULL || roots_gmp == NULL) {
		printf("Error: auxiliary arrays for fixed residue appraoch could not be allocated\n");
		exit(-1);
	}
	for (i = 0; i < numberRoots; i++) {
		mpz_init_set_si(roots_gmp[i], roots[i]);
	}
	for (i = 0; i < batchLength; i++){
		mpz_init(batch[i]);
		mpz_init(modArray[i]);
		mpz_init(smoothnessArray[i]);
	}
	for (i = 0; i < 2 * batchLength - 1; i++){
		mpz_init(productTreeArray[i]);
	}
	numberBatches = (testSize + batchLength - 1) / batchLength;	//ceil(testSize / batchLength)
	mpz_mul_ui(temp, C, numberBatches * batchLength);
	mpz_add(temp, temp, start);
	gmp_printf("number of rounds: %u, end: %Zd\n", numberBatches, temp);
	srand(time(NULL));
	if (chooseFixedResidues > 2) {
		multipleResidues = (mpz_t *) malloc((10) * sizeof(mpz_t));
		if (multipleResidues == NULL) {
			printf("Error: multipleResidues could not be allocated\n");
			exit(-1);
		}
		for (i = 0; i < 10; i++) {
			mpz_init(multipleResidues[i]);
		}
	}
	if (maxNumberPrimeResidues > 0) {
		if (numberResidues_ui > 0) {
			if(chooseFixedResidues == 0) { //choose fixedResidue = 0
				mpz_set(fixedResidue, residues[0]);
			}
			if(chooseFixedResidues == 1) { //choose non-trivial residue
				mpz_set(fixedResidue, residues[1]);
			}
			if(chooseFixedResidues == 2) { //choose random residue
				mpz_set(fixedResidue, residues[rand() % numberResidues_ui]);
			}
			if(chooseFixedResidues > 2) { //choose multiple residues
				mpz_set(multipleResidues[0], residues[0]);
				mpz_set(multipleResidues[1], residues[1]);
				mpz_set(multipleResidues[2], residues[numberResidues_ui / 2]);
				mpz_set(multipleResidues[3], residues[numberResidues_ui / 2 + 1]);
				mpz_set(multipleResidues[4], residues[numberResidues_ui - 1]);
				for (i = 5; i < 10; i++) {
					mpz_set(multipleResidues[i], residues[rand() % numberResidues_ui]);
				}
			}
		}else {
			chosenPrimeResidues = (unsigned long *) malloc(numberFactors * sizeof(long));
			if (chosenPrimeResidues == NULL) {
				printf("Error: chosenPrimeResidues could not be allocated\n");
				exit(-1);
			}
			if(chooseFixedResidues == 0) { //choose fixedResidue = 0
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i]];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &fixedResidue);
			}
			if(chooseFixedResidues == 1) { //choose non-trivial residue
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + 1];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &fixedResidue);
			}
			if(chooseFixedResidues == 2) { //choose random residue
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + (rand() % numbersPrimeResidues[i])];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &fixedResidue);
			}
			if(chooseFixedResidues > 2) { //choose multiple residues
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i]];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[0]));
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + 1];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[1]));
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + (numbersPrimeResidues[i] / 2)];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[2]));
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + (numbersPrimeResidues[i] / 2) + 1];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[3]));
				for (i = 0; i < numberFactors; i++) {
					chosenPrimeResidues[i] = primeResidues[indices[i] + (numbersPrimeResidues[i] - 1)];
				}
				CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[4]));
				for (j = 5; j < 10; j++) {
					for (i = 0; i < numberFactors; i++) {
						chosenPrimeResidues[i] = primeResidues[indices[i] + (rand() % numbersPrimeResidues[i])];
					}
					CRT_ui(numberFactors, chosenPrimeResidues, factors, C, &(multipleResidues[j]));
				}
			}
		}
	}else {
		chosenPositions = (mpz_t *) malloc(numberFactors * sizeof(mpz_t));
		if (chosenPositions == NULL) {
			printf("Error: chosenPositions could not be allocated\n");
			exit(-1);
		}
		if(chooseFixedResidues == 0) { //choose fixedResidue = 0
			for (i = 0; i < numberFactors; i++) {
				mpz_init_set_ui(chosenPositions[i], 0);
			}
			if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &fixedResidue) == 0) {
				printf("Error: could not find residue. Maybe reached last residue.\n");
				exit(-1);
			}
		}
		if(chooseFixedResidues == 1) { //choose non-trivial residue
			for (i = 0; i < numberFactors; i++) {
				mpz_init_set_ui(chosenPositions[i], 1);
			}
			if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &fixedResidue) == 0) {
				printf("Error: could not find residue. Maybe reached last residue.\n");
				exit(-1);
			}
		}
		if(chooseFixedResidues == 2) { //choose random residue
			for (i = 0; i < numberFactors; i++) {
				mpz_init_set_ui(chosenPositions[i], (rand() % 2));
			}
			if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &fixedResidue) == 0) {
				printf("Error: could not find residue. Maybe reached last residue.\n");
				exit(-1);
			}
		}
		if(chooseFixedResidues > 2) { //choose multiple residues
			for (i = 0; i < numberFactors; i++) {
				mpz_init_set_ui(chosenPositions[i], 0);
			}
			if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &(multipleResidues[0])) == 0) {
				printf("Error: could not find residue. Maybe reached last residue.\n");
				exit(-1);
			}
			for (i = 0; i < numberFactors; i++) {
				mpz_init_set_ui(chosenPositions[i], 1);
			}
			if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &(multipleResidues[1])) == 0) {
				printf("Error: could not find residue. Maybe reached last residue.\n");
				exit(-1);
			}
			for (j = 2; j < 5; j++) {
				for (i = 0; i < numberFactors; i++) {
					mpz_init_set_ui(chosenPositions[i], (rand() % 2));
				}
				if (noListResidue (factors, numberFactors, degree, polyA, C, chosenPositions, &(multipleResidues[j])) == 0) {
					printf("Error: could not find residue. Maybe reached last residue.\n");
					exit(-1);
				}
			}
		}
	}
	finalPreResidue = clock();
	printf("pre-computation for fixed residue took %.3f s\n", (double)(finalPreResidue - initialPreResidue) / CLOCKS_PER_SEC);


//Test single iterations
	printf("\nStarting single tests\n");
	
	//fixed interval
	printf("\nFixed interval approach:\n");
	
	printf("test different sieves for a single chunk:\n");
	//single chunk regular sieve for smoothness
	initialSieve = clock();
	mpz_sub_ui(currentStart, start, maxRoot);
	findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, currentStart, size, smoothInterval, smoothNumbers);
	finalSieve = clock();
	printf("single Sieve: %.3f ms\n", (double)(finalSieve - initialSieve) / CLOCKS_PER_SEC * 1000);
	//single chunk log sieve for smoothness
	initialLogSieve = clock();
	mpz_sub_ui(currentStart, start, maxRoot);
	findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, 0, tolerance0);
	finalLogSieve = clock();
	printf("single Log Sieve: %.3f ms\n", (double)(finalLogSieve - initialLogSieve) / CLOCKS_PER_SEC * 1000);
	//single chunk prime truncated log sieve for smoothness
	initialPrimeTruncLogSieve = clock();
	mpz_sub_ui(currentStart, start, maxRoot);
	findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
	finalPrimeTruncLogSieve = clock();
	printf("single Trunc Log Sieve: %.3f ms\n", (double)(finalPrimeTruncLogSieve - initialPrimeTruncLogSieve) / CLOCKS_PER_SEC * 1000);
	//single chunk power truncated log sieve for smoothness
	initialPowerTruncLogSieve = clock();
	mpz_sub_ui(currentStart, start, maxRoot);
	findPowerTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, minExponents, tolerancePower);
	finalPowerTruncLogSieve = clock();
	printf("single Power Trunc Log Sieve: %.3f ms\n", (double)(finalPowerTruncLogSieve - initialPowerTruncLogSieve) / CLOCKS_PER_SEC * 1000);
	
	printf("test ways to find twin smooth integers:\n");
	if (numberResidues_ui > 0) {
		//new approach single chunk (test only residues)
		initialNew = clock();
		numberSmoothIntsModC = 0;
		mpz_sub_ui(currentStart, start, maxRoot);
		if (mpz_cmp_ui(C, size) >= 0) {
			checkAllResiduesFixedInterval(smoothNumbers, numberResidues_ui, residues, currentStart, size, numberRoots, roots, maxRoot, C, relevantSteps, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		}else{
			checkAllResiduesFixedInterval_SmallC(smoothNumbers, numberResidues_ui, residues, currentStart, size, numberRoots, roots, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		}
		finalNew = clock();
		printf("single new: %.3f ms\n", (double)(finalNew - initialNew) / CLOCKS_PER_SEC * 1000);
	}
	//naive approach single chunk (test all smooth integers in chunk)
	initialNaive = clock();
	numberSmoothIntsModC = 0;
	mpz_sub_ui(currentStart, start, maxRoot);
	findSmoothIntsModC(currentStart, size, smoothNumbers, numberRoots, roots, maxRoot, degree, polyA, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
	finalNaive = clock();
	printf("single naive: %.3f ms\n", (double)(finalNaive - initialNaive) / CLOCKS_PER_SEC * 1000);
	
	//fixed residue
	printf("\nFixed residue approach:\n");
	//single batch using FrKlMoWi
	if (chooseFixedResidues <=2) {
		gmp_printf("chosen fixedResidue: %Zd\n", fixedResidue);
		initialFrKlMoWi = clock();
		numberSmoothIntsModC = 0;
		checkFixedResidueFrKlMoWi(productPrimes, fixedResidue, start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		finalFrKlMoWi = clock();
		printf("single FrKlMoWi: %.3f ms\n", (double)(finalFrKlMoWi - initialFrKlMoWi) / CLOCKS_PER_SEC * 1000);
	} else {
		if (maxNumberPrimeResidues > 0) {
			for (i = 0; i < 10; i++) {
				if (i == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[i]);
				initialFrKlMoWi = clock();
				numberSmoothIntsModC = 0;
				checkFixedResidueFrKlMoWi(productPrimes, multipleResidues[i], start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				finalFrKlMoWi = clock();
				printf("single FrKlMoWi: %.3f ms\n", (double)(finalFrKlMoWi - initialFrKlMoWi) / CLOCKS_PER_SEC * 1000);
			}
		} else {
			for (i = 0; i < 5; i++) {
				if (i == 2) printf("random residues:\n");			
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[i]);
				initialFrKlMoWi = clock();
				numberSmoothIntsModC = 0;
				checkFixedResidueFrKlMoWi(productPrimes, multipleResidues[i], start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				finalFrKlMoWi = clock();
				printf("single FrKlMoWi: %.3f ms\n", (double)(finalFrKlMoWi - initialFrKlMoWi) / CLOCKS_PER_SEC * 1000);
			}
		}
	}
	printf("\n");
	//single batch using Bernstein
	if (chooseFixedResidues <=2) {
		gmp_printf("chosen fixedResidue: %Zd\n", fixedResidue);
		initialBernstein = clock();
		numberSmoothIntsModC = 0;
		checkFixedResidueBernstein(productPrimes, fixedResidue, start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		finalBernstein = clock();
		printf("single Bernstein: %.3f ms\n", (double)(finalBernstein - initialBernstein) / CLOCKS_PER_SEC * 1000);
	} else {
		if (maxNumberPrimeResidues > 0) {
			for (i = 0; i < 10; i++) {
				if (i == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[i]);
				initialBernstein = clock();
				numberSmoothIntsModC = 0;
				checkFixedResidueBernstein(productPrimes, multipleResidues[i], start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				finalBernstein = clock();
				printf("single Bernstein: %.3f ms\n", (double)(finalBernstein - initialBernstein) / CLOCKS_PER_SEC * 1000);
			}
		} else {
			for (i = 0; i < 5; i++) {
				if (i == 2) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[i]);
				initialBernstein = clock();
				numberSmoothIntsModC = 0;
				checkFixedResidueBernstein(productPrimes, multipleResidues[i], start, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				finalBernstein = clock();
				printf("single Bernstein: %.3f ms\n", (double)(finalBernstein - initialBernstein) / CLOCKS_PER_SEC * 1000);
			}
		}
	}
	
//Test longer searches
	printf("\nStarting full tests\n");
	
	//fixed interval using primeTruncLogSieve
	printf("\nFixed interval approach:\n");
	if (numberResidues_ui > 0) {
		//new approach whole interval
		initialNewLong = time(NULL);
		numberSmoothIntsModC = 0;
		mpz_sub_ui(currentStart, start, maxRoot);
		if (mpz_cmp_ui(C, size) >= 0) {
			while (mpz_cmp(currentStart, end) < 0) {
				findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
				checkAllResiduesFixedInterval(smoothNumbers, numberResidues_ui, residues, currentStart, size, numberRoots, roots, maxRoot, C, relevantSteps, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				mpz_add_ui(currentStart, currentStart, size - maxRoot);
			}
		}else{
			while (mpz_cmp(currentStart, end) < 0) {
				findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
				checkAllResiduesFixedInterval_SmallC(smoothNumbers, numberResidues_ui, residues, currentStart, size, numberRoots, roots, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
				mpz_add_ui(currentStart, currentStart, size - maxRoot);
			}
		}
		finalNewLong = time(NULL);
		printf("whole interval new: %lu s\n", finalNewLong - initialNewLong);
		printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
		for (i = 0; i < numberSmoothIntsModC; i++) {
			gmp_printf("%Zd\n", smoothIntsModC[i]);	
		}
	}
	//naive appraoch whole interval
	initialNaiveLong = time(NULL);
	numberSmoothIntsModC = 0;
	mpz_sub_ui(currentStart, start, maxRoot);
	while (mpz_cmp(currentStart, end) < 0) {
		findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
		findSmoothIntsModC(currentStart, size, smoothNumbers, numberRoots, roots, maxRoot, degree, polyA, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		mpz_add_ui(currentStart, currentStart, size - maxRoot);
	}
	finalNaiveLong = time(NULL);
	printf("whole interval naive: %lu s\n", finalNaiveLong - initialNaiveLong);
	printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
	for (i = 0; i < numberSmoothIntsModC; i++) {
		gmp_printf("%Zd\n", smoothIntsModC[i]);	
	}
	
	//fixed residue
	printf("\nFixed residue approach:\n");
	//all batches using FrKlMoWi
	if (chooseFixedResidues <=2) {
		gmp_printf("chosen fixedResidue: %Zd\n", fixedResidue);
		initialFrKlMoWiLong = time(NULL);
		numberSmoothIntsModC = 0;
		mpz_set(currentStart, start);
		for (i = 0; i < numberBatches; i++) {
			checkFixedResidueFrKlMoWi(productPrimes, fixedResidue, currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
			mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
		}
		finalFrKlMoWiLong = time(NULL);
		printf("all batches FrKlMoWi: %lu s\n", finalFrKlMoWiLong - initialFrKlMoWiLong);
		printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
		for (i = 0; i < numberSmoothIntsModC; i++) {
			gmp_printf("%Zd\n", smoothIntsModC[i]);	
		}
	} else {
		if (maxNumberPrimeResidues > 0) {
			for (j = 0; j < 10; j++) {
				if (j == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[j]);
				initialFrKlMoWiLong = time(NULL);
				numberSmoothIntsModC = 0;
				mpz_set(currentStart, start);
				for (i = 0; i < numberBatches; i++) {
					checkFixedResidueFrKlMoWi(productPrimes, multipleResidues[j], currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
					mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
				}
				finalFrKlMoWiLong = time(NULL);
				printf("all batches FrKlMoWi: %lu s\n", finalFrKlMoWiLong - initialFrKlMoWiLong);
				printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
				for (i = 0; i < numberSmoothIntsModC; i++) {
					gmp_printf("%Zd\n", smoothIntsModC[i]);	
				}
			}
		} else {
			for (j = 0; j < 5; j++) {
				if (j == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[j]);
				initialFrKlMoWiLong = time(NULL);
				numberSmoothIntsModC = 0;
				mpz_set(currentStart, start);
				for (i = 0; i < numberBatches; i++) {
					checkFixedResidueFrKlMoWi(productPrimes, multipleResidues[j], currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
					mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
				}
				finalFrKlMoWiLong = time(NULL);
				printf("all batches FrKlMoWi: %lu s\n", finalFrKlMoWiLong - initialFrKlMoWiLong);
				printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
				for (i = 0; i < numberSmoothIntsModC; i++) {
					gmp_printf("%Zd\n", smoothIntsModC[i]);	
				}
			}
		}
	}
	printf("\n");
	//all batches using Bernstein
	if (chooseFixedResidues <=2) {
		gmp_printf("chosen fixedResidue: %Zd\n", fixedResidue);
		initialBernsteinLong = time(NULL);
		numberSmoothIntsModC = 0;
		mpz_set(currentStart, start);
		for (i = 0; i < numberBatches; i++) {
			checkFixedResidueBernstein(productPrimes, fixedResidue, currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
			mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
		}
		finalBernsteinLong = time(NULL);
		printf("all batches Bernstein: %lu s\n", finalBernsteinLong - initialBernsteinLong);
		printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
		for (i = 0; i < numberSmoothIntsModC; i++) {
			gmp_printf("%Zd\n", smoothIntsModC[i]);	
		}
	} else {
		if (maxNumberPrimeResidues > 0) {
			for (j = 0; j < 10; j++) {
				if (j == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[j]);
				initialBernsteinLong = time(NULL);
				numberSmoothIntsModC = 0;
				mpz_set(currentStart, start);
				for (i = 0; i < numberBatches; i++) {
					checkFixedResidueBernstein(productPrimes, multipleResidues[j], currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
					mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
				}
				finalBernsteinLong = time(NULL);
				printf("all batches Bernstein: %lu s\n", finalBernsteinLong - initialBernsteinLong);
				printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
				for (i = 0; i < numberSmoothIntsModC; i++) {
					gmp_printf("%Zd\n", smoothIntsModC[i]);	
				}
			}
		} else {
			for (j = 0; j < 5; j++) {
				if (j == 5) printf("random residues:\n");
				gmp_printf("chosen fixedResidue: %Zd\n", multipleResidues[j]);
				initialBernsteinLong = time(NULL);
				numberSmoothIntsModC = 0;
				mpz_set(currentStart, start);
				for (i = 0; i < numberBatches; i++) {
					checkFixedResidueBernstein(productPrimes, multipleResidues[j], currentStart, batchLength, batch, productTreeArray, modArray, smoothnessArray, numberRoots, roots_gmp, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
					mpz_addmul_ui(currentStart, C, batchLength);	//currentStart += C * batchLength
				}
				finalBernsteinLong = time(NULL);
				printf("all batches Bernstein: %lu s\n", finalBernsteinLong - initialBernsteinLong);
				printf("found %u elements l s.t. a(l) / C, b(l) / C are smooth integers:\n", numberSmoothIntsModC);
				for (i = 0; i < numberSmoothIntsModC; i++) {
					gmp_printf("%Zd\n", smoothIntsModC[i]);	
				}
			}
		}
	}

//clear and free
	printf("\nClear and free ... "); fflush(stdout);
	//general parameters
	mpz_clear(start);
	mpz_clear(temp);
	//PTE solution parameters
	mpz_clear(C);
	//parameters for fixed interval
	mpz_clear(maxSizePrimePower);
	//general variables
	for (i = 0; i < 128; i++) {
		mpz_clear(logTable[i]);
	}
	free(logTable);
	mpz_clear(currentStart);
	free(smoothPrimes);
	for (i = 0; i < numberFactors; i++) {
		mpz_clear(factors[i]);
	}
	free(factors);
	mpz_clear(factorSum);
	if (maxNumberPrimeResidues > 0) {
		free(numbersPrimeResidues);
		free(primeResidues);
		free(indices);
		if (numberResidues_ui > 0) {
			for (i = 0; i < numberResidues_ui; i++) {
				mpz_clear(residues[i]);
			}
			free(residues);
		}else {
			free(chosenPrimeResidues);
		}
	}else {
		for (i = 0; i < numberFactors; i++) {
			mpz_clear(chosenPositions[i]);
		}
		free(chosenPositions);
	}
	for (i = 0; i < maxNumberResults; i++) {
		mpz_clear(smoothIntsModC[i]);
	}
	free(smoothIntsModC);
	//varariables for fixed interval
	mpz_clear(totalSize);
	mpz_clear(end);
	for (i = 0; i < size; i++) {
		mpz_clear(smoothInterval[i]);
	}
	free(smoothInterval);
	free(smoothNumbers);
	free(logSmoothPrimes);
	free(maxExponents);
	free(minExponents);
	//varariables for fixed residue
	for (i = 0; i < numberRoots; i++) {
		mpz_clear(roots_gmp[i]);
	}
	free(roots_gmp);
	for (i = 0; i < batchLength; i++){
		mpz_clear(batch[i]);
		mpz_clear(modArray[i]);
		mpz_clear(smoothnessArray[i]);
	}
	free(batch);
	free(modArray);
	free(smoothnessArray);
	for (i = 0; i < 2 * batchLength - 1; i++){
		mpz_clear(productTreeArray[i]);
	}
	free(productTreeArray);
	mpz_clear(productPrimes);
	mpz_clear(fixedResidue);
	if (chooseFixedResidues > 2) {
		for (i = 0; i < 10; i++) {
			mpz_clear(multipleResidues[i]);
		}
		free(multipleResidues);
	}
	printf("done\n");
  return 0;
}
