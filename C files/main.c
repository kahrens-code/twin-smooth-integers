#include<stdio.h>
#include<math.h>
#include<time.h>

#include"parameters.h"
#include"smoothness.h"
#include"naive.h"
#include"new.h"

#define COMPARE	//compare different approaches and print running times to screen

int main() {
	
//general parameters
	unsigned short maxNumberResults = 100;	//fix the size for the arry of results
	unsigned int size = pow(2, 20);		//size of chunks
	unsigned short truncPrime = 1;		//number of small primes to be omitted in the primeTruncLogSieve
	unsigned short truncPower = 10;		//power of 2 up to which prime powers are omitted in the powerTruncLogSieve
	unsigned long start = pow(2, 42);		//start of search interval I
	unsigned int totalSize = 1;			//size of the whole search interval
	unsigned long end = start + totalSize;	//end of search interval I
	unsigned int smoothnessBound = pow(2, 16);	//smoothness bound
	unsigned long maxSizePrimePower = (unsigned long) smoothnessBound;	//define kind of smoothness:
		//powersmooth: maxSizePrimePower = smoothnessBound (restrictive)
		//smooth: maxSizePrimePower = start + size (inefficient)
		//compromise: e.g. maxSizePrimePower = size (for small smoothness bounds)
	
//PTE solution parameters
	unsigned short degree;		//degree of a(x) and b(x)
	unsigned short numberRoots;	//can be less than 2*degree if there are double roots
	short *roots;			//all roots of a(x) and b(x) without multiplicity
	short *polyA;			//roots of a(x) with multiplicity
	short *polyB;			//roots of b(x) with multiplicity
	unsigned int C;		//difference C=|a(x)-b(x)|

//set PTE specific parameters	
	unsigned short startParameter = 1;	//0: not changed, >0: set to specific start value
	setParameters(&n6ex5, &degree, &numberRoots, &roots, &polyA, &polyB, &C, &start, startParameter, totalSize, &end);

//	start = 24924000000000;
//	end  =  24925000000000;
//	start = 23000000000000;
//	end  =  23750000000000;
//	start = 23750000000000;
//	end  =  24500000000000;
//	start = 24500000000000;
//	end  =  25250000000000;
//	start = 25250000000000;
//	end  =  26000000000000;

//variables
	unsigned long currentStart;		//start of current chunk
	unsigned int numberSmoothPrimes;	//number of primes below the smoothness bound
	unsigned int *smoothPrimes;		//pointer to array of such primes
	unsigned short *maxExponents;		//highest power of each prime to be included in the search
	unsigned short *minExponents;		//smallest power of each prime above 2^truncPower
	unsigned long *logTable;		//array with thresholds for rounded logarithms
	unsigned short *logSmoothPrimes;	//rounded logarithms of smooth primes
	unsigned short tolerance;		//size of allowed non-smooth factor (2^tolerance) for prime truncated log sieve
	unsigned int surplusSmooth;		//approximate number of additional smooth integers per chunk when using primeTruncLogSieve
	#ifdef COMPARE
	unsigned short tolerance0;		//size of allowed non-smooth factor (2^tolerance) for not truncated log sieve
	unsigned short tolerancePower;	//size of allowed non-smooth factor (2^tolerance) for power truncated log sieve
	unsigned int surplusSmooth0;		//approximate number of additional smooth integers per chunk when using logSieve
	unsigned int surplusSmoothPower;	//approximate number of additional smooth integers per chunk when using powerTruncLogSieve
	#endif
	unsigned long *smoothInterval;	//auxiliary array for regular sieve
	unsigned short *smoothNumbers;	//array of 1s and 0s for smooth and non-smooth integers, respectively
	unsigned int numberResidues;		//number |R| of residue classes r mod C s.t. a(r)/C and b(r)/C are integers
	unsigned int *residues;		//pointer to array of such classes
	unsigned int relevantSteps;		//number of "consecutive classes" needed to cover a chunk
	unsigned short maxRoot;		//largest root of the polynomials a(x) and b(x)
	unsigned short numberSmoothIntsModC;	//number of l in I s.t. a(l)/C and b(l)/C are smooth integers
	unsigned long *smoothIntsModC;	//pointer to array of such elements

//timing
	clock_t initialPre, finalPre;
	clock_t initialSieve, finalSieve;
	clock_t initialLogSieve, finalLogSieve;
	clock_t initialTruncLogSieve, finalTruncLogSieve;
	clock_t initialTruncLogSieve2, finalTruncLogSieve2;
	clock_t initialNaive, finalNaive;
	clock_t initialNew, finalNew;
	time_t initialNaiveLong, finalNaiveLong;
	time_t initialNewLong, finalNewLong;

//pre-computation
	preSmoothness(smoothnessBound, &numberSmoothPrimes, &smoothPrimes, maxSizePrimePower, &maxExponents, numberRoots, roots, &maxRoot, &logTable, &logSmoothPrimes, start, size, &smoothInterval, &smoothNumbers, maxNumberResults, &smoothIntsModC, truncPrime, truncPower, &minExponents, &tolerance, &surplusSmooth);
	#ifdef COMPARE
	findTolerance (numberSmoothPrimes, smoothPrimes, maxExponents, logTable, logSmoothPrimes, start, size, smoothInterval, smoothNumbers, 0, &tolerance0, &surplusSmooth0);
	findPowerTolerance (numberSmoothPrimes, smoothPrimes, maxExponents, logTable, logSmoothPrimes, start, size, smoothInterval, smoothNumbers, minExponents, &tolerancePower, &surplusSmoothPower);
	#endif
//pre-computation for new approach
	initialPre = clock();
	findResidues(numberSmoothPrimes, smoothPrimes, degree, polyA, C, &numberResidues, &residues);
	relevantSteps = findRelevantSteps(numberResidues, residues, size, C);
	finalPre = clock();
	#ifdef COMPARE
	//single chunk regular sieve for smoothness
		initialSieve = clock();
		currentStart = start - maxRoot;
		findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, currentStart, size, smoothInterval, smoothNumbers);
		finalSieve = clock();
	//single chunk log sieve for smoothness
		initialLogSieve =clock();
		currentStart = start - maxRoot;
		findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, 0, tolerance0);
		finalLogSieve = clock();
	//single chunk prime truncated log sieve for smoothness
		initialTruncLogSieve = clock();
		currentStart = start - maxRoot;
		findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
		finalTruncLogSieve = clock();
	//single chunk truncated log sieve2 for smoothness
		initialTruncLogSieve2 = clock();
		currentStart = start - maxRoot;
		findPowerTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, minExponents, tolerance0);
		finalTruncLogSieve2 = clock();

	//naive approach single chunk
		initialNaive = clock();
		numberSmoothIntsModC = 0;
		currentStart = start - maxRoot;
		findSmoothIntsModC(currentStart, size, smoothNumbers, numberRoots, roots, maxRoot, degree, polyA, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		finalNaive = clock();
	//new approach single chunk
		initialNew = clock();
		numberSmoothIntsModC = 0;
		currentStart = start - maxRoot;
		if (C >= size) {
			checkMOD(smoothNumbers, numberResidues, residues, currentStart, size, numberRoots, roots, maxRoot, C, relevantSteps, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		}else{
			checkMod (smoothNumbers, numberResidues, residues, start, size, numberRoots, roots, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		}
		finalNew = clock();
	#endif
free(minExponents);
//naive appraoch whole interval
	initialNaiveLong = time(NULL);
	numberSmoothIntsModC = 0;
 	currentStart = start - maxRoot;
	while (currentStart < end) {
		//findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, currentStart, size, smoothInterval, smoothNumbers); 
		findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
		findSmoothIntsModC(currentStart, size, smoothNumbers, numberRoots, roots, maxRoot, degree, polyA, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
		currentStart += size - maxRoot;
	}
	finalNaiveLong = time(NULL);

//new approach whole interval
	initialNewLong = time(NULL);
	numberSmoothIntsModC = 0;
	currentStart = start - maxRoot;
	if (C >= size) {
		while (currentStart < end) {
			//findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, currentStart, size, smoothInterval, smoothNumbers);
			findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
			checkMOD(smoothNumbers, numberResidues, residues, currentStart, size, numberRoots, roots, maxRoot, C, relevantSteps, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
			currentStart += size - maxRoot;
		}
	}else{
		while (currentStart < end) {
			//findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, currentStart, size, smoothInterval, smoothNumbers);
			findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, currentStart, size, smoothNumbers, truncPrime, tolerance);
			checkMod (smoothNumbers, numberResidues, residues, start, size, numberRoots, roots, maxRoot, C, &numberSmoothIntsModC, smoothIntsModC, maxNumberResults);
			currentStart += size - maxRoot;
		}
	}
	finalNewLong = time(NULL);

//print results on the screen
	#ifdef COMPARE
		printf("running times:\n");
		printf("pre-computation (new): %f s\n", (double)(finalPre-initialPre)/CLOCKS_PER_SEC);
		printf("single Sieve: %.3f ms\n", (double)(finalSieve-initialSieve)/CLOCKS_PER_SEC*1000);
		printf("single Log Sieve: %.3f ms\n", (double)(finalLogSieve-initialLogSieve)/CLOCKS_PER_SEC*1000);
		printf("single Trunc Log Sieve: %.3f ms\n", (double)(finalTruncLogSieve-initialTruncLogSieve)/CLOCKS_PER_SEC*1000);
		printf("single Power Trunc Log Sieve: %.3f ms\n", (double)(finalTruncLogSieve2-initialTruncLogSieve2)/CLOCKS_PER_SEC*1000);
		printf("single Naive: %.3f ms\n", (double)(finalNaive-initialNaive)/CLOCKS_PER_SEC*1000);
		printf("single New: %.3f ms\n", (double)(finalNew-initialNew)/CLOCKS_PER_SEC*1000);
		printf("total Naive: %lu s\n", finalNaiveLong-initialNaiveLong);
		printf("total New: %lu s\n", finalNewLong-initialNewLong);
	#endif
	printf("logSieve:\n");
	printf("tolerance to include all smooth integers: %u\n", tolerance);
	printf("roughly %u more smooth integers per chunk\n", surplusSmooth);
	#ifdef COMPARE
	printf("not truncated log sieve: tolerance: %u, surplus: %u\n", tolerance0, surplusSmooth0);
	printf("power truncated log siev: tolerance: %u, surplus: %u\n", tolerancePower, surplusSmoothPower);
	#endif
	printf("parameters:\n");
	printf("sBound: %.2f, start: %.2f, end: %.2f, size: %.2f\n", log2(smoothnessBound), log2(start), log2(end), log2(size));
	printf("degree: %u, C: %u\n", degree, C);
	printf("a(x)=");
	for (unsigned short i=0; i<degree; i++) {
		printf("(x-%d)", polyA[i]);
	}printf("\n");
	printf("b(x)=");
	for (unsigned short i=0; i<degree; i++) {
		printf("(x-%d)", polyB[i]);
	}printf("\n");
	printf("results:\n");
	printf("smoothPrimes: %u", numberSmoothPrimes);
	printf(", residues: %u", numberResidues);
	printf(", relevantSteps: %u\n", relevantSteps);
	printf("smoothIntegersModC: %u\n", numberSmoothIntsModC);
	if (numberSmoothIntsModC>0) {
		printf("Elements l in I s.t. a(l)/C,b(l)/C are smooth integers:\n");
		for (unsigned int i=0; i<numberSmoothIntsModC; i++) {
			printf("%lu\n", smoothIntsModC[i]);
		}
	}

//free allocated memory
	free(smoothPrimes);
	free(logSmoothPrimes);
	free(maxExponents);
	free(minExponents);
	free(smoothInterval);
	free(smoothNumbers);
	free(residues);
	free(smoothIntsModC);
	
	return 0;
}
