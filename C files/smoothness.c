#include"smoothness.h"

//find primes below smoothness bound
void findSmoothPrimes (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes) {
	
	unsigned int i;		//has to be same type as smoothnessBound
	unsigned long j;	//has to be of size at least 2 * smoothnessBound +1
	unsigned int *allNumbers = (unsigned int *) calloc(smoothnessBound, sizeof(int));
	if (allNumbers == NULL) {
		printf("Error: allNumbers could not be allocated\n");
		exit(-1);
	}

	*numberSmoothPrimes = 0;
	for (i=1; i<smoothnessBound; i++) {		//this is a sieve of Eratosthenes
		if (allNumbers[i] == 0) {
			(*numberSmoothPrimes)++;
			for (j=(unsigned long)2*i+1; j<smoothnessBound; j+=i+1) {
				allNumbers[j] = 1;
			}
		}
	}
	*smoothPrimes = (unsigned int *) malloc(*numberSmoothPrimes * sizeof(int));
	if (*smoothPrimes == NULL) {
                printf("Error: smoothPrimes could not be allocated\n");
                exit(-1);
        }
	j = 0;
	for (i=1; i<smoothnessBound; i++) { //transfer it to a shorter list
		if (allNumbers[i] == 0){
			(*smoothPrimes)[j]=i+1;
			j++;
		}
	}

	free(allNumbers);
}

//find highest power for each prime to be in cluded in the sieving for smooth integers
//define kind of smoothness
//powersmooth: maxSizePrimePower = smoothnessBound (restrictive)
//smooth: maxSizePrimePower = end (inefficient)
//compromise: e.g. maxSizePrimePower = size (for small smoothness bounds)
void findMaxExponents (unsigned long maxSizePrimePower, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short **maxExponents) {

	unsigned int i;
	unsigned short j;
	unsigned long q;
	*maxExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	if (maxExponents == NULL) {
		printf("Error: maxExponents could not be allocated\n");
		exit(-1);
	}

	for (i=0; i<numberSmoothPrimes; i++) {		//find max e with p^e < bound
		j = 1;
		q = smoothPrimes[i] * smoothPrimes[i];
		while (q < maxSizePrimePower) {
			j++;
			q *= smoothPrimes[i];
		}
		(*maxExponents)[i] = j;
	}
}

//calculate the rounded logarithm of all primes below the smoothness bound
void findLogSmoothPrimes (unsigned long *logTable, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *logSmoothPrimes) {

	unsigned int i;
	unsigned short l = 0;

	for (i=0; i<numberSmoothPrimes; i++) {
		while (smoothPrimes[i] > logTable[l]) {
			l++;		//smoothness bound << 2^64, so l<64
		}
		logSmoothPrimes[i] = l;
	}
}

//regular sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned long start, unsigned int size, unsigned long *interval, unsigned short *smoothNumbers) {
	
	unsigned int i;
	unsigned int e;
	unsigned int q;
	unsigned int step;

	if (start + size < start) {
		printf("Error: interval exceeds 64Bit values (findSmoothNumbers)\n");
		exit(-1);
	}
	for (i=0; i<size; i++) {
		interval[i] = 1;
	}
	for (i=0; i<numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		q = 1;
		for (e=1; e<=maxExponents[i]; e++) {
			q *= smoothPrimes[i];
			step = q - (start % q);
			if (step == q) {
				step = 0;
			}
			while (step < size) {
				interval[step] *= smoothPrimes[i];
				step += q;
			}
		}
	}
	for (i=0; i<size; i++) {		//translate it into 1 "smooth" or 0 "not smooth"
		if (start + i == interval[i]) {
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
	}
}

//truncated log sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned long *logTable, unsigned long start, unsigned int size, unsigned short *smoothNumbers, unsigned short truncation, unsigned short tolerance) {

	unsigned int i;
	unsigned int e;
	unsigned int q;
	unsigned int step;

	if (start + size < start) {
		printf("Error: interval exceeds 64Bit values (findTruncLogSmoothNumbers)\n");
		exit(-1);
	}
	for (i=0; i<size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i=truncation; i<numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		q = 1;
		for (e=1; e<=maxExponents[i]; e++) {
			q *= smoothPrimes[i];
			step = q - (start % q);
			if (step == q) {
				step = 0;
			}
			while (step < size) {
				smoothNumbers[step] += logSmoothPrimes[i];
				step += q;
			}
		}
	}
	step = round(log2(start));
	for (i=0; i<size; i++) {		//find rounded log2 of l = start + i
		if (start + i > logTable[step]) {
			step++;
			if (step == 64) {
				printf("Warning: log(start + i) > 63.5 (i = %u). Remaining elements are set to 64\n", i);
				while (i < size) {
					if (64 < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
						smoothNumbers[i] = 1;
					}else{
						smoothNumbers[i] = 0;
					}
					i++;
				}
				break;
			}
		}
		if (step < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
	}
}

//truncated log sieve 2
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findTruncLogSmoothNumbers2 (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned long *logTable, unsigned long start, unsigned int size, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short tolerance) {

	unsigned int i;
	unsigned int e;
	unsigned int q;
	unsigned int step;

	if (start + size < start) {
		printf("Error: interval exceeds 64Bit values (findTruncLogSmoothNumbers2)\n");
		exit(-1);
	}
	for (i=0; i<size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i=1; i<numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		q = 1;
		for (e=minExponents[i]; e<=maxExponents[i]; e++) {
			q *= smoothPrimes[i];
			step = q - (start % q);
			if (step == q) {
				step = 0;
			}
			while (step < size) {
				smoothNumbers[step] += logSmoothPrimes[i];
				step += q;
			}
		}
	}
	step = round(log2(start));
	for (i=0; i<size; i++) {	//find rounded log2 of l=start+i
		if (start + i > logTable[step]) {
			step++;
			if (step == 64) {
				printf("Warning: log(start + i) > 63.5 (i = %u). Remaining elements are set to 64\n", i);
				while (i < size) {
					if (64 < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
						smoothNumbers[i] = 1;
					}else{
						smoothNumbers[i] = 0;
					}
					i++;
				}
				break;
			}
		}
		if (step < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
	}
}

//find tolerance s.t. the truncated logarithmic sieving does not exclude any smooth integers
void findTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned long *logTable, unsigned short *logSmoothPrimes, unsigned long start, unsigned int size, unsigned long *smoothInterval, unsigned short *smoothNumbers, unsigned long *smoothIntsModC, unsigned short truncation, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i = 0;
	unsigned short *compareSmoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	if (compareSmoothNumbers == NULL) {
		printf("Error: memory for preSmoothness could not be allocated\n");
		exit(-1);
	}
	*tolerance = 0;
	*surplusSmooth = 0;
	findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, start, size, smoothInterval, smoothNumbers);
	findTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncation, *tolerance);
	while (i < size) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			(*tolerance)++;
			findTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncation, *tolerance);
		}else{
			i++;
		}
	}
	(*tolerance)++;
	findTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncation, *tolerance);
	
	for (i=0; i<size; i++) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			printf("logSieve misses smooth value: %lu\n", start+i);
		}
		if (compareSmoothNumbers[i] > smoothNumbers[i]) {
			(*surplusSmooth)++;
		}
	}
	free(compareSmoothNumbers);
}

//pre-computation for smoothness sieving
void preSmoothness (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes, unsigned long maxSizePrimePower, unsigned short **maxExponents, unsigned short numberRoots, short *roots, unsigned short *maxRoot, unsigned long **logTable, unsigned short **logSmoothPrimes, unsigned long start, unsigned int size, unsigned long **smoothInterval, unsigned short **smoothNumbers, unsigned short maxNumberResults, unsigned long **smoothIntsModC, unsigned short truncation, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i;
	
	findSmoothPrimes(smoothnessBound, numberSmoothPrimes, smoothPrimes);
	findMaxExponents(maxSizePrimePower, *numberSmoothPrimes, *smoothPrimes, maxExponents);
	
	*maxRoot = 0;
	for (i=1; i<numberRoots; i++) {
		if (roots[i] > *maxRoot) {
			*maxRoot = roots[i];
		}
	}
	static unsigned long logTable64[64];
	for (i=0; i<64; i++) {
		logTable64[i] = floor(pow(2, i + 0.5));
	}
	*logTable = logTable64;
	*logSmoothPrimes = (unsigned short *) malloc(*numberSmoothPrimes * sizeof(short));
	*smoothInterval = (unsigned long *) malloc(size * sizeof(long));
	*smoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	*smoothIntsModC = (unsigned long *) calloc(maxNumberResults, sizeof(long));
	
	if (*logSmoothPrimes == NULL || *smoothInterval == NULL || *smoothNumbers == NULL || *smoothIntsModC == NULL) {
		printf("Error: memory for preSmoothness could not be allocated\n");
		exit(-1);
	}
	findLogSmoothPrimes(*logTable, *numberSmoothPrimes, *smoothPrimes, *logSmoothPrimes);
	
	findTolerance(*numberSmoothPrimes, *smoothPrimes, *maxExponents, *logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, *smoothIntsModC, truncation, tolerance, surplusSmooth);
}
