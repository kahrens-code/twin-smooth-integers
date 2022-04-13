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
	for (i=1; i<smoothnessBound; i++) { 	//transfer it to a shorter list
		if (allNumbers[i] == 0) {
			(*smoothPrimes)[j] = i + 1;
			j++;
		}
	}

	free(allNumbers);
}

//find highest power for each prime to be included in the sieving for smooth integers
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

//find powers for each prime to be in cluded in the sieving for smooth integers
//define kind of smoothness
//powersmooth: maxSizePrimePower = smoothnessBound (restrictive)
//smooth: maxSizePrimePower = end (inefficient, possibly too large for unsigned long)
//compromise: e.g. maxSizePrimePower = size (for small smoothness bounds)
void findExponents (unsigned long maxSizePrimePower, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short truncPower, unsigned short **maxExponents, unsigned short **minExponents) {

	unsigned int i;
	unsigned short j;
	unsigned long q;
	unsigned long threshold;
	*maxExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	*minExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	if (*maxExponents == NULL || *minExponents == NULL) {
		printf("Error: max/min Exponents could not be allocated\n");
		exit(-1);
	}

	threshold = pow2(truncPower);
	for (i=0; i<numberSmoothPrimes; i++) {		//find e with theshold <= p^e < maxSizePrimePower
		j = 1;
		q = smoothPrimes[i];
		while (q < threshold) {
			j++;
			q *= smoothPrimes[i];
		}
		(*minExponents)[i] = j;
		q *= smoothPrimes[i];
		while (q < maxSizePrimePower) {
			j++;
			q *= smoothPrimes[i];
		}
		(*maxExponents)[i] = j;
	}
}

//calculate the rounded logarithm of all primes below the smoothness bound
void findLogSmoothPrimes (unsigned __int128 *logTable, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short **logSmoothPrimes) {

	unsigned int i;
	unsigned short l = 0;

	*logSmoothPrimes = (unsigned short *) malloc(numberSmoothPrimes * sizeof(short));
	if (*logSmoothPrimes == NULL) {
		printf("Error: logSmoothPrimes could not be allocated\n");
		exit(-1);
	}

	for (i=0; i<numberSmoothPrimes; i++) {
		while (smoothPrimes[i] > logTable[l]) {
			l++;		//smoothness bound << 2^64, so l<64
		}
		(*logSmoothPrimes)[i] = l;
	}
}

//regular sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned __int128 start, unsigned int size, unsigned __int128 *interval, unsigned short *smoothNumbers) {
	
	unsigned int i;
	unsigned short e;
	unsigned long q;
	unsigned long step;

	if (start + size < start) {
		printf("Error: interval exceeds 128bit values (findSmoothNumbers)\n");
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
void findPrimeTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned __int128 *logTable, unsigned __int128 start, unsigned int size, unsigned short *smoothNumbers, unsigned short truncPrime, unsigned short tolerance) {

	unsigned int i;
	unsigned short e;
	unsigned long q;
	unsigned long step;

	if (start + size < start) {
		printf("Error: interval exceeds 128bit values (findPrimeTruncLogSmoothNumbers)\n");
		exit(-1);
	}
	for (i=0; i<size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i=truncPrime; i<numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
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
	step = round(log2(start)) - 1;		//reduce by 1 to take care of rounding error in __int128 to double conversion
	for (i=0; i<size; i++) {		//find rounded log2 of l = start + i
		if (start + i > logTable[step]) {
			step++;
			if (step == 128) {
				printf("Warning: log2(start + i) > 127.5 (i = %u). Remaining elements are set to 128\n", i);
				while (i < size) {
					if (128 < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
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

//power truncated log sieve
//CAUTION: intervals with elements beyond 2^64 cause overflow and have to be tackled differently
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findPowerTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, unsigned __int128 *logTable, unsigned __int128 start, unsigned int size, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short tolerance) {

	unsigned int i;
	unsigned short e;
	unsigned long q;
	unsigned long step;

	if (start + size < start) {
		printf("Error: interval exceeds 128bit values (findPowerTruncLogSmoothNumbers)\n");
		exit(-1);
	}
	for (i=0; i<size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i=0; i<numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
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
	step = round(log2(start)) - 1;		//reduce by 1 to take care of rounding error in __int128 to double conversion
	for (i=0; i<size; i++) {	//find rounded log2 of l=start+i
		if (start + i > logTable[step]) {
			step++;
			if (step == 128) {
				printf("Warning: log2(start + i) > 127.5 (i = %u). Remaining elements are set to 128\n", i);
				while (i < size) {
					if (128 < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
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

//find tolerance s.t. the prime truncated logarithmic sieving does not exclude any smooth integers
void findTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned __int128 *logTable, unsigned short *logSmoothPrimes, unsigned __int128 start, unsigned int size, unsigned __int128 *smoothInterval, unsigned short *smoothNumbers, unsigned short truncPrime, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i = 0;
	unsigned short *compareSmoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	if (compareSmoothNumbers == NULL) {
		printf("Error: memory for findTolerance could not be allocated\n");
		exit(-1);
	}
	*tolerance = 0;
	*surplusSmooth = 0;
	findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, start, size, smoothInterval, smoothNumbers);
	findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncPrime, *tolerance);
	while (i < size) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			(*tolerance)++;
			findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncPrime, *tolerance);
		}else{
			i++;
		}
	}
	(*tolerance)++;
	findPrimeTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, truncPrime, *tolerance);
	
	for (i=0; i<size; i++) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			printf("logSieve misses smooth value: ");
			printu128(start + i); printf("\n");
		}
		if (compareSmoothNumbers[i] > smoothNumbers[i]) {
			(*surplusSmooth)++;
		}
	}
	free(compareSmoothNumbers);
}

//find tolerance s.t. the power truncated logarithmic sieving does not exclude any smooth integers
void findPowerTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned __int128 *logTable, unsigned short *logSmoothPrimes, unsigned __int128 start, unsigned int size, unsigned __int128 *smoothInterval, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i = 0;
	unsigned short *compareSmoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	if (compareSmoothNumbers == NULL) {
		printf("Error: memory for findPowerTolerance could not be allocated\n");
		exit(-1);
	}
	*tolerance = 0;
	*surplusSmooth = 0;
	findSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, start, size, smoothInterval, smoothNumbers);
	findPowerTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, minExponents, *tolerance);
	while (i < size) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			(*tolerance)++;
			findPowerTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, minExponents, *tolerance);
		}else{
			i++;
		}
	}
	(*tolerance)++;
	findPowerTruncLogSmoothNumbers(numberSmoothPrimes, smoothPrimes, maxExponents, logSmoothPrimes, logTable, start, size, compareSmoothNumbers, minExponents, *tolerance);
	
	for (i=0; i<size; i++) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			printf("logSieve misses smooth value: ");
			printu128(start + i); printf("\n");
		}
		if (compareSmoothNumbers[i] > smoothNumbers[i]) {
			(*surplusSmooth)++;
		}
	}
	free(compareSmoothNumbers);
}

//pre-computation for smoothness sieving
void preSmoothness (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes, unsigned long maxSizePrimePower, unsigned short **maxExponents, unsigned short numberRoots, short *roots, unsigned short *maxRoot, unsigned __int128 **logTable, unsigned short **logSmoothPrimes, unsigned __int128 start, unsigned int size, unsigned __int128 **smoothInterval, unsigned short **smoothNumbers, unsigned short maxNumberResults, unsigned __int128 **smoothIntsModC, unsigned short truncPrime, unsigned short truncPower, unsigned short **minExponents, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	short i;
	unsigned int j;
	static unsigned __int128 logTable128[128];
	*smoothInterval = (unsigned __int128 *) malloc(size * sizeof(__int128));
	*smoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	*smoothIntsModC = (unsigned __int128 *) calloc(maxNumberResults, sizeof(__int128));
	if (*smoothInterval == NULL || *smoothNumbers == NULL || *smoothIntsModC == NULL) {
		printf("Error: memory for preSmoothness could not be allocated\n");
		exit(-1);
	}
	unsigned __int128 e19 = 10000000000000000000u; //10^19
	unsigned __int128 e38 = e19*e19;
	logTable128[127] = e38 * 2 + e19 * 4061596916800451154 + 5033772477625056927; //floor(2^127.5)
	for (i=126; i>=0; i--) {
		logTable128[i] = logTable128[i+1] / 2;
	}
	*logTable = logTable128;
	*maxRoot = 0;
	for (j=1; j<numberRoots; j++) {
		if (roots[j] > *maxRoot) {
			*maxRoot = roots[j];
		}
	}
	findSmoothPrimes(smoothnessBound, numberSmoothPrimes, smoothPrimes);
	findLogSmoothPrimes(*logTable, *numberSmoothPrimes, *smoothPrimes, logSmoothPrimes);
	if (pow2(truncPower) >= maxSizePrimePower || truncPrime >= *numberSmoothPrimes) {
		printf("Error: too much truncation, nothing left\n");
		exit(-1);
	}
	findExponents(maxSizePrimePower, *numberSmoothPrimes, *smoothPrimes, truncPower, maxExponents, minExponents);
	findTolerance(*numberSmoothPrimes, *smoothPrimes, *maxExponents, *logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, truncPrime, tolerance, surplusSmooth);
}
