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
	for (i = 1; i < smoothnessBound; i++) {		//this is a sieve of Eratosthenes
		if (allNumbers[i] == 0) {
			(*numberSmoothPrimes)++;
			for (j = (unsigned long) 2 * i + 1; j < smoothnessBound; j += i + 1) {
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
	for (i = 1; i < smoothnessBound; i++) { 	//transfer it to a shorter list
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
void findMaxExponents (mpz_t maxSizePrimePower, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short **maxExponents) {

	unsigned int i;
	unsigned short j;
	mpz_t q; mpz_init(q);
	*maxExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	if (maxExponents == NULL) {
		printf("Error: maxExponents could not be allocated\n");
		exit(-1);
	}

	for (i = 0; i < numberSmoothPrimes; i++) {		//find max e with p^e < bound
		j = 0;
		mpz_set_ui(q, smoothPrimes[i]);
		while (mpz_cmp(q, maxSizePrimePower) < 0) {
			j++;
			mpz_mul_ui(q, q, smoothPrimes[i]);
		}
		(*maxExponents)[i] = j;
	}
	
	mpz_clear(q);
}

//find powers for each prime to be included in the sieving for smooth integers
//define kind of smoothness
//powersmooth: maxSizePrimePower = smoothnessBound (restrictive)
//smooth: maxSizePrimePower = end (inefficient)
//compromise: e.g. maxSizePrimePower = size (for small smoothness bounds)
void findExponents (mpz_t maxSizePrimePower, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short truncPower, unsigned short **maxExponents, unsigned short **minExponents) {

	unsigned int i;
	unsigned short j;
	mpz_t q; mpz_init(q);
	mpz_t threshold; mpz_init_set_ui(threshold, 1);
	*maxExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	*minExponents = (unsigned short *) calloc(numberSmoothPrimes, sizeof(short));
	if (*maxExponents == NULL || *minExponents == NULL) {
		printf("Error: max/min Exponents could not be allocated\n");
		exit(-1);
	}

	mpz_mul_2exp(threshold, threshold, truncPower);
	for (i = 0; i < numberSmoothPrimes; i++) {	//find e with theshold <= p^e < maxSizePrimePower
		j = 1;
		mpz_set_ui(q, smoothPrimes[i]);
		while (mpz_cmp(q, threshold) < 0) {
			j++;
			mpz_mul_ui(q, q, smoothPrimes[i]);
		}
		(*minExponents)[i] = j;
		mpz_mul_ui(q, q, smoothPrimes[i]);
		while (mpz_cmp(q, maxSizePrimePower) < 0) {
			j++;
			mpz_mul_ui(q, q, smoothPrimes[i]);
		}
		(*maxExponents)[i] = j;
	}
	
	mpz_clear(q);
	mpz_clear(threshold);
}

//calculate the rounded logarithm of all primes below the smoothness bound
void findLogSmoothPrimes (mpz_t *logTable, unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short **logSmoothPrimes) {

	unsigned int i;
	unsigned short l = 0;

	*logSmoothPrimes = (unsigned short *) malloc(numberSmoothPrimes * sizeof(short));
	if (*logSmoothPrimes == NULL) {
		printf("Error: logSmoothPrimes could not be allocated\n");
		exit(-1);
	}

	for (i=0; i<numberSmoothPrimes; i++) {
		while (mpz_cmp_ui(logTable[l], smoothPrimes[i]) < 0) {
			l++;		//smoothness bound << 2^64, so l<64
		}
		(*logSmoothPrimes)[i] = l;
	}
}

//regular sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t start, unsigned int size, mpz_t *interval, unsigned short *smoothNumbers) {
	
	unsigned int i;
	unsigned short e;
	mpz_t q; mpz_init(q);
	mpz_t step; mpz_init(step);

	for (i = 0; i < size; i++) {
		mpz_set_ui(interval[i], 1);
	}
	for (i = 0; i < numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		mpz_set_ui(q, 1);
		for (e = 1; e <= maxExponents[i]; e++) {
			mpz_mul_ui(q, q, smoothPrimes[i]);
			mpz_mod(step, start, q);
			mpz_sub(step, q, step);	//step = q - (start % q);
			if (mpz_cmp(step, q) == 0) {
				mpz_set_ui(step, 0);
			}
			while (mpz_cmp_ui(step, size) < 0) {
				mpz_mul_ui(interval[mpz_get_ui(step)], interval[mpz_get_ui(step)], smoothPrimes[i]);
				mpz_add(step, step, q);
			}
		}
	}
	mpz_set(step, start);
	for (i = 0; i < size; i++) {		//translate it into 1 "smooth" or 0 "not smooth"
		if (mpz_cmp(step, interval[i]) == 0) {
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
		mpz_add_ui(step, step, 1);
	}
	
	mpz_clear(q);
	mpz_clear(step);
}

//truncated log sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findPrimeTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, mpz_t *logTable, mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short truncPrime, unsigned short tolerance) {

	unsigned int i;
	unsigned short e;
	mpz_t q; mpz_init(q);
	mpz_t step; mpz_init(step);
	unsigned int logStep;

	for (i = 0; i < size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i = truncPrime; i < numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		mpz_set_ui(q, 1);
		for (e = 1; e <= maxExponents[i]; e++) {
			mpz_mul_ui(q, q, smoothPrimes[i]);
			mpz_mod(step, start, q);
			mpz_sub(step, q, step);	//step = q - (start % q);
			if (mpz_cmp(step, q) == 0) {
				mpz_set_ui(step, 0);
			}
			while (mpz_cmp_ui(step, size) < 0) {
				smoothNumbers[mpz_get_ui(step)] += logSmoothPrimes[i];
				mpz_add(step, step, q);
			}
		}
	}
	logStep = (unsigned int) mpz_sizeinbase(start, 2) - 1;		//reduce by 1 to take care of surplus in sizeinbase
	for (i = 0; i < size; i++) {		//find rounded log2 of l = start + i
		mpz_add_ui(step, start, i);
		while (mpz_cmp(step, logTable[logStep]) > 0) {
			logStep++;
			if (logStep == 128) {
				printf("Warning: log2(start + %u) > %f. Remaining elements are set to %u\n", i, 127.5, 128);
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
		if (logStep < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
	}
	
	mpz_clear(q);
	mpz_clear(step);
}

//power truncated log sieve
//CAUTION: just returns an array of 0s and 1s, start value is NOT included
void findPowerTruncLogSmoothNumbers (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, unsigned short *logSmoothPrimes, mpz_t *logTable, mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short tolerance) {

	unsigned int i;
	unsigned short e;
	mpz_t q; mpz_init(q);
	mpz_t step; mpz_init(step);
	unsigned int logStep;

	for (i = 0; i < size; i++) {
		smoothNumbers[i] = 0;
	}
	for (i = 0; i < numberSmoothPrimes; i++) {		//use a sieve to find smooth elements
		mpz_set_ui(q, 1);
		for (e = minExponents[i]; e <= maxExponents[i]; e++) {
			mpz_mul_ui(q, q, smoothPrimes[i]);
			mpz_mod(step, start, q);
			mpz_sub(step, q, step);	//step = q - (start % q);
			if (mpz_cmp(step, q) == 0) {
				mpz_set_ui(step, 0);
			}
			while (mpz_cmp_ui(step, size) < 0) {
				smoothNumbers[mpz_get_ui(step)] += logSmoothPrimes[i];
				mpz_add(step, step, q);
			}
		}
	}
	logStep = (unsigned int) mpz_sizeinbase(start, 2) - 1;			//reduce by 1 to take care of surplus in sizeinbase
	for (i = 0; i < size; i++) {	//find rounded log2 of l=start+i
		mpz_add_ui(step, start, i);
		while (mpz_cmp(step, logTable[logStep]) > 0) {
			logStep++;
			if (logStep == 128) {
				printf("Warning: log2(start + %u) > %f. Remaining elements are set to %u\n", i, 127.5, 128);
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
		if (logStep < smoothNumbers[i] + tolerance) {   //translate into 1 "smooth" or 0 "not smooth"
			smoothNumbers[i] = 1;
		}else{
			smoothNumbers[i] = 0;
		}
	}
	
	mpz_clear(q);
	mpz_clear(step);
}

//find tolerance s.t. the prime truncated logarithmic sieving does not exclude any smooth integers
void findTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t *logTable, unsigned short *logSmoothPrimes, mpz_t start, unsigned int size, mpz_t *smoothInterval, unsigned short *smoothNumbers, unsigned short truncPrime, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i = 0;
	mpz_t step; mpz_init(step);
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
	
	for (i = 0; i < size; i++) {
		if (compareSmoothNumbers[i] < smoothNumbers[i]) {
			mpz_add_ui(step, start, i);
			gmp_printf("logSieve misses smooth value: %Zd\n", step);
		}
		if (compareSmoothNumbers[i] > smoothNumbers[i]) {
			(*surplusSmooth)++;
		}
	}
	
	mpz_clear(step);
	free(compareSmoothNumbers);
}

//find tolerance s.t. the power truncated logarithmic sieving does not exclude any smooth integers
void findPowerTolerance (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short *maxExponents, mpz_t *logTable, unsigned short *logSmoothPrimes, mpz_t start, unsigned int size, mpz_t *smoothInterval, unsigned short *smoothNumbers, unsigned short *minExponents, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	unsigned int i = 0;
	mpz_t step; mpz_init(step);
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
			mpz_add_ui(step, start, i);
			gmp_printf("logSieve misses smooth value: %Zd\n", step);
		}
		if (compareSmoothNumbers[i] > smoothNumbers[i]) {
			(*surplusSmooth)++;
		}
	}
	
	mpz_clear(step);
	free(compareSmoothNumbers);
}

//general pre-computation for smoothness
void preSmoothness (unsigned int smoothnessBound, unsigned int *numberSmoothPrimes, unsigned int **smoothPrimes, mpz_t **logTable, unsigned short maxNumberResults, mpz_t **smoothIntsModC, unsigned short numberRoots, short *roots, unsigned short *maxRoot) {
	
	int i;
	*logTable = (mpz_t *) malloc(128 * sizeof(mpz_t));
	*smoothIntsModC = (mpz_t *) calloc(maxNumberResults, sizeof(mpz_t));
	if (*logTable == NULL || *smoothIntsModC == NULL) {
		printf("Error: memory for preSmoothness could not be allocated\n");
		exit(-1);
	}
	mpz_init_set_str((*logTable)[127], "240615969168004511545033772477625056927", 10);	//floor(2^127.5)
	for (i = 126; i >= 0; i--) {
		mpz_init((*logTable)[i]);
		mpz_fdiv_q_2exp((*logTable)[i], (*logTable)[i + 1], 1);
	}
	for (i = 0; i < maxNumberResults; i++) {
		mpz_init((*smoothIntsModC)[i]);
	}
	*maxRoot = 0;
	for (i = 1; i < numberRoots; i++) {
		if (roots[i] > *maxRoot) {
			*maxRoot = roots[i];
		}
	}
	findSmoothPrimes(smoothnessBound, numberSmoothPrimes, smoothPrimes);
}

//pre-computation for smoothness sieving
void preSievingCompare (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t maxSizePrimePower, unsigned short **maxExponents, mpz_t *logTable, unsigned short **logSmoothPrimes, mpz_t start, unsigned int size, mpz_t **smoothInterval, unsigned short **smoothNumbers, unsigned short truncPrime, unsigned short truncPower, unsigned short **minExponents, unsigned short *tolerance0, unsigned short *tolerance, unsigned short *tolerancePower, unsigned int *surplusSmooth0, unsigned int *surplusSmooth, unsigned int *surplusSmoothPower) {
	
	*smoothInterval = (mpz_t *) malloc(size * sizeof(mpz_t));
	*smoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	if (*smoothInterval == NULL || *smoothNumbers == NULL) {
		printf("Error: memory for preSieving could not be allocated\n");
		exit(-1);
	}
	mpz_t truncStart; mpz_init_set_ui(truncStart, 1);
	mpz_mul_2exp(truncStart, truncStart, (unsigned int) truncPower);
	if (mpz_cmp(truncStart, maxSizePrimePower) >= 0 || truncPrime >= numberSmoothPrimes) {
		printf("Error: too much truncation, nothing left\n");
		exit(-1);
	}
	mpz_clear(truncStart);
	findLogSmoothPrimes(logTable, numberSmoothPrimes, smoothPrimes, logSmoothPrimes);
	findExponents(maxSizePrimePower, numberSmoothPrimes, smoothPrimes, truncPower, maxExponents, minExponents);
	findTolerance(numberSmoothPrimes, smoothPrimes, *maxExponents, logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, 0, tolerance0, surplusSmooth0);
	findTolerance(numberSmoothPrimes, smoothPrimes, *maxExponents, logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, truncPrime, tolerance, surplusSmooth);
	findPowerTolerance (numberSmoothPrimes, smoothPrimes, *maxExponents, logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, *minExponents, tolerancePower, surplusSmoothPower);
}

//pre-computation for smoothness prime truncated log sieving
void preSievingPrimeTruncLog (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t maxSizePrimePower, unsigned short **maxExponents, mpz_t *logTable, unsigned short **logSmoothPrimes, mpz_t start, unsigned int size, mpz_t **smoothInterval, unsigned short **smoothNumbers, unsigned short truncPrime, unsigned short *tolerance, unsigned int *surplusSmooth) {
	
	*smoothInterval = (mpz_t *) malloc(size * sizeof(mpz_t));
	*smoothNumbers = (unsigned short *) malloc(size * sizeof(short));
	if (*smoothInterval == NULL || *smoothNumbers == NULL) {
		printf("Error: memory for preSieving could not be allocated\n");
		exit(-1);
	}
	if (truncPrime >= numberSmoothPrimes) {
		printf("Error: too much truncation, nothing left\n");
		exit(-1);
	}
	unsigned short truncPower = 10;
	unsigned short *minExponents;
	findLogSmoothPrimes(logTable, numberSmoothPrimes, smoothPrimes, logSmoothPrimes);
	findExponents(maxSizePrimePower, numberSmoothPrimes, smoothPrimes, truncPower, maxExponents, &minExponents);
	findTolerance(numberSmoothPrimes, smoothPrimes, *maxExponents, logTable, *logSmoothPrimes, start, size, *smoothInterval, *smoothNumbers, truncPrime, tolerance, surplusSmooth);
	free(minExponents);
}
