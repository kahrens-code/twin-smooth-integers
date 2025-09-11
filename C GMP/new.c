#include"new.h"

//compare function for sorting
int cmpfnc (const void *a, const void *b) {
	return mpz_cmp(*(mpz_t *)a, *(mpz_t *)b);
}

//Chinese remainder theorem
void CRT (
	unsigned int number,	//number of different moduli
	mpz_t *values, 	//array of residues of each modulus
	mpz_t *n,	//array of moduli
	mpz_t NN,	//product of moduli
	mpz_t *result	//pointer to set to output of CRT
	) {
	
	mpz_t *N = (mpz_t *) malloc(number * sizeof(mpz_t)); 
	mpz_t *M = (mpz_t *) malloc(number * sizeof(mpz_t)); 
	if (N == NULL || M == NULL) {
		printf("Error: CRT arrays could not be allocated\n");
		exit(-1);
	}
	unsigned int i;
	mpz_t r, oldR; mpz_inits(r, oldR, NULL);
	mpz_t s, oldS; mpz_inits(s, oldS, NULL);
	mpz_t t, oldT; mpz_inits(t, oldT, NULL);
	mpz_t quot; mpz_init(quot);
	mpz_t y; mpz_init_set_ui(y, 0);

	for (i = 0; i < number; i++) {	//find M[i] s.t. M[i] * (NN / n[i]) = 1 mod n[i]
		mpz_inits(N[i], M[i], NULL);
		mpz_divexact(N[i], NN, n[i]);
		mpz_set(oldR, N[i]);
		mpz_set(r, n[i]);
		mpz_set_ui(oldS, 1);
		mpz_set_ui(s, 0);
		mpz_set_ui(oldT, 0);
		mpz_set_ui(t, 1);
		while (mpz_cmp_ui(r, 0) != 0) {	//extended Euclidean algorithm with output s * N[i] + t * n[i] = gcd(N[i], n[i])
			mpz_fdiv_q(quot, oldR, r);
			mpz_swap(r, oldR);	//temp = oldR; oldR = r; r = temp;
			mpz_submul(r, quot, oldR);	//r = r - quot * oldR;
			mpz_swap(s, oldS);	//temp = oldS; oldS = s; s = temp;
			mpz_submul(s, quot, oldS);	//s = s - quot * oldS;
			mpz_swap(t, oldT);	//temp = oldT; oldT = t; t = temp;
			mpz_submul(t, quot, oldT);	//t = t - quot * oldT;
		}
		mpz_set(M[i], oldS);
	}
	for (i = 0; i < number; i++) {
		mpz_mul(quot, values[i], M[i]);
		mpz_mul(quot, quot, N[i]);
		mpz_add(y, y, quot);
		mpz_mod(y, y, NN);	//y = (y + values[i] * M[i] * N[i]) % NN;
		mpz_clears(M[i], N[i], NULL);
	}

	mpz_set(*result, y);
	
	mpz_clears(r, oldR, NULL);
	mpz_clears(s, oldS, NULL);
	mpz_clears(t, oldT, NULL);
	mpz_clear(quot);
	mpz_clear(y);
	free(N);
	free(M);
}

//Chinese remainder theorem
void CRT_ui (
	unsigned int number,	//number of different moduli
	unsigned long *values, 	//array of residues of each modulus
	mpz_t *n,	//array of moduli
	mpz_t NN,	//product of moduli
	mpz_t *result	//pointer to set to output of CRT
	) {
	
	mpz_t *N = (mpz_t *) malloc(number * sizeof(mpz_t)); 
	mpz_t *M = (mpz_t *) malloc(number * sizeof(mpz_t)); 
	if (N == NULL || M == NULL) {
		printf("Error: CRT arrays could not be allocated\n");
		exit(-1);
	}
	unsigned int i;
	mpz_t r, oldR; mpz_inits(r, oldR, NULL);
	mpz_t s, oldS; mpz_inits(s, oldS, NULL);
	mpz_t t, oldT; mpz_inits(t, oldT, NULL);
	mpz_t quot; mpz_init(quot);
	mpz_t y; mpz_init_set_ui(y, 0);

	for (i = 0; i < number; i++) {	//find M[i] s.t. M[i] * (NN / n[i]) = 1 mod n[i]
		mpz_inits(N[i], M[i], NULL);
		mpz_divexact(N[i], NN, n[i]);
		mpz_set(oldR, N[i]);
		mpz_set(r, n[i]);
		mpz_set_ui(oldS, 1);
		mpz_set_ui(s, 0);
		mpz_set_ui(oldT, 0);
		mpz_set_ui(t, 1);
		while (mpz_cmp_ui(r, 0) != 0) {	//extended Euclidean algorithm with output s * N[i] + t * n[i] = gcd(N[i], n[i])
			mpz_fdiv_q(quot, oldR, r);
			mpz_swap(r, oldR);	//temp = oldR; oldR = r; r = temp;
			mpz_submul(r, quot, oldR);	//r = r - quot * oldR;
			mpz_swap(s, oldS);	//temp = oldS; oldS = s; s = temp;
			mpz_submul(s, quot, oldS);	//s = s - quot * oldS;
			mpz_swap(t, oldT);	//temp = oldT; oldT = t; t = temp;
			mpz_submul(t, quot, oldT);	//t = t - quot * oldT;
		}
		mpz_set(M[i], oldS);
	}
	for (i = 0; i < number; i++) {
		mpz_mul_ui(quot, M[i], values[i]);
		mpz_mul(quot, quot, N[i]);
		mpz_add(y, y, quot);
		mpz_mod(y, y, NN);	//y = (y + values[i] * M[i] * N[i]) % NN;
		mpz_clears(M[i], N[i], NULL);
	}

	mpz_set(*result, y);
	
	mpz_clears(r, oldR, NULL);
	mpz_clears(s, oldS, NULL);
	mpz_clears(t, oldT, NULL);
	mpz_clear(quot);
	mpz_clear(y);
	free(N);
	free(M);
}

void factorIntoSmoothPrimes (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, mpz_t *logTable, mpz_t input, mpz_t **factors, unsigned int *numberFactors, mpz_t *factorSum) {
	
//find log2(input) as rough estimate for number of factors
	unsigned int logInput = 0;
	
	while (mpz_cmp(logTable[logInput], input) <= 0) {
		logInput++;
		if (logInput == 128) {
			printf("Warning: log2(input) > %f. Treated as %u.\n", 127.5, 128);
			break;
		}
	}
	
//trial division with smooth primes
	mpz_t *tempFactors = (mpz_t *) malloc(logInput * sizeof(mpz_t));
	if (tempFactors == NULL) {
		printf("Error: tempFactors could not be allocated\n");
		exit(-1);
	}
	unsigned short newFactor = 0;
	unsigned int i = 0;
	mpz_t decompose; mpz_init_set(decompose, input);
	*numberFactors = 0;
	mpz_set_ui(*factorSum, 0);
	mpz_init_set_ui(tempFactors[*numberFactors], 1);
	
	while (mpz_cmp_ui(decompose, 1) > 0) {
		if (mpz_divisible_ui_p(decompose, smoothPrimes[i])) {
			mpz_divexact_ui(decompose, decompose, smoothPrimes[i]);
			mpz_mul_ui(tempFactors[*numberFactors], tempFactors[*numberFactors], smoothPrimes[i]);
			newFactor = 1;
		}else {
			i++;
			if (i == numberSmoothPrimes) {
				printf("Error: input not smooth\n");
				exit(-1);
			}
			if (newFactor == 1) {
				mpz_add(*factorSum, *factorSum, tempFactors[*numberFactors]);
				(*numberFactors)++;
				if (*numberFactors == logInput) {
					printf("Error: input has more than log2(input) prime factors");
					exit(-1);
				}
				mpz_init_set_ui(tempFactors[*numberFactors], 1);
				newFactor = 0;
			}
		}
	}
	(*numberFactors)++;
	
	mpz_clear(decompose);
	
//write factors into a shorter array
	*factors =  (mpz_t *) malloc(*numberFactors * sizeof(mpz_t));
	for (i = 0; i < *numberFactors; i++) {
		mpz_init_set((*factors)[i], tempFactors[i]);
		mpz_clear(tempFactors[i]);
	}
	
	free(tempFactors);
}

//find residue classes such that poly(res) = 0 mod prime power factors of C	
void findPrimeResidues (mpz_t *factors, unsigned int numberFactors, unsigned long maxNumberPrimeResidues, unsigned short degree, short *poly, mpz_t *numberResidues, unsigned long **primeResidues, unsigned long **numbersPrimeResidues, unsigned long **indices) {
	
	unsigned int i;
	unsigned long j;
	unsigned int k;
	unsigned long currentFactor;
	unsigned long totalNumberPrimeResidues;
	mpz_t eval; mpz_init(eval);
	mpz_t temp; mpz_init(temp);
	unsigned long *tempPrimeResidues = (unsigned long *) calloc(maxNumberPrimeResidues, sizeof(long));
	*numbersPrimeResidues = (unsigned long *) calloc(numberFactors, sizeof(long));
	*indices = (unsigned long *) calloc(numberFactors + 1, sizeof(long));
	mpz_t *poly_gmp = (mpz_t *) malloc(degree * sizeof(mpz_t));
	if (tempPrimeResidues == NULL || numbersPrimeResidues == NULL || indices == NULL || poly_gmp == NULL) {
		printf("Error: one of the arrays in findPrimeResidues could not be allocated\n");
		exit(-1);
	}
	for (i = 0; i < degree; i++) {
		mpz_init_set_si(poly_gmp[i], poly[i]);
	}
	totalNumberPrimeResidues = 0;
	mpz_set_ui(*numberResidues, 1);
	for (i = 0; i < numberFactors; i++) {	//find solutions for a(x) = 0 modulo prime power factors of C
		currentFactor = mpz_get_ui(factors[i]);
		for (j = 0; j < currentFactor; j++) {
			mpz_set_ui(eval, 1);
			for (k = 0; k < degree; k++) {	//compute a(j) % currentFactor
				mpz_ui_sub(temp, j, poly_gmp[k]);
				mpz_mul(eval, eval, temp);
				mpz_mod_ui(eval, eval, currentFactor);	//eval = (eval * (j - poly[k])) % currentFactor; use % to keep size reasonable
			}
			if (mpz_cmp_ui(eval, 0) == 0) {
				tempPrimeResidues[totalNumberPrimeResidues] = j;
				totalNumberPrimeResidues++;
				(*numbersPrimeResidues)[i]++;
			}
		}
		mpz_mul_ui(*numberResidues, *numberResidues, (*numbersPrimeResidues)[i]);
		(*indices)[i + 1] = (*indices)[i] + (*numbersPrimeResidues)[i];	//keep track of which elements in primeResidues correspond to which factor; the last entry is the total number of primeResidues
	}
//write primeResidues into a shorter array and clear/free tempPrimeResidues
	*primeResidues = (unsigned long *) malloc((*indices)[numberFactors] * sizeof(long));
	for (i = 0; i < (*indices)[numberFactors]; i++) {
		(*primeResidues)[i] = tempPrimeResidues[i];
	}
	
	for (i = 0; i < degree; i++) {
		mpz_clear(poly_gmp[i]);
	}
	mpz_clears(eval, temp, NULL);
	free(tempPrimeResidues);
	free(poly_gmp);
}

//create ordered list of residue classes such that poly(res) = 0 mod C	 	
void createResidueList (mpz_t *factors, unsigned int numberFactors, unsigned long numberResidues_ui, unsigned long *primeResidues, unsigned long *numbersPrimeResidues, unsigned long *indices, mpz_t C, mpz_t **residues) {
	unsigned long i;
	unsigned long j;
	unsigned long *tupelCRT = (unsigned long *) calloc(numberFactors, sizeof(long));
	unsigned long *currentPrimeResidues = (unsigned long *) malloc(numberFactors * sizeof(long));
	*residues = (mpz_t *) malloc(numberResidues_ui * sizeof(mpz_t));
	if (*residues == NULL || tupelCRT == NULL || currentPrimeResidues == NULL) {
                printf("Error: residue lists could not be allocated\n");
                exit(-1);
        }
	
	for (i = 0; i < numberFactors; i++) {
		currentPrimeResidues[i] = primeResidues[indices[i] + tupelCRT[i]];
	}	
	for (i = 0; i < numberResidues_ui; i++) {	//use CRT to combine solutions mod prime powers to solutions mod C
		mpz_init((*residues)[i]);
		CRT_ui(numberFactors, currentPrimeResidues, factors, C, &((*residues)[i]));
		for (j = 0; j < numberFactors; j++) {
			if (tupelCRT[j] + 1 < numbersPrimeResidues[j]) {
				tupelCRT[j]++;
				currentPrimeResidues[j] = primeResidues[indices[j] + tupelCRT[j]];
				break;
			}else {
				tupelCRT[j] = 0;
				currentPrimeResidues[j] = primeResidues[indices[j] + tupelCRT[j]];
			}
		}
	}
	qsort(*residues, numberResidues_ui, sizeof(mpz_t), cmpfnc);	//order residue classes by smallest non-negative representative
	
	free(tupelCRT);
	free(currentPrimeResidues);
}

int noListResidue (mpz_t *factors, unsigned int numberFactors, unsigned short degree, short *poly, mpz_t C, mpz_t *currentPositions, mpz_t *residue) {
	
	unsigned int i;
	mpz_t j; mpz_init(j);
	unsigned int k;
	mpz_t eval; mpz_init(eval);
	mpz_t temp; mpz_init(temp);
	/*mpz_t *tempPrimeResidues = (mpz_t *) calloc(numberFactors, sizeof(mpz_t));
	if (tempPrimeResidues == NULL) {
		printf("Error: one of the arrays in findPrimeResidues could not be allocated\n");
		exit(-1);
	}*/
	mpz_t *poly_gmp = (mpz_t *) malloc(degree * sizeof(mpz_t));
	if (poly_gmp == NULL) {
		printf("Error: poly_gmp could not be allocated\n");
		exit(-1);
	}
	for (i = 0; i < degree; i++) {
		mpz_init_set_si(poly_gmp[i], poly[i]);
	}
	unsigned short foundResidue = 1;
	
	for (i = 0; i < numberFactors; i++) {	//find solutions for a(x) = 0 modulo prime power factors of C
		for (mpz_set(j, currentPositions[i]); 1; mpz_add_ui(j, j, 1)) {
			if (mpz_cmp(j, factors[i]) == 0) {	//walk through all combinations of residues mod prime power factors of C
				mpz_set_ui(j, 0);
				if (i + 1 < numberFactors) {
					mpz_add_ui(currentPositions[i + 1], currentPositions[i + 1], 1);
				}else {
					foundResidue = 0;	//checked all combinations
				}
			}
			mpz_set_ui(eval, 1);
			for (k = 0; k < degree; k++) {	//compute a(j) % factors[i]
				mpz_sub(temp, j, poly_gmp[k]);
				mpz_mul(eval, eval, temp);
				mpz_mod(eval, eval, factors[i]);	//eval = (eval * (j - poly[k])) % factors[i]; use % to keep size reasonable
			}
			if (mpz_cmp_ui(eval, 0) == 0) {
				mpz_set(currentPositions[i], j);
				break;
			}
		}
	}
	if (foundResidue == 1) {
		CRT(numberFactors, currentPositions, factors, C, residue);
	}
	
	for (i = 0; i < degree; i++) {
		mpz_clear(poly_gmp[i]);
	}
	mpz_clears(j, eval, temp, NULL);
	free(poly_gmp);
	return(foundResidue);
}

//for C < size:
//check if the residue classes have smooth representatives in a fixed interval [start, start + size]
void checkAllResiduesFixedInterval_SmallC (unsigned short *smoothNumbers, unsigned long numberResidues_ui, mpz_t *residues, mpz_t start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults) {
	
	unsigned long C_ui;
	if (mpz_cmp_ui(C, size) < 0) {
		C_ui = mpz_get_ui(C); //C < size < unsigned long
	}else {
		printf("Error: C is larger than size");
		exit(-1);
	}
	unsigned long startModC = mpz_fdiv_ui(start, C_ui);
	unsigned long i;
	unsigned long j;
	mpz_t temp; mpz_init(temp);
	unsigned short k;
	unsigned short smooth;

	for (i = 0; i < numberResidues_ui; i++) {
		mpz_sub_ui(temp, residues[i], startModC);	//find smalles representative in interval
		while (mpz_cmp_ui(temp, maxRoot) < 0) {	//prevent pointer out of array
			mpz_add(temp, temp, C);
		}
		if (mpz_cmp_ui(temp, size) < 0) {
			j = mpz_get_ui(temp);	//temp < size < unsigned long
		} else {
			continue;
		}
		while (j < size) {
			if (smoothNumbers[j] == 1) {
				smooth = 1;
				for (k = 1; k < numberRoots; k++) {
					if (smoothNumbers[j - roots[k]] == 0) {
						smooth = 0;
						break;
					}
				}
				if (smooth == 1) {
					mpz_add_ui(smoothIntsModC[*numberSmoothIntsModC], start, j);
					//gmp_printf("smooth int mod C: %Zd (new)\n", smoothIntsModC[*numberSmoothIntsModC]);
					(*numberSmoothIntsModC)++;
					if (*numberSmoothIntsModC == maxNumberResults) {
						printf("Error: reached max number of %u results\n", maxNumberResults);
						exit(-1);
					}
				}
			}
			j += C_ui;
		}
	}
	
	mpz_clear(temp);
}

//for C >= size:
//relevant steps to cover interval of size "size"
unsigned int findRelevantSteps (unsigned long numberResidues_ui, mpz_t *residues, unsigned int size, mpz_t C) {

	unsigned int relevantSteps = 1;
	unsigned long i = 0;
	mpz_t temp; mpz_init(temp);

	while (i < numberResidues_ui && relevantSteps + i < numberResidues_ui) {	//handle wrap-around
		mpz_sub(temp, residues[i + relevantSteps], residues[i]);
		if (mpz_cmp_ui(temp, size) < 0) {
			relevantSteps++;
			if (relevantSteps == numberResidues_ui) {
				return relevantSteps;
			}
		}else{
			i++;
		}
	}
	while (i < numberResidues_ui) {
		mpz_add(temp, residues[(i + relevantSteps) % numberResidues_ui], C);
		mpz_sub(temp, temp, residues[i]);
		if (mpz_cmp_ui(temp, size) < 0) {
			relevantSteps++;
			if (relevantSteps == numberResidues_ui) {
				return relevantSteps;
			}
		}else{
			i++;
		}
	}
	
	mpz_clear(temp);

	return relevantSteps;
}

//check if the residue classes have smooth representatives in a fixed interval [start, start + size]
void checkAllResiduesFixedInterval (unsigned short *smoothNumbers, unsigned long numberResidues_ui, mpz_t *residues, mpz_t start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, mpz_t C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults) {

	mpz_t startModC; mpz_init(startModC); mpz_mod(startModC, start, C);
	unsigned long i;
	unsigned long j;
	unsigned short k;
	mpz_t temp; mpz_init(temp);
	unsigned short smooth;
	unsigned long searchL = 0;
	unsigned long searchR = numberResidues_ui;
	unsigned long searchPivot = numberResidues_ui / 2;	//floor(numberResidues_ui / 2)
	unsigned long lastResidue;

	while (searchR - searchL > 0) {	//find smallest representative of R in chunk
		if (mpz_cmp(startModC, residues[searchPivot]) > 0) {
			searchL = searchPivot + 1;
		}else{
			searchR = searchPivot;
		}
		searchPivot = (searchL + searchR) / 2;	//floor((searchL + searchR) / 2)
	}
	lastResidue = searchPivot + relevantSteps;
	if (lastResidue > numberResidues_ui) {	//handle wrap-around
		for (i = 0; i < (lastResidue % numberResidues_ui); i++) {	//check if representatives are smooth
			mpz_sub(temp, residues[i], startModC);
			mpz_add(temp, temp, C); //temp = residues[i] - startModC + C;
			if (mpz_cmp_ui(temp, size) < 0 && mpz_cmp_ui(temp, maxRoot) >= 0) {
				j = mpz_get_ui(temp);
				if (smoothNumbers[j] == 1) {
					smooth = 1;
					for (k = 1; k < numberRoots; k++) {
						if (smoothNumbers[j - roots[k]] == 0) {
							smooth = 0;
							break;
						}
					}
					if (smooth == 1) {
						mpz_add_ui(smoothIntsModC[*numberSmoothIntsModC], start, j);
						//gmp_printf("smooth int mod C: %Zd (new)\n", smoothIntsModC[*numberSmoothIntsModC]);
						(*numberSmoothIntsModC)++;
						if (*numberSmoothIntsModC == maxNumberResults) {
							printf("Error: reached max number of %u results\n", maxNumberResults);
							exit(-1);
						}
					}
				}
			}
		}
		lastResidue = numberResidues_ui;
	}
	for (i = searchPivot; i < lastResidue; i++) {		//check if representatives are smooth
		mpz_sub(temp, residues[i], startModC);	//j = residues[i] - startModC;
		if (mpz_cmp_ui(temp, size) < 0 && mpz_cmp_ui(temp, maxRoot) >= 0) {
			j = mpz_get_ui(temp);
			if (smoothNumbers[j] == 1) {
				smooth = 1;
				for (k = 1; k < numberRoots; k++) {
					if (smoothNumbers[j - roots[k]] == 0) {
						smooth = 0;
						break;
					}
				}
				if (smooth == 1) {
					mpz_add_ui(smoothIntsModC[*numberSmoothIntsModC], start, j);
					//gmp_printf("smooth int mod C: %Zd (new)\n", smoothIntsModC[*numberSmoothIntsModC]);
					(*numberSmoothIntsModC)++;
					if (*numberSmoothIntsModC == maxNumberResults) {
						printf("Error: reached max number of %u results\n", maxNumberResults);
						exit(-1);
					}
				}
			}
		}
	}
	
	mpz_clears(startModC, temp, NULL);
}

//for numberResidues > unsigned long
//find largest power of two such that the product of all elements in an array of that length is still smaller than then product of the smooth primes
void findBatchSize (mpz_t start, mpz_t productPrimes, unsigned int *arrayLength) {
	
	unsigned int currentPower, nextPower;
	mpz_t tempProd; mpz_init_set_ui(tempProd, 1); 
	unsigned int log2elements = (unsigned int) mpz_sizeinbase(start, 2) - 1;
	unsigned int log2productPrimes = (unsigned int) mpz_sizeinbase(productPrimes, 2) - 1;
	*arrayLength = log2productPrimes / log2elements;	//floor(log2productPrimes / log2elements) is a first estimate for array length 
	nextPower = (unsigned int) floor(log(*arrayLength) / log(2));	//find suitable power of two
	mpz_mul_2exp(tempProd, tempProd, nextPower);	//mpz_ui_pow_ui(tempProd, 2, nextPower);
	nextPower = mpz_get_ui(tempProd);
	mpz_pow_ui(tempProd, start, nextPower);
	while (mpz_cmp(tempProd, productPrimes) <= 0) {
		currentPower = nextPower;
		nextPower *= 2;
		mpz_mul(tempProd, tempProd, tempProd);
	}
	*arrayLength = currentPower;
	/*unsigned short assertFlag = 0;
	unsigned int i;
	for (i = 2; i < 1048577; i *= 2) {
		if (i == *arrayLength) {
			assertFlag=1;
			break;
		}
	}
	assert(assertFlag == 1);*/
	mpz_clear(tempProd);
}

void checkFixedResidueFrKlMoWi (mpz_t productPrimes, mpz_t residue, mpz_t start, unsigned int batchLength, mpz_t *batch, mpz_t *productTreeArray, mpz_t *modArray, mpz_t *smoothnessArray, unsigned short numberRoots, mpz_t *roots_gmp, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults) {
	
	//initialize and set
	unsigned int i;
	unsigned short j;
	mpz_t temp; mpz_init(temp);
	unsigned int numberSmoothElements;
	unsigned int segmentSize = numberRoots - 1;
	unsigned int numberSmoothSegments;
	unsigned int *smoothSegments;
	unsigned int maxSegmentsRound;
	unsigned int numberRounds;
	unsigned int lastRoundSegments;
	unsigned int batchCounter;
	unsigned int currentSegment;
	unsigned int shortBatchLength = 1024;
	
	//set up batch
	mpz_mod(temp, start, C);
	mpz_sub(temp, start, temp);
	mpz_add(temp, temp, residue);	//find representative
	if (mpz_cmp(temp, start) < 0) {	//ensure representative >= start
		mpz_add(temp, temp, C);
	}
	if (mpz_cmp_ui(temp, maxRoot) < 0) {	//ensure (temp - root) is non-negative
		mpz_add(temp, temp, C);
	}
	mpz_set(batch[0], temp);
	for (i = 1; i < batchLength; i++) {
		mpz_add(batch[i], batch[i - 1], C);
	}
	
	//test batch (first factor l = l - 0)
	smoothBatchFrKlMoWiProduct(productPrimes, batchLength, batch, productTreeArray, modArray, &numberSmoothElements, smoothnessArray);
	
	//test other factors (l - roots[i]) for smooth l)
	//*numberSmoothIntsModC = 0;	//disable for long test
	maxSegmentsRound = batchLength / segmentSize;	//floor(batchLength / segmentSize)
	smoothSegments = (unsigned int*) malloc(maxSegmentsRound * sizeof(int));
	numberRounds = (numberSmoothElements + maxSegmentsRound - 1) / maxSegmentsRound;	//ceil(numberSmoothElements / maxSegmentsRound)
	lastRoundSegments = numberSmoothElements % maxSegmentsRound;
	currentSegment = 0;
	while (numberRounds > 1) {
		batchCounter = 0;
		for (i = 0; i < maxSegmentsRound; i++) {
			for (j = 1; j < numberRoots; j++) {	//roots[0] = 0 was tested before
				mpz_sub(batch[batchCounter], smoothnessArray[currentSegment], roots_gmp[j]);
				batchCounter++;
			}
			currentSegment++;
		}
		while (batchCounter < batchLength) {	//fill batch with ones (need length a power of two)
			mpz_set_ui(batch[batchCounter], 1);
			batchCounter++;
		}
		smoothBatchFrKlMoWiProductSegments(productPrimes, batchLength, batch, maxSegmentsRound, segmentSize, productTreeArray, modArray, &numberSmoothSegments, smoothSegments);
		for (i = 0; i < numberSmoothSegments; i++) {
			mpz_set(smoothIntsModC[*numberSmoothIntsModC], smoothnessArray[smoothSegments[i]]);
			(*numberSmoothIntsModC)++;
		}
		numberRounds--;
	}
	
	//last round
	batchCounter = 0;
	while (currentSegment < numberSmoothElements) {
		for (j = 1; j < numberRoots; j++) {	//roots[0] = 0 was tested before
			mpz_sub(batch[batchCounter], smoothnessArray[currentSegment], roots_gmp[j]);
			batchCounter++;
		}
		currentSegment++;
	}
	while (shortBatchLength < batchCounter) {
		shortBatchLength *= 2;
	}
	while (batchCounter < shortBatchLength) {	//fill batch with ones (need length a power of two)
		mpz_set_ui(batch[batchCounter], 1);
		batchCounter++;
	}
	smoothBatchFrKlMoWiProductSegments(productPrimes, shortBatchLength, batch, lastRoundSegments, segmentSize, productTreeArray, modArray, &numberSmoothSegments, smoothSegments);
	for (i = 0; i < numberSmoothSegments; i++) {
		mpz_set(smoothIntsModC[*numberSmoothIntsModC], smoothnessArray[smoothSegments[i]]);
		(*numberSmoothIntsModC)++;
	}
	
	//clear and free
	mpz_clear(temp);
	free(smoothSegments);
}

void checkFixedResidueBernstein (mpz_t productPrimes, mpz_t residue, mpz_t start, unsigned int batchLength, mpz_t *batch, mpz_t *productTreeArray, mpz_t *modArray, mpz_t *smoothnessArray, unsigned short numberRoots, mpz_t *roots_gmp, unsigned short maxRoot, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults) {
	
	//initialize and set
	unsigned int i;
	unsigned short j;
	mpz_t temp; mpz_init(temp);
	unsigned int numberSmoothElements;
	unsigned int segmentSize = numberRoots - 1;
	unsigned int numberSmoothSegments;
	unsigned int *smoothSegments;
	unsigned int maxSegmentsRound;
	unsigned int numberRounds;
	unsigned int lastRoundSegments;
	unsigned int batchCounter;
	unsigned int currentSegment;
	unsigned int shortBatchLength = 1024;
	mpz_t largestElement; mpz_init(largestElement);
	
	//set up batch
	mpz_mod(temp, start, C);
	mpz_sub(temp, start, temp);
	mpz_add(temp, temp, residue);	//find representative
	if (mpz_cmp(temp, start) < 0) {	//ensure representative >= start
		mpz_add(temp, temp, C);
	}
	if (mpz_cmp_ui(temp, maxRoot) < 0) {	//ensure (temp - root) is non-negative
		mpz_add(temp, temp, C);
	}
	mpz_set(batch[0], temp);
	for (i = 1; i < batchLength; i++) {
		mpz_add(batch[i], batch[i - 1], C);
	}
	
	//test batch (first factor l = l - 0)
	smoothBatchBernsteinProduct(productPrimes, batchLength, batch, productTreeArray, modArray, &numberSmoothElements, smoothnessArray);
	
	//test other factors (l - roots[i]) for smooth l)
	//*numberSmoothIntsModC = 0;	//disable for long test
	maxSegmentsRound = batchLength / segmentSize;	//floor(batchLength / segmentSize)
	smoothSegments = (unsigned int*) malloc(maxSegmentsRound * sizeof(int));
	numberRounds = (numberSmoothElements + maxSegmentsRound - 1) / maxSegmentsRound;	//ceil(numberSmoothElements / maxSegmentsRound)
	lastRoundSegments = numberSmoothElements % maxSegmentsRound;
	currentSegment = 0;
	while (numberRounds > 1) {
		batchCounter = 0;
		for (i = 0; i < maxSegmentsRound; i++) {
			for (j = 1; j < numberRoots; j++) {	//roots[0] = 0 was tested before
				mpz_sub(batch[batchCounter], smoothnessArray[currentSegment], roots_gmp[j]);
				batchCounter++;
			}
			currentSegment++;
		}
		mpz_sub(largestElement, smoothnessArray[currentSegment - 1], roots_gmp[1]);
		while (batchCounter < batchLength) {	//fill batch with ones (need length a power of two)
			mpz_set_ui(batch[batchCounter], 1);
			batchCounter++;
		}
		smoothBatchBernsteinProductSegments(productPrimes, batchLength, batch, largestElement, maxSegmentsRound, segmentSize, productTreeArray, modArray, &numberSmoothSegments, smoothSegments);
		for (i = 0; i < numberSmoothSegments; i++) {
			mpz_set(smoothIntsModC[*numberSmoothIntsModC], smoothnessArray[smoothSegments[i]]);
			(*numberSmoothIntsModC)++;
		}
		numberRounds--;
	}
	
	//last round
	batchCounter = 0;
	while (currentSegment < numberSmoothElements) {
		for (j = 1; j < numberRoots; j++) {	//roots[0] = 0 was tested before
			mpz_sub(batch[batchCounter], smoothnessArray[currentSegment], roots_gmp[j]);
			batchCounter++;
		}
		currentSegment++;
	}
	mpz_sub(largestElement, smoothnessArray[currentSegment - 1], roots_gmp[1]);
	while (shortBatchLength < batchCounter) {
		shortBatchLength *= 2;
	}
	while (batchCounter < shortBatchLength) {	//fill batch with ones (need length a power of two)
		mpz_set_ui(batch[batchCounter], 1);
		batchCounter++;
	}
	smoothBatchBernsteinProductSegments(productPrimes, shortBatchLength, batch, largestElement, lastRoundSegments, segmentSize, productTreeArray, modArray, &numberSmoothSegments, smoothSegments);
	for (i = 0; i < numberSmoothSegments; i++) {
		mpz_set(smoothIntsModC[*numberSmoothIntsModC], smoothnessArray[smoothSegments[i]]);
		(*numberSmoothIntsModC)++;
	}
	
	//clear and free
	mpz_clear(temp);
	free(smoothSegments);
}
