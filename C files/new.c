#include"new.h"

//compare function for sorting
int cmpfnc (const void *a, const void *b) {
	return ((*(unsigned int *)a > *(unsigned int *)b)-(*(unsigned int *)a < *(unsigned int *)b));
}

//Chinese remainder theorem for tupels
unsigned int CRTTupel (unsigned int number, unsigned short *tupel, unsigned short *indices, unsigned int *values, unsigned int *n, unsigned int NN) {
	long N[number];
	long M[number];
	unsigned int i;
	long r, oldR;
	long s, oldS;
	long t, oldT;
	long temp;
	long quot;
	long y=0;
	unsigned int x;

	for (i=0; i<number; i++) {
		N[i] = NN/n[i];
		oldR = N[i];
		r = n[i];
		oldS = 1;
		s = 0;
		oldT = 0;
		t = 1;
		while (r != 0) {
			quot = floor((double) oldR / r);
			temp = oldR;
			oldR = r;
			r = temp - quot * r;
			temp = oldS;
			oldS = s;
			s = temp - quot * s;
			temp = oldT;
			oldT = t;
			t = temp - quot * t;
		}
		M[i] = oldS;
	}
	for (i=0; i<number; i++) {
		y = (y + values[indices[i] + tupel[i]] * M[i] * N[i]) % NN;
	}

	while (y < 0) {
		y += NN;
	}
	x = y;

	return x;
}

void findResidues (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short degree, short *poly, unsigned int C, unsigned int *numberResidues, unsigned int **residues) {

	unsigned int i = 0;
	unsigned int j = C;

//factorise C	
	unsigned int factors[(unsigned int)ceil(log(C))];	//use calloc if C becomes too big
	unsigned int numberFactors = 0;
	unsigned short newFactor = 0;

	factors[0] = 1;
	while (j > 1) {
		if (j % smoothPrimes[i] == 0) {
			j /= smoothPrimes[i];
			factors[numberFactors] *= smoothPrimes[i];
			newFactor = 1;
		}else {
			i++;
			if (i == numberSmoothPrimes) {
				printf("Error: C not smooth\n");
				exit(-1);
			}
			if (newFactor == 1) {
				numberFactors++;
				if (numberFactors == (unsigned int)ceil(log(C))) {
					printf("Error: C has more than log(C) prime factors");
					exit(-1);
				}
				factors[numberFactors] = 1;
				newFactor = 0;
			}
		}
	}
	numberFactors++;
	
//find residue classes such that poly(res) = 0 mod C	
	unsigned int *primeResidues = (unsigned int *) calloc(C, sizeof(int));
	unsigned int *numbersPrimeResidues = (unsigned int *) calloc(numberFactors, sizeof(int));
	unsigned short *indices = (unsigned short *) calloc(numberFactors + 1, sizeof(short));
	unsigned short *tupelCRT = (unsigned short *) calloc(numberFactors, sizeof(short));
	if (primeResidues == NULL || numbersPrimeResidues == NULL || indices == NULL || tupelCRT == NULL) {
		printf("Error: One of the arrays could not be allocated\n");
		exit(-1);
	}
	long eval;
	unsigned int k;
	unsigned int r = 0;

	*numberResidues = 1;
	for (i=0; i<numberFactors; i++) {	//find solutions for a(x) = 0 modulo prime power factors of C
		for (j=0; j<factors[i]; j++) {
			eval = 1;
			for (k=0; k<degree; k++) {
				eval = (eval * ((long)j - poly[k])) % (long)factors[i];	//use % to prevent overflow
			}
			if (eval == 0) {
				primeResidues[r] = j;
				r++;
				numbersPrimeResidues[i]++;
			}
		}
		*numberResidues *= numbersPrimeResidues[i];
		indices[i + 1] = indices[i] + numbersPrimeResidues[i];
	}
	*residues = (unsigned int *) malloc(*numberResidues * sizeof(int));
	if (*residues == NULL) {
                printf("Error: residues could not be allocated\n");
                exit(-1);
        }
	for (i=0; i<*numberResidues; i++) {		//use CRT to combine solutions mod prime powers to solutions mod C
		(*residues)[i] = CRTTupel(numberFactors, tupelCRT, indices, primeResidues, factors, C);
		for (j=0; j<numberFactors; j++) {
			if (tupelCRT[j] + 1 < numbersPrimeResidues[j]) {
				tupelCRT[j]++;
				for (k=0; k<j; k++) {
					tupelCRT[k] = 0;
				}
				break;
			}
		}
	}
	qsort(*residues, *numberResidues, sizeof(int), cmpfnc);	//order residue classes by smallest non-negative representative
	
	free(indices);
	free(tupelCRT);
	free(numbersPrimeResidues);
	free(primeResidues);
}

//for C < size:
//check if the residue classes have smooth representatives
void checkMod (unsigned short *smoothNumbers, unsigned int numberResidues, unsigned int *residues, unsigned long start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned int C, unsigned short *numberSmoothIntsModC, unsigned long *smoothIntsModC, unsigned short maxNumberResults) {
	
	unsigned int startModC = start % C;
	unsigned int i;
	long j;
	unsigned int k;
	unsigned short smooth;

	for (i=0; i<numberResidues; i++) {
		j = (long)residues[i] - startModC;	//find smalles representative in interval
		while (j < maxRoot) {			//prevent pointer out of array
			j += C;
		}
		while (j < size) {
			if (smoothNumbers[j] == 1) {
				smooth = 1;
				for (k=1; k<numberRoots; k++) {
					if (smoothNumbers[j - roots[k]] == 0) {
						smooth = 0;
						break;
					}
				}
				if (smooth == 1) {
					//printf("smooth int mod C: %lu (new)\n", start + j);
					smoothIntsModC[*numberSmoothIntsModC] = start+j;
					(*numberSmoothIntsModC)++;
					if (*numberSmoothIntsModC == maxNumberResults) {
						printf("Error: reached max number of %u results\n", maxNumberResults);
						exit(-1);
					}
				}
			}
			j += C;
		}
	}
}

//for C >= size:
//relevant steps to cover chunk
unsigned int findRelevantSteps (unsigned int numberResidues, unsigned int *residues, unsigned int size, unsigned int C) {

	unsigned int relevantSteps = 1;
	unsigned int i = 0;

	while ( i <numberResidues && relevantSteps + i < numberResidues) {	//handle wrap-around
		if (residues[i + relevantSteps] - residues[i] < size) {
			relevantSteps++;
			if (relevantSteps == numberResidues) {
				return numberResidues;
			}
		}else{
			i++;
		}
	}
	while (i<numberResidues) {
		if (residues[(i + relevantSteps) % numberResidues] + C - residues[i] < size) {
			relevantSteps++;
			if (relevantSteps == numberResidues) {
				return numberResidues;
			}
		}else{
			i++;
		}
	}

	return relevantSteps;
}

//check if the residue classes have smooth representatives
void checkMOD (unsigned short *smoothNumbers, unsigned int numberResidues, unsigned int *residues, unsigned long start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned int C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, unsigned long *smoothIntsModC, unsigned short maxNumberResults) {

	unsigned int startModC = start % C;
	unsigned int i;
	long j;
	unsigned int k;
	unsigned short smooth;
	unsigned int searchL = 0;
	unsigned int searchR = numberResidues;
	unsigned int searchPivot = floor((double)numberResidues / 2);
	unsigned int lastResidue;

	while (searchR-searchL > 0) {		//find smallest representative of R in chunk
		if (startModC > residues[searchPivot]) {
			searchL = searchPivot+1;
		}else{
			searchR = searchPivot;
		}
		searchPivot = floor((double)(searchL + searchR) / 2);
	}
	lastResidue = searchPivot + relevantSteps;
	if (lastResidue > numberResidues) {				//handle wrap-around
		for (i=0; i<(lastResidue%numberResidues); i++) {	//check if representatives are smooth
			j = residues[i] - startModC + C;
			if (j < size && maxRoot <= j) {
				if (smoothNumbers[j] == 1) {
					smooth = 1;
					for (k=1; k<numberRoots; k++) {
						if (smoothNumbers[j - roots[k]] == 0) {
							smooth = 0;
							break;
						}
					}
					if (smooth == 1) {
						//printf("smooth int mod C: %lu (new)\n", start + j);
						smoothIntsModC[*numberSmoothIntsModC] = start + j;
						(*numberSmoothIntsModC)++;
						if (*numberSmoothIntsModC == maxNumberResults) {
							printf("Error: reached max number of %u results\n", maxNumberResults);
							exit(-1);
						}
					}
				}
			}
		}
	}
	if (numberResidues < lastResidue) {
		lastResidue = numberResidues;
	}
	for (i=searchPivot; i<lastResidue; i++) {	//check if representatives are smooth
		j = residues[i] - startModC;
		if (j < size && maxRoot <= j) {
			if (smoothNumbers[j] == 1) {
				smooth = 1;
				for (k=1; k<numberRoots; k++) {
					if (smoothNumbers[j - roots[k]] == 0) {
						smooth = 0;
						break;
					}
				}
				if (smooth == 1) {
					//printf("smooth int mod C: %lu (new)\n", start+j);
					smoothIntsModC[*numberSmoothIntsModC] = start + j;
					(*numberSmoothIntsModC)++;
					if (*numberSmoothIntsModC == maxNumberResults) {
						printf("Error: reached max number of %u results\n", maxNumberResults);
						exit(-1);
					}
				}
			}
		}
	}
}
