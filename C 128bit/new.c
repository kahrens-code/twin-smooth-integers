#include"new.h"

//compare function for sorting
int cmpfnc (const void *a, const void *b) {
	return ((*(unsigned long *)a > *(unsigned long *)b)-(*(unsigned long *)a < *(unsigned long *)b));
}

//Chinese remainder theorem for tupels
unsigned long CRTTupel (unsigned int number, unsigned long *tupel, unsigned long *indices, unsigned long *values, unsigned long *n, unsigned long NN) {
	
	__int128 N[number];
	__int128 M[number];
	unsigned int i;
	__int128 r, oldR;
	__int128 s, oldS;
	__int128 t, oldT;
	__int128 temp;
	__int128 quot;
	__int128 y=0;
	unsigned long x;

	for (i=0; i<number; i++) {
		N[i] = NN/n[i];
		oldR = N[i];
		r = n[i];
		oldS = 1;
		s = 0;
		oldT = 0;
		t = 1;
		while (r != 0) {
			quot = oldR / r;
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

void findResidues (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short degree, short *poly, unsigned long C, unsigned long *numberResidues, unsigned long **residues) {

	unsigned long i = 0;
	unsigned long j = C;

//factorise C	
	unsigned long factors[(unsigned long)ceil(log(C))];		//use calloc if C becomes too big
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
	unsigned long *primeResidues = (unsigned long *) calloc(C, sizeof(long));
	//unsigned long *primeResidues = (unsigned long *) calloc(102000000, sizeof(long)); //use this for the provided example of degree 7, if the other allocation fails
	unsigned long *numbersPrimeResidues = (unsigned long *) calloc(numberFactors, sizeof(long));
	unsigned long *indices = (unsigned long *) calloc(numberFactors + 1, sizeof(long));
	unsigned long *tupelCRT = (unsigned long *) calloc(numberFactors, sizeof(long));
	if (primeResidues == NULL || numbersPrimeResidues == NULL || indices == NULL || tupelCRT == NULL) {
		printf("Error: One of the arrays in findResidues could not be allocated\n");
		exit(-1);
	}
	__int128 eval;
	unsigned long k;
	unsigned long r = 0;

	*numberResidues = 1;
	for (i=0; i<numberFactors; i++) {		//find solutions for a(x) = 0 modulo prime power factors of C
		for (j=0; j<factors[i]; j++) {
			eval = 1;
			for (k=0; k<degree; k++) {
				eval = (eval * ((__int128)j - poly[k])) % (__int128)factors[i];		//use % to prevent overflow
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
	*residues = (unsigned long *) malloc(*numberResidues * sizeof(long));
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
	qsort(*residues, *numberResidues, sizeof(long), cmpfnc);	//order residue classes by smallest non-negative representative
	
	free(indices);
	free(tupelCRT);
	free(numbersPrimeResidues);
	free(primeResidues);
}

//for C < size:
//check if the residue classes have smooth representatives
void checkMod (unsigned short *smoothNumbers, unsigned long numberResidues, unsigned long *residues, unsigned __int128 start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned long C, unsigned short *numberSmoothIntsModC, unsigned __int128 *smoothIntsModC, unsigned short maxNumberResults) {
	
	unsigned int startModC = start % C;
	unsigned int i;
	long j;
	unsigned short k;
	unsigned short smooth;

	for (i=0; i<numberResidues; i++) {
		j = (long)residues[i] - startModC;	//find smalles representative in interval
		while (j < maxRoot) {				//prevent pointer out of array
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
					/*printf("smooth int mod C: ");
					printu128(start + i);
					printf(" (new)\n");*/
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
unsigned int findRelevantSteps (unsigned long numberResidues, unsigned long *residues, unsigned int size, unsigned long C) {

	unsigned int relevantSteps = 1;
	unsigned long i = 0;

	while ( i < numberResidues && relevantSteps + i < numberResidues) {	//handle wrap-around
		if (residues[i + relevantSteps] - residues[i] < size) {
			relevantSteps++;
			if (relevantSteps == numberResidues) {
				return numberResidues;
			}
		}else{
			i++;
		}
	}
	while (i < numberResidues) {
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
void checkMOD (unsigned short *smoothNumbers, unsigned long numberResidues, unsigned long *residues, unsigned __int128 start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned long C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, unsigned __int128 *smoothIntsModC, unsigned short maxNumberResults) {

	unsigned long startModC = start % C;
	unsigned long i;
	unsigned int j;
	unsigned short k;
	unsigned short smooth;
	unsigned long searchL = 0;
	unsigned long searchR = numberResidues;
	unsigned long searchPivot = numberResidues / 2;
	unsigned long lastResidue;

	while (searchR-searchL > 0) {		//find smallest representative of R in chunk
		if (startModC > residues[searchPivot]) {
			searchL = searchPivot+1;
		}else{
			searchR = searchPivot;
		}
		searchPivot = (searchL + searchR) / 2;
	}
	lastResidue = searchPivot + relevantSteps;
	if (lastResidue > numberResidues) {		//handle wrap-around
		for (i=0; i<(lastResidue%numberResidues); i++) {		//check if representatives are smooth
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
						/*printf("smooth int mod C: ");
						printu128(start + i);
						printf(" (new)\n");*/
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
	for (i=searchPivot; i<lastResidue; i++) {		//check if representatives are smooth
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
					/*printf("smooth int mod C: ");
					printu128(start + i);
					printf(" (new)\n");*/
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
