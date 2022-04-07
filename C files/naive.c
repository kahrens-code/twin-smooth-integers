#include"naive.h"

void findSmoothIntsModC (unsigned long start, unsigned int size, unsigned short *smoothNumbers, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned short degree, short *poly, unsigned int C, unsigned short *numberSmoothIntsModC, unsigned long *smoothIntsModC, unsigned short maxNumberResults) {

	unsigned int i;
	unsigned int j;
	unsigned short smooth;
	long eval;
	unsigned int intModC;

	for (i=maxRoot; i<size; i++) {
		if (smoothNumbers[i] == 1) {		//check if l = start + i is smooth
			smooth = 1;
			for (j=1; j<numberRoots; j++) {
				if (smoothNumbers[i - roots[j]] == 0) {	//then check other factors l - a_i
					smooth = 0;
					break;
				}
			}
			if (smooth == 1) {		//check if a(l) = 0 mod C
				eval = 1;
				intModC = (start + i) % C;
				for (j=0; j<degree; j++) {
					eval = (eval * (intModC - poly[j])) % C;	//use mod C to prevent overflow
				}
				if (eval == 0) {		//write the found values into an array
					//printf("smooth int mod C: %lu (naive)\n", start + i);
					smoothIntsModC[*numberSmoothIntsModC] = start + i;
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
