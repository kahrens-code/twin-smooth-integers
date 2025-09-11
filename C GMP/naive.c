#include"naive.h"

//open: modulo with negatives

void findSmoothIntsModC (mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned short degree, short *poly, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults) {

	unsigned int i;
	unsigned short j;
	unsigned short smooth;
	mpz_t eval; mpz_init(eval);
	mpz_t tempEval; mpz_init(tempEval);
	mpz_t intModC; mpz_init(intModC);

	for (i = maxRoot; i < size; i++) {
		if (smoothNumbers[i] == 1) {							//check if l = start + i is smooth
			smooth = 1;
			for (j = 1; j < numberRoots; j++) {
				if (smoothNumbers[i - roots[j]] == 0) {			//then check other factors l - a_i
					smooth = 0;
					break;
				}
			}
			if (smooth == 1) {									//check if a(l) = 0 mod C
				mpz_set_ui(eval, 1);
				mpz_add_ui(intModC, start, i);
				mpz_mod(intModC, intModC, C);	//intModC = (start + i) % C;
				for (j = 0; j < degree; j++) {
					mpz_set_si(tempEval, poly[j]);
					mpz_sub(tempEval, intModC, tempEval);	//tempEval = intModC - poly[j]
					mpz_mul(eval, eval, tempEval);
					mpz_mod(eval, eval, C);	//eval = (eval * (intModC - poly[j])) % C; use mod C to keep size reasonable
				}
				if (mpz_cmp_ui(eval, 0) == 0) {								//write the found values into an array
					//gmp_printf("smooth int mod C: %Zd (naive)\n");
					mpz_add_ui(smoothIntsModC[*numberSmoothIntsModC], start, i);
					(*numberSmoothIntsModC)++;
					if (*numberSmoothIntsModC == maxNumberResults) {
						printf("Error: reached max number of %u results\n", maxNumberResults);
						exit(-1);
					}
				}
			}
		}
	}
	
	mpz_clears(eval, tempEval, intModC, NULL);
}
