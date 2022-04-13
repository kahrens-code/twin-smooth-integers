#include"u128.h"

void printu128 (unsigned __int128 u128) {
	
	unsigned __int128 e19 = 10000000000000000000u; //10^19
	unsigned __int128 e38 = e19*e19;
	if (u128 <= 18446744073709551615u) {
		printf("%lu", (unsigned long)u128);
	}else{
		if (u128 < e38) {
			printf("%lu%019lu", (unsigned long)(u128 / e19), (unsigned long)(u128 % e19));
		}else{
			printf("%lu%019lu%019lu", (unsigned long)(u128 / e38), (unsigned long)((u128 % e38) / e19), (unsigned long)(u128 % e19));
		}
	}
}

//integer powers of 2
//only used during pre-computation
unsigned __int128 pow2(unsigned short exponent) {
	
	unsigned short i;
	unsigned __int128 result = 1;
	
	if (exponent > 127) {
		printf("Error: 2^%u is bigger than 2^128\n", exponent);
		exit(-1);
	}
	for (i=0; i<exponent; i++) {
		result *= 2;
	}
	
	return result;
}
