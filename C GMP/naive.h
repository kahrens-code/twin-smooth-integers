#ifndef NAIVE_H
#define NAIVE_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gmp.h>


void findSmoothIntsModC (mpz_t start, unsigned int size, unsigned short *smoothNumbers, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned short degree, short *poly, mpz_t C, unsigned short *numberSmoothIntsModC, mpz_t *smoothIntsModC, unsigned short maxNumberResults);

#endif
