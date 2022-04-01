#ifndef NEW_H
#define NEW_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//find residue classes r that solve poly(r)=0 mod C
void findResidues (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short degree, short *poly, unsigned int C, unsigned int *numberResidues, unsigned int **residues);

//check relevant elements for C<size
void checkMod (unsigned short *smoothNumbers, unsigned int numberResidues, unsigned int *residues, unsigned long start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned int C, unsigned short *numberSmoothIntsModC, unsigned long *smoothIntsModC, unsigned short maxNumberResults);

//for C>=size:
//find relevant steps to cover the interval
unsigned int findRelevantSteps (unsigned int numberResidues, unsigned int *residues, unsigned int size, unsigned int C);

//check relevant elements for C>=size
void checkMOD (unsigned short *smoothNumbers, unsigned int numberResidues, unsigned int *residues, unsigned long start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned int C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, unsigned long *smoothIntsModC, unsigned short maxNumberResults);

#endif
