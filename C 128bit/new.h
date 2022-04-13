#ifndef NEW_H
#define NEW_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

//find residue classes r that solve poly(r)=0 mod C
void findResidues (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned short degree, short *poly, unsigned long C, unsigned long *numberResidues, unsigned long **residues);

//check relevant elements for C<size
void checkMod (unsigned short *smoothNumbers, unsigned long numberResidues, unsigned long *residues, unsigned __int128 start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned long C, unsigned short *numberSmoothIntsModC, unsigned __int128 *smoothIntsModC, unsigned short maxNumberResults);

//for C>=size:
//find relevant steps to cover the interval
unsigned int findRelevantSteps (unsigned long numberResidues, unsigned long *residues, unsigned int size, unsigned long C);

//check relevant elements for C>=size
void checkMOD (unsigned short *smoothNumbers, unsigned long numberResidues, unsigned long *residues, unsigned __int128 start, unsigned int size, unsigned short numberRoots, short *roots, unsigned short maxRoot, unsigned long C, unsigned int relevantSteps, unsigned short *numberSmoothIntsModC, unsigned __int128 *smoothIntsModC, unsigned short maxNumberResults);

#endif
