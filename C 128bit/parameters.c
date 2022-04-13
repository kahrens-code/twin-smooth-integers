#include"parameters.h"

//declare available examples
PTEsolution n6ex1, n6ex2, n6ex3, n6ex4, n6ex5;	//examples with degree n=6
PTEsolution n7ex1;								//example with degree n=7
PTEsolution n8ex1;								//example with degree n=8

//set all parameters accordingly to the chosen example PTE solution
//option for startParameter:
//0: 						start stays unchanged
//1, ..., PTE.numerPrimes:	set start to corresponding value
void setParameters (PTEsolution *PTE, unsigned short *degree, unsigned short *numberRoots, short **roots, short **polyA, short **polyB, unsigned long *C, unsigned __int128 *start, unsigned short startParameter, unsigned __int128 totalSize, unsigned __int128 *end) {

//define available examples
	n6ex1.degree = 6;
	n6ex1.numberRoots = 9;
	static short roots61[9] = {0, 1, 3, 5, 8, 11, 13, 15, 16};
	n6ex1.roots = roots61;
	static short poly61A[6] = {0, 3, 5, 11, 13, 16};
	n6ex1.polyA = poly61A;
	static short poly61B[6] = {1, 1, 8, 8, 15, 15};
	n6ex1.polyB = poly61B;
	n6ex1.C = 14400;
	n6ex1.numberPrimes = 4;
	static unsigned __int128 startPrimes61[4] = {5170314000000, 6781477000000, 9244655000000, 19052682000000};
	n6ex1.startPrimes = startPrimes61;

	n6ex2.degree = 6;
	n6ex2.numberRoots = 9;
	static short roots62[9] = {0, 2, 7, 8, 15, 22, 23, 28, 30}; 
	n6ex2.roots = roots62;
	static short poly62A[6] = {0, 7, 8, 22, 23, 30};
	n6ex2.polyA = poly62A;
	static short poly62B[6] = {2, 2, 15, 15, 28, 28};
	n6ex2.polyB = poly62B;
	n6ex2.C = 705600;
	n6ex2.numberPrimes = 2;
	static unsigned __int128 startPrimes62[2] = {22687888000000, 26042586000000};
	n6ex2.startPrimes = startPrimes62;

	n6ex3.degree = 6;
	n6ex3.numberRoots = 9;
	static short roots63[9] = {0, 2, 5, 16, 21, 26, 37, 40, 42};
	n6ex3.roots = roots63;
	static short poly63A[6] = {0, 5, 16, 26, 37, 42};
	n6ex3.polyA = poly63A;
	static short poly63B[6] = {2, 2, 21, 21, 40, 40};
	n6ex3.polyB = poly63B;
	n6ex3.C = 2822400;
	n6ex3.numberPrimes = 2;
	static unsigned __int128 startPrimes63[2] = {18675743000000, 36144284000000};
	n6ex3.startPrimes = startPrimes63;

	n6ex4.degree = 6;
	n6ex4.numberRoots = 9;
	static short roots64[9] = {0, 3, 7, 33, 40, 47, 73, 77, 80};
	n6ex4.roots = roots64;
	static short poly64A[6] = {0, 7, 33, 47, 73, 80};
	n6ex4.polyA = poly64A;
	static short poly64B[6] = {3, 3, 40, 40, 77, 77};
	n6ex4.polyB = poly64B;
	n6ex4.C = 85377600;
	n6ex4.numberPrimes = 2;
	static unsigned __int128 startPrimes64[2] = {19798693000000, 32519458000000};
	n6ex4.startPrimes = startPrimes64;

	n6ex5.degree = 6;
	n6ex5.numberRoots = 9;
	static short roots65[9] = {0, 4, 11, 24, 35, 46, 59, 66, 70};
	n6ex5.roots = roots65;
	static short poly65A[6] = {0, 11, 24, 46, 59, 70};
	n6ex5.polyA = poly65A;
	static short poly65B[6] = {4, 4, 35, 35, 66, 66};
	n6ex5.polyB = poly65B;
	n6ex5.C = 85377600;
	n6ex5.numberPrimes = 1;
	static unsigned __int128 startPrimes65[1] = {24924102000000};
	n6ex5.startPrimes = startPrimes65;

	n7ex1.degree = 7;
	n7ex1.numberRoots = 14;
	static short roots71[14] = {0, 1, 11, 18, 19, 30, 39, 50 ,56, 68 ,70, 79, 81, 84};
	n7ex1.roots = roots71;
	static short poly71A[7] = {0, 18, 19, 50, 56, 79, 81};
	n7ex1.polyA = poly71A;
	static short poly71B[7] = {1, 11, 30, 39, 68, 70, 84};
	n7ex1.polyB = poly71B;
	n7ex1.C = 5145940800;
	n7ex1.numberPrimes = 0;
	//static unsigned long startPrimes71[] = {};
	//n7ex1.startPrimes = startPrimes71;


	n8ex1.degree = 8;
	n8ex1.numberRoots = 16;
	static short roots81[16] = {0, 1, 2, 4, 9, 11, 20, 23, 27, 30, 39, 41, 46, 48, 49, 50};
	n8ex1.roots = roots81;
	static short poly81A[8] = {0, 4, 9, 23, 27, 41, 46, 50};
	n8ex1.polyA = poly81A;
	static short poly81B[8] = {1, 2, 11, 20 ,30, 39, 48, 49};
	n8ex1.polyB = poly81B;
	n8ex1.C = 1210809600;
	n8ex1.numberPrimes = 0;
	//static unsigned long startPrimes81[] = {};
	//n8ex1.startPrimes = startPrimes81;

//set parameters to the ones of the chosen example	
	*degree = PTE->degree;
	*numberRoots = PTE->numberRoots;
	*roots = PTE->roots;
	*polyA = PTE->polyA;
	*polyB = PTE->polyB;
	*C = PTE->C;
	if (startParameter > PTE->numberPrimes) {
		printf("Error: startParameter is larger than PTE->numberPrimes (%u). start is unchanged.\n",PTE->numberPrimes);
	}else{
		if (startParameter > 0) {
			*start = PTE->startPrimes[startParameter - 1];
			*end = *start + totalSize;
		}
	}
}
