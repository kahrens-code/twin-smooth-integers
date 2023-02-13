#include"parameters.h"

//declare available examples
PTEsolution n6ex1, n6ex2, n6ex3, n6ex4, n6ex5;	//examples with degree n=6
PTEsolution n7ex1, n7ex2, n7ex3, n7ex4, n7ex5, n7ex6, n7ex7, n7ex8, n7ex9, n7ex10;	//example with degree n=7
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
	static short roots71[14] = {0, 1, 11, 18, 19, 30, 39, 50, 56, 68, 70, 79, 81, 84};
	n7ex1.roots = roots71;
	static short poly71A[7] = {0, 18, 19, 50, 56, 79, 81};
	n7ex1.polyA = poly71A;
	static short poly71B[7] = {1, 11, 30, 39, 68, 70, 84};
	n7ex1.polyB = poly71B;
	n7ex1.C = 5145940800; //33 bit
	n7ex1.numberPrimes = 0;
	//static unsigned long startPrimes71[] = {};
	//n7ex1.startPrimes = startPrimes71;
	
	n7ex2.degree = 7;
	n7ex2.numberRoots = 14;
	static short roots72[14] = {0, 1, 13, 18, 27, 38, 44, 58, 64, 75, 84, 89, 101, 102};
	n7ex2.roots = roots72;
	static short poly72A[7] = {0, 18, 27, 58, 64, 89, 101};
	n7ex2.polyA = poly72A;
	static short poly72B[7] = {1, 13, 38, 44, 75, 84, 102};
	n7ex2.polyB = poly72B;
	n7ex2.C = 13967553600; //34 bit
	n7ex2.numberPrimes = 0;
	//static unsigned long startPrimes72[] = {};
	//n7ex2.startPrimes = startPrimes72;
	
	n7ex3.degree = 7;
	n7ex3.numberRoots = 14;
	static short roots73[14] = {0, 4, 11, 24, 31, 52, 57, 74, 106, 119, 126, 137, 147, 150};
	n7ex3.roots = roots73;
	static short poly73A[7] = {0, 24, 31, 74, 106, 137, 147};
	n7ex3.polyA = poly73A;
	static short poly73B[7] = {4, 11, 52, 57, 119, 126, 150};
	n7ex3.polyB = poly73B;
	n7ex3.C = 293318625600; //39 bit
	n7ex3.numberPrimes = 0;
	//static unsigned long startPrimes73[] = {};
	//n7ex3.startPrimes = startPrimes73;
	
	n7ex4.degree = 7;	//from http://eslpower.org/k123456.htm
	n7ex4.numberRoots = 14;
	static short roots74[14] = {0, 3, 9, 14, 43, 46, 133, 141, 156, 175, 176, 193, 199, 204};
	n7ex4.roots = roots74;
	static short poly74A[7] = {0, 14, 43, 141, 156, 193, 199};
	n7ex4.polyA = poly74A;
	static short poly74B[7] = {3, 9, 46, 133, 175, 176, 204};
	n7ex4.polyB = poly74B;
	n7ex4.C = 1037896675200; //40 bit
	n7ex4.numberPrimes = 0;
	//static unsigned long startPrimes74[] = {};
	//n7ex4.startPrimes = startPrimes74;
	
	n7ex5.degree = 7;
	n7ex5.numberRoots = 14;
	static short roots75[14] = {0, 1, 47, 59, 68, 87, 126, 142, 181, 200, 209, 221, 267, 268};
	n7ex5.roots = roots75;
	static short poly75A[7] = {0, 59, 68, 142, 181, 221, 267};
	n7ex5.polyA = poly75A;
	static short poly75B[7] = {1, 47, 87, 126, 200, 209, 268};
	n7ex5.polyB = poly75B;
	n7ex5.C = 5771633313600; //43 bit
	n7ex5.numberPrimes = 0;
	//static unsigned long startPrimes75[] = {};
	//n7ex5.startPrimes = startPrimes75;
	
	n7ex6.degree = 7;
	n7ex6.numberRoots = 14;
	static short roots76[14] = {0, 23, 39, 108, 118, 187, 316, 362, 491, 560, 570, 639, 655, 678};
	n7ex6.roots = roots76;
	static short poly76A[7] = {0, 108, 114, 300, 336, 474, 486};
	n7ex6.polyA = poly76A;
	static short poly76B[7] = {6, 66, 180, 234, 408, 420, 504};
	n7ex6.polyB = poly76B;
	n7ex6.C = 1440534083788800; //51 bit
	n7ex6.numberPrimes = 0;
	//static unsigned long startPrimes76[] = {};
	//n7ex6.startPrimes = startPrimes76;
	
	n7ex7.degree = 7;
	n7ex7.numberRoots = 14;
	static short roots77[14] = {0, 23, 39, 108, 118, 187, 316, 362, 491, 560, 570, 639, 655, 678};
	n7ex7.roots = roots77;
	static short poly77A[7] = {0, 108, 118, 362, 491, 639, 655};
	n7ex7.polyA = poly77A;
	static short poly77B[7] = {23, 39, 187, 316, 560, 570, 678};
	n7ex7.polyB = poly77B;
	n7ex7.C = 11471328290822400; //54 bit
	n7ex7.numberPrimes = 0;
	//static unsigned long startPrimes77[] = {};
	//n7ex7.startPrimes = startPrimes77;
	
	n7ex8.degree = 7;	//parametric solution from ??? for -3
	n7ex8.numberRoots = 14;
	static short roots78[14] = {0, 19, 53, 116, 148, 245, 291, 397, 443, 540, 572, 635, 669, 688};
	n7ex8.roots = roots78;
	static short poly78A[7] = {0, 116, 148, 397, 443, 635, 669};
	n7ex8.polyA = poly78A;
	static short poly78B[7] = {19, 53, 245, 291, 540, 572, 688};
	n7ex8.polyB = poly78B;
	n7ex8.C = 15256916548473600; //54 bit
	n7ex8.numberPrimes = 0;
	//static unsigned long startPrimes78[] = {};
	//n7ex8.startPrimes = startPrimes78;
	
	n7ex9.degree = 7;
	n7ex9.numberRoots = 14;
	static short roots79[14] = {0, 28, 53, 111, 240, 323, 415, 521, 613, 696, 825, 883, 908, 936};
	n7ex9.roots = roots79;
	static short poly79A[7] = {0, 111, 240, 521, 613, 883, 908};
	n7ex9.polyA = poly79A;
	static short poly79B[7] = {28, 53, 323, 415, 696, 825, 936};
	n7ex9.polyB = poly79B;
	n7ex9.C = 106911286818336000; //57 bit
	n7ex9.numberPrimes = 0;
	//static unsigned long startPrimes79[] = {};
	//n7ex9.startPrimes = startPrimes79;
	
	n7ex10.degree = 7;
	n7ex10.numberRoots = 14;
	static short roots710[14] = {0, 17, 123, 222, 235, 422, 440, 650, 668, 855, 868, 967, 1073, 1090};
	n7ex10.roots = roots710;
	static short poly710A[7] = {0, 222, 235, 650, 668, 967, 1073};
	n7ex10.polyA = poly710A;
	static short poly710B[7] = {17, 123, 422, 440, 855, 868, 1090};
	n7ex10.polyB = poly710B;
	n7ex10.C = 314073647406288000; //59 bit
	n7ex10.numberPrimes = 0;
	//static unsigned long startPrimes710[] = {};
	//n7ex10.startPrimes = startPrimes710;

	n8ex1.degree = 8;
	n8ex1.numberRoots = 16;
	static short roots81[16] = {0, 1, 2, 4, 9, 11, 20, 23, 27, 30, 39, 41, 46, 48, 49, 50};
	n8ex1.roots = roots81;
	static short poly81A[8] = {0, 4, 9, 23, 27, 41, 46, 50};
	n8ex1.polyA = poly81A;
	static short poly81B[8] = {1, 2, 11, 20, 30, 39, 48, 49};
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
