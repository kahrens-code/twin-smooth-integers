#include"batchSmoothness.h"

//product of mpz_t array as binary tree
void treeProduct (unsigned int arrayLength, mpz_t *array, mpz_t *product) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays 
	mpz_t *tempArray1 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	mpz_t *tempArray2 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempArray1 == NULL || tempArray2 == NULL) {
		printf("Error: tempArrays could not be allocated\n");
		exit(-1);
	}
	mpz_t *PtrArrayA = array;
	mpz_t *PtrArrayB = tempArray2;
	mpz_t *PtrSwap = tempArray1;
	unsigned int length = arrayLength;
	unsigned int currentNew = 0;
	unsigned int currentOld = 0;
	unsigned int i;
	
	for (i = 0; i < length; i++){
		mpz_init(tempArray1[i]);
		mpz_init(tempArray2[i]);
	}
	
//compute product
	while (length > 1) {
		currentNew = 0;
		currentOld = 0;
		while (currentOld < length - 1) {
			mpz_mul(PtrArrayB[currentNew], PtrArrayA[currentOld], PtrArrayA[currentOld + 1]);	//multiply two entries of the old level for one entry of the new level
			currentNew++;
			currentOld += 2;
		}
		if (currentOld < length) {	//add the last element of the old level in case the number of entries is odd
			mpz_set(PtrArrayB[currentNew], PtrArrayA[currentOld]);
			currentNew++;
		}
		PtrArrayA = PtrArrayB;	//swap old and new arrays
		PtrArrayB = PtrSwap;
		PtrSwap = PtrArrayA;
		length = currentNew;
	}
	mpz_set(*product, PtrArrayA[0]);
	
//clear and free
	for (i = 0; i < length; i++){
		mpz_clear(tempArray1[i]);
		mpz_clear(tempArray2[i]);
	}
	free(tempArray1);
	free(tempArray2);
}

//product of unsigned int array as binary tree
void treeProduct_ui (unsigned int arrayLength, unsigned int *array, mpz_t *product) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays 
	mpz_t *tempArray1 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	mpz_t *tempArray2 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempArray1 == NULL || tempArray2 == NULL) {
		printf("Error: tempArrays could not be allocated\n");
		exit(-1);
	}
	mpz_t *PtrArrayA = tempArray1;
	mpz_t *PtrArrayB = tempArray2;
	mpz_t *PtrSwap = tempArray1;
	unsigned int length = arrayLength;
	unsigned int currentNew = 0;
	unsigned int currentOld = 0;
	unsigned int i;
	for (i = 0; i < length; i++){
		mpz_init(tempArray1[i]);
		mpz_init(tempArray2[i]);
	}
	
//compute product
	if (length > 1) {	//first loop also converts ui into mpz
		while (currentOld < length - 1) {
			mpz_set_ui(PtrArrayA[currentNew], array[currentOld]);
			currentOld++;
			mpz_mul_ui(PtrArrayA[currentNew], PtrArrayA[currentNew], array[currentOld]);	//multiply two entries of the old level for one entry of the new level
			currentNew++;
			currentOld++;
		}
		if (currentOld < length) {	//add the last element of the old level in case the number of entries is odd
			mpz_set_ui(PtrArrayA[currentNew], array[currentOld]);
			currentNew++;
		}
		length = currentNew;
	}
	else {
		mpz_set_ui(PtrArrayA[0], array[0]);
	}
	while (length > 1) {
		currentNew = 0;
		currentOld = 0;
		while (currentOld < length - 1) {
			mpz_mul(PtrArrayB[currentNew], PtrArrayA[currentOld], PtrArrayA[currentOld + 1]);	//multiply two entries of the old level for one entry of the new level
			currentNew++;
			currentOld += 2;
		}
		if (currentOld < length) {	//add the last element of the old level in case the number of entries is odd
			mpz_set(PtrArrayB[currentNew], PtrArrayA[currentOld]);
			currentNew++;
		}
		PtrArrayA = PtrArrayB;	//swap old and new arrays
		PtrArrayB = PtrSwap;
		PtrSwap = PtrArrayA;
		length = currentNew ;
	}
	mpz_set(*product, PtrArrayA[0]);
	
//clear and free
	for (i = 0; i < length; i++){
		mpz_clear(tempArray1[i]);
		mpz_clear(tempArray2[i]);
	}
	free(tempArray1);
	free(tempArray2);
}

//naive product of mpz_t array
void naiveProduct (unsigned int arrayLength, mpz_t *array, mpz_t *product) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	unsigned int i;
	mpz_set_ui(*product, 1);
	for (i = 0; i < arrayLength; i++) {
		mpz_mul(*product, *product, array[i]);
	}
}

//naive product of unsigned int array
void naiveProduct_ui (unsigned int arrayLength, unsigned int *array, mpz_t *product) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	unsigned int i;
	mpz_set_ui(*product, 1);
	for (i = 0; i < arrayLength; i++) {
		mpz_mul_ui(*product, *product, array[i]);
	}
}

//product of mpz_t array as binary tree
//save intermediate steps
//only for arrays of length a power of two
void treeProductSave (unsigned int arrayLength, mpz_t *array, mpz_t *tree) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays
	mpz_t *tempArray1 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	mpz_t *tempArray2 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempArray1 == NULL || tempArray2 == NULL) {
		printf("Error: tempArrays could not be allocated\n");
		exit(-1);
	}
	mpz_t *PtrArrayA = array;
	mpz_t *PtrArrayB = tempArray2;
	mpz_t *PtrSwap = tempArray1;
	unsigned int length = arrayLength;
	unsigned int currentNew = 0;
	unsigned int currentOld = 0;
	unsigned int currentTree = 0;
	unsigned int i;
	for (i = 0; i < arrayLength; i++){
		mpz_init(tempArray1[i]);
		mpz_init(tempArray2[i]);
	}
	
//compute product and save intermediate steps
	while (length > 1) {
		currentNew = 0;
		currentOld = 0;
		while (currentOld < length - 1) {
			mpz_mul(PtrArrayB[currentNew], PtrArrayA[currentOld], PtrArrayA[currentOld + 1]);	//multiply two entries of the old level for one entry of the new level
			mpz_set(tree[currentTree], PtrArrayA[currentOld]);	//save current level in tree
			currentOld++;
			currentTree++;
			mpz_set(tree[currentTree], PtrArrayA[currentOld]);
			currentOld++;
			currentTree++;
			currentNew++;
		}
		PtrArrayA = PtrArrayB;	//swap old and new arrays
		PtrArrayB = PtrSwap;
		PtrSwap = PtrArrayA;
		length = currentNew ;
	}
	mpz_set(tree[currentTree], PtrArrayA[0]);	//add final element (the procuct) to the tree
	
//clear and free
	for (i = 0; i < arrayLength; i++){
		mpz_clear(tempArray1[i]);
		mpz_clear(tempArray2[i]);
	}
	free(tempArray1);
	free(tempArray2);
}

//input modulo all elements in array
void naiveMod (mpz_t input, unsigned int arrayLength, mpz_t *array, mpz_t *modArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	unsigned int i;
	for (i = 0; i < arrayLength; i++) {
		mpz_mod(modArray[i], input, array[i]);
	}
}

//input modulo all leaves in tree as binary tree
//only for arrays of length a power of two
void treeMod (mpz_t input, unsigned int arrayLength, mpz_t *tree, mpz_t *modArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays
	mpz_t *tempArray1 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	mpz_t *tempArray2 = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempArray1 == NULL || tempArray2 == NULL) {
		printf("Error: tempArrays could not be allocated\n");
		exit(-1);
	}
	
	mpz_t *PtrArrayA = tempArray1;
	mpz_t *PtrArrayB = tempArray2;
	mpz_t *PtrSwap = tempArray1;
	unsigned int level = 0;
	unsigned int length = 1;
	unsigned int maxLevel = (unsigned int) (log(arrayLength)/log(2));
	unsigned int currentNew = arrayLength - 1;
	unsigned int currentOld = arrayLength - 1;
	unsigned int currentTree = 2 * arrayLength - 2;
	unsigned int i;
	for (i = 0; i < arrayLength; i++){
		mpz_init(tempArray1[i]);
		mpz_init(tempArray2[i]);
	}

//compute modulo
	mpz_mod(PtrArrayA[arrayLength - 1], input, tree[currentTree]);	//input mod product of array
	for (level = 1; level < maxLevel + 1; level++){
		currentNew = arrayLength - 1;
		currentOld = arrayLength - 1;
		for (i = 0; i < length; i++){
			currentTree--;	//we go "backwards" to return to the origignal order of the array
			mpz_mod(PtrArrayB[currentNew], PtrArrayA[currentOld], tree[currentTree]);	//compute the current level modulo two elements of the next level
			currentNew--;
			currentTree--;
			mpz_mod(PtrArrayB[currentNew], PtrArrayA[currentOld], tree[currentTree]);
			currentNew--;
			currentOld--;
		}
		PtrArrayA = PtrArrayB;	//swap old and new arrays
		PtrArrayB = PtrSwap;
		PtrSwap = PtrArrayA;
		length *= 2;
	}
	
	//*modArray = PtrArrayA;	//set output to the final level (leaves) of the tree

//clear and free
	for (i = 0; i < arrayLength; i++){
		mpz_set(modArray[i], PtrArrayA[i]);
		mpz_clear(PtrArrayA[i]);
		mpz_clear(PtrArrayB[i]);
	}
	free(PtrArrayA);
	free(PtrArrayB);
}

//test elements in array for smoothness given a list of all smooth primes FrKlMoWi 
void smoothBatchFrKlMoWiList (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned int arrayLength, mpz_t *array, unsigned int *numberSmoothElements, mpz_t **smoothnessArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays
	mpz_t productPrimes; mpz_init(productPrimes);
	mpz_t *productTreeArray = (mpz_t *) malloc((2 * arrayLength - 1) * sizeof(mpz_t));
	if (productTreeArray == NULL) {
		printf("Error: productTreeArray could not be allocated\n");
		exit(-1);
	}
	mpz_t *modArray = (mpz_t *) malloc((arrayLength) * sizeof(mpz_t));
	if (modArray == NULL) {
		printf("Error: modArray could not be allocated\n");
		exit(-1);
	}
	mpz_t *tempSmoothArray = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempSmoothArray == NULL) {
		printf("Error: tempSmoothArray could not be allocated\n");
		exit(-1);
	}
	unsigned int i;
	for (i = 0; i < arrayLength; i++){
		mpz_init(modArray[i]);
		mpz_init(tempSmoothArray[i]);
		mpz_init(productTreeArray[i]);
	}
	for (i = arrayLength; i < 2 * arrayLength - 1; i++){
		mpz_init(productTreeArray[i]);
	}
	
//compute products and modulo
	treeProduct_ui(numberSmoothPrimes, smoothPrimes, &productPrimes);	//compute product of smooth primes
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute FrKlMoWi smoothness algorithm
	*numberSmoothElements = 0;
	for (i = 0; i < arrayLength; i++) {
		mpz_set(tempSmoothArray[i], array[i]);
		while (mpz_cmp_ui(modArray[i], 1) > 0) {	//Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
			//mpz_mul(modArray[i], modArray[i], modArray[i]); //possible speedup
			mpz_gcd(modArray[i], modArray[i], tempSmoothArray[i]); 
			mpz_divexact(tempSmoothArray[i], tempSmoothArray[i], modArray[i]);
		}
		if (mpz_cmp_ui(tempSmoothArray[i], 1) == 0 || mpz_cmp_ui(modArray[i], 0) == 0) {
			mpz_set(tempSmoothArray[*numberSmoothElements], array[i]);
			(*numberSmoothElements)++;
		}
	}
	
//write smooth elements in shorter array
	if (*numberSmoothElements > 0) {
		*smoothnessArray = (mpz_t *) malloc(*numberSmoothElements * sizeof(mpz_t));
		if (smoothnessArray == NULL) {
			printf("Error: smoothnessArray could not be allocated\n");
			exit(-1);
		}
		for (i = 0; i < *numberSmoothElements; i++) {
			mpz_init_set((*smoothnessArray)[i], tempSmoothArray[i]); //remember to clear them in main function!
		}
	}
	
//clear and free
	mpz_clear(productPrimes);
	for (i = 0; i < arrayLength; i++){
		mpz_clear(modArray[i]);
		mpz_clear(tempSmoothArray[i]);
		mpz_clear(productTreeArray[i]);
	}
	for (i = arrayLength; i < 2 * arrayLength - 1; i++){
		mpz_clear(productTreeArray[i]);
	}
	free(productTreeArray);
	free(modArray);
	free(tempSmoothArray);
}

//test elements in array for smoothness given the product of all smooth primes FrKlMoWi
void smoothBatchFrKlMoWiProduct (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothElements, mpz_t *smoothnessArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	unsigned int i;
	
//compute products and modulo
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute FrKlMoWi smoothness algorithm
	*numberSmoothElements = 0;
	for (i = 0; i < arrayLength; i++) {
		mpz_set(smoothnessArray[i], array[i]);
		while (mpz_cmp_ui(modArray[i], 1) > 0) {	//Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
			//mpz_mul(modArray[i], modArray[i], modArray[i]); //possible speedup
			mpz_gcd(modArray[i], modArray[i], smoothnessArray[i]); 
			mpz_divexact(smoothnessArray[i], smoothnessArray[i], modArray[i]);
		}
		if (mpz_cmp_ui(smoothnessArray[i], 1) == 0 || mpz_cmp_ui(modArray[i], 0) == 0) {
			mpz_set(smoothnessArray[*numberSmoothElements], array[i]);
			(*numberSmoothElements)++;
		}
	}
}

//test segments in array for smoothness given the product of all smooth primes FrKlMoWi
void smoothBatchFrKlMoWiProductSegments (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, unsigned int numberSegments, unsigned int segmentSize, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothSegments, unsigned int *smoothSegments) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
	unsigned int i;
	unsigned int j;
	unsigned int current;
	unsigned short smooth;
	mpz_t temp; mpz_init(temp);
	
//compute products and modulo
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute FrKlMoWi smoothness algorithm
	*numberSmoothSegments = 0;
	for (i = 0; i < numberSegments; i++) {
		smooth = 1;
		current = i * segmentSize;
		for (j = 0; j < segmentSize; j++) {
			mpz_set(temp, array[current]);
			while (mpz_cmp_ui(modArray[current], 1) > 0) {	//Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
				//mpz_mul(modArray[current], modArray[current], modArray[current]); //possible speedup
				mpz_gcd(modArray[current], modArray[current], temp); 
				mpz_divexact(temp, temp, modArray[current]);
			}
			if (mpz_cmp_ui(temp, 1) != 0 && mpz_cmp_ui(modArray[current], 0) != 0) {
				smooth = 0;
				break;
			}
			current++;
		}
		if (smooth == 1) {
			smoothSegments[*numberSmoothSegments] = i;
			(*numberSmoothSegments)++;
		}
	}
	
//clear and free
	mpz_clear(temp);
}

//test elements in array for smoothness given a list of all smooth primes Bernstein 
void smoothBatchBernsteinList (unsigned int numberSmoothPrimes, unsigned int *smoothPrimes, unsigned int arrayLength, mpz_t *array, unsigned int *numberSmoothElements, mpz_t **smoothnessArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
//initialize arrays
	mpz_t productPrimes; mpz_init(productPrimes);
	mpz_t *productTreeArray = (mpz_t *) malloc((2 * arrayLength - 1) * sizeof(mpz_t));
	if (productTreeArray == NULL) {
		printf("Error: productTreeArray could not be allocated\n");
		exit(-1);
	}
	mpz_t *modArray = (mpz_t *) malloc((arrayLength) * sizeof(mpz_t));
	if (modArray == NULL) {
		printf("Error: modArray could not be allocated\n");
		exit(-1);
	}
	mpz_t *tempSmoothArray = (mpz_t *) malloc(arrayLength * sizeof(mpz_t));
	if (tempSmoothArray == NULL) {
		printf("Error: tempSmoothArray could not be allocated\n");
		exit(-1);
	}
	unsigned int i;
	mpz_t powerOfTwo; mpz_init_set_ui(powerOfTwo, 2);
	for (i = 0; i < arrayLength; i++){
		mpz_init(modArray[i]);
		mpz_init(tempSmoothArray[i]);
		mpz_init(productTreeArray[i]);
	}
	for (i = arrayLength; i < 2 * arrayLength - 1; i++){
		mpz_init(productTreeArray[i]);
	}
	
	//compute products and modulo
	treeProduct_ui(numberSmoothPrimes, smoothPrimes, &productPrimes);	//compute product of smooth primes
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute Bernstein smoothness algorithm
	*numberSmoothElements = 0;
	while (mpz_cmp(array[arrayLength - 1], powerOfTwo) > 0) { //last element in array is the largest. Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
		for (i = 0; i < arrayLength; i++) {
			//if (mpz_cmp_ui(array[i], powerOfTwo) > 0) {	//elements have roughly the same size
				mpz_mul(modArray[i], modArray[i], modArray[i]);
				mpz_mod(modArray[i], modArray[i], array[i]);
			//}
		}
		mpz_mul(powerOfTwo, powerOfTwo, powerOfTwo);
	}
	for (i = 0; i < arrayLength; i++) {
		if (mpz_cmp_ui(modArray[i], 0) == 0) {
			mpz_set(tempSmoothArray[*numberSmoothElements], array[i]);
			(*numberSmoothElements)++;
		}
	}

//write smooth elements in shorter array
	if (*numberSmoothElements > 0) {
		*smoothnessArray = (mpz_t *) malloc(*numberSmoothElements * sizeof(mpz_t));
		if (smoothnessArray == NULL) {
			printf("Error: smoothnessArray could not be allocated\n");
			exit(-1);
		}
		for (i = 0; i < *numberSmoothElements; i++) {
			mpz_init_set((*smoothnessArray)[i], tempSmoothArray[i]); //remember to clear them in main function!
		}
	}
	
//clear and free
	mpz_clear(productPrimes);
	mpz_clear(powerOfTwo);
	for (i = 0; i < arrayLength; i++){
		mpz_clear(modArray[i]);
		mpz_clear(tempSmoothArray[i]);
		mpz_clear(productTreeArray[i]);
	}
	for (i = arrayLength; i < 2 * arrayLength - 1; i++){
		mpz_clear(productTreeArray[i]);
	}
	free(productTreeArray);
	free(modArray);
	free(tempSmoothArray);
}

//test elements in array for smoothness given the product of all smooth primes Bernstein
void smoothBatchBernsteinProduct (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothElements, mpz_t *smoothnessArray) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	unsigned int i;
	mpz_t powerOfTwo; mpz_init_set_ui(powerOfTwo, 2);
	
	//compute product and modulo
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute Bernstein smoothness algorithm
	*numberSmoothElements = 0;
	while (mpz_cmp(array[arrayLength - 1], powerOfTwo) > 0) { //last element in array is the largest. Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
		for (i = 0; i < arrayLength; i++) {
			//if (mpz_cmp_ui(array[i], powerOfTwo) > 0) {	//elements have roughly the same size
				mpz_mul(modArray[i], modArray[i], modArray[i]);
				mpz_mod(modArray[i], modArray[i], array[i]);
			//}
		}
		mpz_mul(powerOfTwo, powerOfTwo, powerOfTwo);
	}
	for (i = 0; i < arrayLength; i++) {
		if (mpz_cmp_ui(modArray[i], 0) == 0) {
			mpz_set(smoothnessArray[*numberSmoothElements], array[i]);
			(*numberSmoothElements)++;
		}
	}
	
//clear and free
	mpz_clear(powerOfTwo);
}

//test segments in array for smoothness given the product of all smooth primes Bernstein
void smoothBatchBernsteinProductSegments (mpz_t productPrimes, unsigned int arrayLength, mpz_t *array, mpz_t largestElement, unsigned int numberSegments, unsigned int segmentSize, mpz_t *productTreeArray, mpz_t *modArray, unsigned int *numberSmoothSegments, unsigned int *smoothSegments) {
	
	if (arrayLength == 0) {
		printf("Error: Array is empty\n");
		exit(-1);
	}
	
	unsigned int i;
	unsigned int j;
	unsigned int current;
	unsigned short smooth;
	mpz_t powerOfTwo; mpz_init_set_ui(powerOfTwo, 2);
	
	//compute product and modulo
	treeProductSave(arrayLength, array, productTreeArray);	//compute product of array and save intermediate steps
	treeMod(productPrimes, arrayLength, productTreeArray, modArray);	//compute product of smooth primes modulo elements in array using tree
	
//compute Bernstein smoothness algorithm
	while (mpz_cmp(largestElement, powerOfTwo) > 0) { //Optional: set maximal number of repititions for efficiency. This results in a different (stronger) type of smoothness (cf. powersmooth for sieve).
		for (i = 0; i < arrayLength; i++) {
			//if (mpz_cmp_ui(array[i], powerOfTwo) > 0) {	//elements have roughly the same size
				mpz_mul(modArray[i], modArray[i], modArray[i]);
				mpz_mod(modArray[i], modArray[i], array[i]);
			//}
		}
		mpz_mul(powerOfTwo, powerOfTwo, powerOfTwo);
	}
	
	*numberSmoothSegments = 0;
	for (i = 0; i < numberSegments; i++) {
		smooth = 1;
		current = i * segmentSize;
		for (j = 0; j < segmentSize; j++) {
			if (mpz_cmp_ui(modArray[current], 0) != 0) {
				smooth = 0;
				break;
			}
			current++;
		}
		if (smooth == 1) {
			smoothSegments[*numberSmoothSegments] = i;
			(*numberSmoothSegments)++;
		}
	}
	
//clear and free
	mpz_clear(powerOfTwo);
}
