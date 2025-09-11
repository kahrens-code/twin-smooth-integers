# Twin Smooth Integers

Sieving for large twin smooth integers using single solutions to the Prouhet-Tarry-Escott (PTE) problem.
Optimised for single solutions by doing the modulo computation first.

Ideal PTE solutions provide polynomials a(x) and b(x) such that |a(x) - b(x)| = C is an interger.

To search for twin smooth integers (i.e. n + 1 and n - 1 are smooth) one can define a search interval I and use one of the following approaches.

Naive approach:

    For every element l in I:
        if a(l) and b(l) are smooth:
            if a(l) / C and b(l) / C are integers:
                print l

The tree approach by Costello, Meyer and Naehrig combined many PTE solutions and used common factors of the polynomials.

New approach (for single solutions):

    For every element l in I:
        if a(l) / C and b(l) / C are integers:
            if a(l) and b(l) are smooth:
                print l
                
This decreases the number of look-ups for smoothness and allows to do a large part of the "modulo C" calculations in a pre-computation, that can be used for all search intervals.

For large integers this becomes inefficient and we can use a third approach with batch smoothness tests. This requires large integers and is only implemented in C with GMP.

Fixed residue approach:

	Choose r such that a(r) and b(r) are integers
	For every integer k between start and end:
		if a(r + k C) and b(r + k C) are smooth:
                print r + k C

There is a Sage implementation for explanatory reasons, two C implementations with maximal variable sizes of 64 and 128 bit, respectively, and a C implementation using GMP to allow for arbitrary sizes. For explanations of the individual steps please see the commentary in the code or the paper.
The list of used and/or referenced PTE solutions can be found in the parameter files.

More details and references can be found in the paper "Sieving for Large Twin Smooth Integers Using Single Solutions to Prouhet-Tarry-Escott" https://eprint.iacr.org/2023/219/.
