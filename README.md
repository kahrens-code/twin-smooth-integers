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

Optimised approach (for single solutions):

    For every element l in I:
        if a(l) / C and b(l) / C are integers:
            if a(l) and b(l) are smooth:
                print l
                
This decreases the number of look-ups for smoothness and allows to do a large part of the "modulo C" calculations in a pre-computation, that can be used for all search intervals.

There is a Sage implementation for explanatory reasons and two C implementations with maximal variable sizes of 64 and 128 bit respectively. For explanation of the individual steps please see the commentary in the code.
The list of used and/or referenced PTE solutions can be found in the parameters.c file for the 128 bit C implementation.

More details and references can be found in the paper https://eprint.iacr.org/2023/219/.
