# Twin Smooth Integers

Sieving for large twin smooth integers using single solutions to the Prouhet-Tarry-Escott (PTE) problem.
Optimized for single solutions by doing the modulo calculus first.

The Sieve of Eratosthenes marks all integers in an interval as smooth (1) or non-smooth (0).
The Naive Approach gives a list of integers l producing twin smooth integers a(l)/C,b(l)/C by checking every smooth integer in the interval.
The Optimised Approach gives a list of integers l producing twin smooth integers a(l)/C,b(l)/C by checking only those integers l such that a(l) and b(l) are 0 modulo C.

More details can be found in the paper ***comming soon***.
