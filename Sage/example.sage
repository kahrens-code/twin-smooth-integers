#Setup
load("PTESolutions.sage")
load("sieveEratosthenes.sage")
load("naiveApproach.sage")
load("newApproach.sage")

#define the PTE solution
a = a1;
b = b1;
C = C1;
roots = roots1;
residues = findResidues(a,C);
print("a(x) =", a(x));
print("b(x) =", b(x));
print("C =", C);

#mark all integers in the interval [2^40,2^40+2^20) as 2^16-power-smooth (1) or not-2^16-power-smooth (0)
smoothnessBound = 2^16;
start = 2^40;
smoothPrimesList = smoothPrimes(smoothnessBound);
smoothList40 = sieve(start,2^20,smoothPrimesList);
print("smoothness bound:", smoothnessBound);
print("start:", start);
print("interval size: 2^20");

#use the naive or the new algorithm
print("naive:");
print(checkAll(smoothList40,roots,a,C));
print("new:");
print(checkRes(smoothList40,residues,C,roots));
