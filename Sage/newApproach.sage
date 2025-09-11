#New approach to find twin smooth primes using a single PTE solution


#find those residue classes [r] such that a(r+k*C)/C and b(r+k*C)/C are integers

#polyA is one of the two polynomials corresponding to the PTE solution
#modC is the difference C=|a(x)-b(x)| between the two polynomials correspondingto the PTE solution
#the output is an ordered list of representatives from the interval [0,C)

def findResidues(polyA,modC):
    facC=factor(modC);numFac=len(facC);
    primePowers=[facC[i][0]^facC[i][1] for i in range(numFac)];     #find the prime factorisation of modC
    primeRes=[];
    for q in primePowers:                                           #for each prime factor
        partialRes=[];                                              #make a list of relvant residue classes mod q
        for i in range(q):                                          #where q=p^e is the maximal power of that prime dividing modC
            if ZZ(polyA(i))%q==0: partialRes.append(i)
        primeRes.append(partialRes)                                 #and add it to a list of lists
    numComb=prod(len(primeRes[i]) for i in range(numFac))           #find the number of all possible combinations with one entry of
    resList=[[] for i in range(numComb)]                            #each sublist and prepare a list to write them down
    mult1=numComb; mult2=1; numCurrent=1;                           #make sure every combination is created
    for i in range(numFac): 
        mult2=mult2*numCurrent;
        numCurrent=len(primeRes[i]);
        mult1=ZZ(mult1/numCurrent); 
        for j in range(mult2): 
            for k in range(numCurrent): 
                for l in range(mult1): 
                    resList[mult1*numCurrent*j+mult1*k+l].append(primeRes[i][k])
    residues=[CRT(resList[i],primePowers) for i in range(numComb)]  #finally use the CRT to get the relevant residues mod modC
    return(sorted(residues))                                        #and order them by "<"


#check if representatives are power-smooth

#smoothList is a list of integers marked as smooth (1) or non-smooth (0)
#residues is a list of residue classes [r] such that a(r+k*C)/C and b(r+k*C)/C are integers or rather of the corresponding integers r
#modC is the difference C=|a(x)-b(x)| between the two polynomials correspondingto the PTE solution
#roots is a list of all elements in the PTE solution or equivalentely all zeros of both of the corresponding polynomials
#the output is a list of integers l for which a(l)/C and b(l)/C are twin smooth integers

#if C<size there are multiple representatives of each residue class in an interval of length "size"
def checkRes(smoothList,residues,modC,roots):
    startModC=smoothList[0][0]%modC; size=len(smoothList); minimum=max(roots);
    values=[]; numRoots=len(roots);
    for resid in residues:
        current=resid-startModC;                                    #smoothList[0][0]-startModC is 0 mod modC 
        while current<minimum: current=current+modC;                #prevent current-roots[r]<0
        while current<size:
            if smoothList[current][1]==1:                           #check if smoothList[current][0] is smooth
                for r in range(1,numRoots):                         #roots[0] is 0 (already checked)
                    if smoothList[current-roots[r]][1]==0: break;   #check if all factors of a(smoothList[current][1]) are smooth
                else: values.append(smoothList[current][0]);        #add smoothList[current][0] to the values, if there was no "break" 
            current=current+modC;                                   #proceede with the next representative
    return(values)

#if C>=size there is at most one representative of each residue class in an interval of length "size"
def checkRES(smoothList,residues,modC,roots):
    startModC=smoothList[0][0]%modC; size=len(smoothList);
    minimum=max(roots); numRes=len(residues);
    values=[]; numRoots=len(roots);
    searchL=0;searchR=numRes;searchPivot=floor(numRes/2);           #find the smallest representative (of all classes) in the interval
    while searchR-searchL>0:                                        #using a bisection method
        if startModC>residues[searchPivot]: searchL=searchPivot+1;
        else: searchR=searchPivot
        searchPivot=floor((searchL+searchR)/2)
    step=searchPivot;current=residues[step]-startModC;
    while minimum>current:                                          #prevent current-roots[r]<0
        step=step+1;
        if step<numRes: current=residues[step]-startModC;
        else:
            step=0;
            startModC=startModC-modC;
            current=-startModC;
    for i in range (step,numRes):                                   #check all representatives starting from the one we found
        current=residues[i]-startModC;
        if current<size:                                            #if we are still in the interval
            if smoothList[current][1]==1:                           #check for smoothness
                for r in range(1,numRoots):
                    if smoothList[current-roots[r]][1]==0: break;
                else: values.append(smoothList[current][0]);
        else: return(values)                                        #else return the list of values producing twin smooth integers
    startModC=startModC-modC;                                       #if we reached the last residue class before the end of the
    for i in range (0,step):                                        #interval, we loop around and continue with r+modC
        current=residues[i]-startModC;                              #=residues[i]-startModC+modC
        if current<size:
            if smoothList[current][1]==1:
                for r in range(1,numRoots):
                    if smoothList[current-roots[r]][1]==0: break;
                else: values.append(smoothList[current][0]);
        else: return(values)


#if we prefere one function for both cases, we can define the following
def checkC(smoothList,residues,modC,roots):
    if modC>=len(smoothList):
        return(checkRES(smoothList,residues,modC,roots))
    else: return(checkRes(smoothList,residues,modC,roots))
