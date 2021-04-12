# install.packages('AlphaSimR')
library(AlphaSimR)

SP = SimParam$
        new(founderPop)$
        addTraitA(1000)$
        setVarE(H2=0.4)

SP$setVarE(H2=0.2)


# founder
founderPop = runMacs(nInd = 50, # number of individuals
                     nChr = 21, # number of chromosomes
                     segSites = 1000, # number of segregating sites
                     inbred=TRUE) # should founder be inbred
Parents = newPop(founderPop)

# 1st year
F1 = randCross(Parents, 200) # 200 ind

# 2nd/3rd year
# head row (HDRW) nursery
HDRW = makeDH(F1, 100) # 100 DH lines
HDRW = setPheno(HDRW, H2=0.1) # visual selection

# 4th year
PYT = selectWithinFam(HDRW, 5) # best 5 in familiy
PYT = setPheno(PYT)

# 5th year
AYT = selectInd(PYT, 100) # best 100 individuals
AYT = setPheno(AYT, reps=4) # select at 4 locations

# 6th year
EYT = selectInd(AYT, 10) # best 10 lines
EYT = setPheno(EYT, reps=16) # select at 16 locations

# 7th year
Variety = selectInd(EYT, 1)

yield=list(Parents=gv(Parents),
            F1=gv(F1),
            HDRW=gv(HDRW),
            PYT=gv(PYT),
            AYT=gv(AYT),
            EYT=gv(EYT),
            Variety=gv(Variety))
boxplot(yield, ylab="GeneticValue")


F1 = randCross(Parents, 200) # 200 ind

PYT = selectWithinFam(HDRW, 5) # best 5 in familiy

AYT = selectInd(PYT, 100) # best 100 individuals