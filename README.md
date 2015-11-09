# XSim

[![Build Status](https://travis-ci.org/reworkhow/GenSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/GenSim.jl)

XSim is a fast and user-friendly tool to simulate sequence data and complicated pedigree structures

####Features

* An efficient CPOS algorithm
* Using founders that are characterized by real genome sequence data
* Complicated pedigree structures among descendants

####Quick-start

```Julia
using XSim

#set genome information
chrLength = 0.1
numChr    = 5
numLoci   = 10
mutRate   = 0.0
locusInt  = chrLength/numLoci
mapPos    = collect(locusInt/2:locusInt:chrLength)
geneFreq  = fill(0.5,numLoci)
qtlMarker = fill(false,numLoci); qtlMarker[sample(1:numLoci)]=true
qtlEffects= randn(numLoci)
XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,qtlMarker,qtlEffects,mutRate)

#generate founders
popSizeFounder = 10
sires = sampleFounders(popSizeFounder)
dams  = sampleFounders(popSizeFounder)

#random mating
ngen,popSize = 10,10
siresF, damsF, genF = sampleRan(popSize, nGen, sires, dams, gen=gen);

#sample following a pedigree
nf = 3  #number of halfsib families 
ns = 1  #number of sires
nd = 2  #number of dams per sire 
no = 1  #number of offspring per dam
ind = collect(1:(ns+ns*nd+ns*nd*no))
sire = int([fill(0,ns+ns*nd),kron([1:ns],ones(nd*no))])
dam  = int([fill(0,ns+ns*nd),kron([1+ns:ns+ns*nd],ones(no))])
pedSim  = [ind sire dam]

pedArray = Array(XSim.PedNode,size(pedSim,1))
for i in 1:size(pedSim,1)
    indi  = pedSim[i,1]
    sirei = pedSim[i,2]
    dami  = pedSim[i,3]
    pedArray[pedSim[i,1]] = XSim.PedNode(indi,sirei,dami)
end

pedFounders = sampleFounders(popSizeFounder)
animals = XSim.samplePed(pedArray,pedFounders);

#concat several cohorts
concatCohorts(siresF,damsF)

#generate populations
ngen,popSize    = 10,10


##
pop1.popSample(ngen,popSize)
pop2 = pop1.popNew(10);
pop3 = popCross(5,pop1,pop2);

#generate genotypes
M = pop3.getGenotypes()
```

####More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.clone("https://github.com/reworkhow/XSim.jl.git")`
* **Documentation**: [available here](http://xsimjl.readthedocs.org/en/latest/)
* **Authors**: [Hao Cheng](http://reworkhow.github.io),[Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan), [Dorian Garrick](http://www.ans.iastate.edu/faculty/index.php?id=dorian)
* **Citing XSim** 

Cheng H, Garrick D, and Fernando R (2015) Xsim: Simulation of descendants from ancestors with sequence data. G3: Genes-Genomes-Genetics, 5(7):1415-1417.
