# GenSim

GenSim is a fast and user-friendly tool to simulate sequence data and complicated pedigree structures


####Features


* An efficient CPOS algorithm
* Using founders that are characterized by real genome sequence data
* Complicated pedigree structures among descendants

####algorithm behind

* A CPOS algorithm is implemented to efficiently simulate sequence data and complicated pedigree structures.

####Quick-start

```Julia
using GenSim

#set genome information
chrLength, numChr, numLoci, mutRate = 1.0, 1, 100, 0.0
locusInt  = chrLength/numLoci
mapPos    = [0:locusInt:(chrLength-0.0001)]
geneFreq  = fill(0.5,numLoci)

GenSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,mutRate)
jsim = GenSim.startPop()

#generate populations
ngen,popSize    = 10,10
jsim.popSample(ngen,popSize)

#generate genotypes
M=jsim.getGenotypes()
```

####Authors and Contributors

* Hao Cheng, Rohan Fernando and Dorian Garrick


[![Build Status](https://travis-ci.org/reworkhow/GenSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/GenSim.jl)
