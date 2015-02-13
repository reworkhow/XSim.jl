# GenSim

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














[![Build Status](https://travis-ci.org/reworkhow/GenSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/GenSim.jl)
