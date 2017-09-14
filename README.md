# XSim

[![Build Status](https://travis-ci.org/reworkhow/XSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/XSim.jl)

XSim is a fast and user-friendly tool to simulate sequence data and complicated pedigree structures

#### Features

* An efficient CPOS algorithm
* Using founders that are characterized by real genome sequence data
* Complicated pedigree structures among descendants

#### Quick-start

```Julia
using XSim
#set genome information
using StatsBase
numChr,numLoci,chrLength,mutRate = 2,10,0.1,0.0
mapPos     = collect(0.005:0.01:0.1)
geneFreq   = fill(0.5,numLoci)
qtlMarker  = fill(false,numLoci)
qtlMarker[sample(1:numLoci)]= true
qtlEffects = randn(numLoci)
XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,qtlMarker,qtlEffects,mutRate)

#generate founders
popSizeFounder = 2
sires = sampleFounders(popSizeFounder)
dams  = sampleFounders(popSizeFounder)

#random mating
ngen,popSize = 5,10
sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);
```

#### More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.add("XSim")`
* **Documentation**: [available here](https://github.com/reworkhow/XSim.jl/wiki)
* **Authors**: Hao Cheng,Rohan Fernando,Dorian Garrick
* **Citing XSim** 

>Cheng H, Garrick D, and Fernando R (2015) XSim: Simulation of descendants from ancestors with sequence data. G3: Genes-Genomes-Genetics, 5(7):1415-1417.
