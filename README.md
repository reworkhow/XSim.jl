# XSim

[![Build Status](https://travis-ci.org/reworkhow/XSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/XSim.jl)

XSim is a fast and user-friendly tool to simulate sequence data and complicated pedigree structures

#### Features

* An efficient CPOS algorithm
* Using founders that are characterized by real genome sequence data
* Complicated pedigree structures among descendants

#### Quick-start

```Julia
#load XSim package
using XSim

#set genome information
chrLength= 0.1  #length of each chromosome 
numChr   = 2    #number of chromosomes
nLoci    = 10   #number of loci for each chromosome
nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
build_genome(numChr,chrLength,nLoci,nQTL) #this genome information will be used for subsequent computaions

#generate founders
popSizeFounder = 2
sires = sampleFounders(popSizeFounder);
dams  = sampleFounders(popSizeFounder);

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
* Run `Pkg.add(PackageSpec(name="XSim", rev="master"))` to get the newest unofficial JWAS. Run `Pkg.free("XSim")` to go back to the official one.

>Cheng H, Garrick D, and Fernando R (2015) XSim: Simulation of descendants from ancestors with sequence data. G3: Genes-Genomes-Genetics, 5(7):1415-1417.
