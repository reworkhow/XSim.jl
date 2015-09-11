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
chrLength, numChr, numLoci, mutRate = 1.0, 1, 100, 0.0
locusInt  = chrLength/numLoci
mapPos    = [0:locusInt:(chrLength-0.0001)];
geneFreq  = fill(0.5,numLoci);

XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,mutRate)
pop1 = startPop()

#generate populations
ngen,popSize    = 10,10

pop1.popSample(ngen,popSize)
pop2 = pop1.popNew(10);
pop3 = popCross(5,pop1,pop2);

#generate genotypes
M = pop3.getGenotypes()
```

####More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.add("JWAS")`
* **Documentation**: [available here](http://xsimjl.readthedocs.org/en/latest/)
* **Authors**: [Hao Cheng](http://reworkhow.github.io),[Rohan Fernando](http://www.ans.iastate.edu/faculty/index.php?id=rohan), [Dorian Garrick](http://www.ans.iastate.edu/faculty/index.php?id=dorian)



