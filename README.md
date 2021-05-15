# XSim

[![Build Status](https://travis-ci.org/reworkhow/XSim.jl.svg?branch=master)](https://travis-ci.org/reworkhow/XSim.jl)

XSim is a fast and user-friendly tool to simulate sequence data and complicated pedigree structures

#### Features

* An efficient CPOS algorithm
* Using founders that are characterized by real genome sequence data
* Complicated pedigree structures among descendants

#### Quick-start

```Julia
# Load XSim package
using XSim

# Set genome information
build_genome(species="Pig")
n_qtl = [5, 10]
Vg    = [ 1 .5
         .5  1]
build_phenome(n_qtl, Vg)

# Generate founders
popSizeFounder = 50
sires = Cohort(popSizeFounder);
dams  = Cohort(popSizeFounder);

#random mating
h2       = [0.5, 0.3]
weights  = [1.0, 0.0]
n        = 100
n_select = 10

f1 = random_mate(sires, dams, n)
f1 = select(f1, n_select, h2=h2, weights=weights)

f2 = self_mate(f1, n)
f2 = select(f2, n_select, h2=h2, weights=weights)

f3 = self_mate(f2, n)
f3 = select(f3, n_select, h2=h2, weights=weights)

sum(get_BVs(f1), dims=1)
sum(get_BVs(f2), dims=1)
sum(get_BVs(f3), dims=1)

```


#### More

* **homepage**: [QTL.rocks](http://QTL.rocks)
* **Installation**: at the Julia REPL, `Pkg.add("XSim")`
* **Documentation**: [available here](https://github.com/reworkhow/XSim.jl/wiki)
* **Authors**: Hao Cheng,Rohan Fernando,Dorian Garrick
* **Citing XSim** 

>Cheng H, Garrick D, and Fernando R (2015) XSim: Simulation of descendants from ancestors with sequence data. G3: Genes-Genomes-Genetics, 5(7):1415-1417.
