using XSim
using Base.Test

# write your own tests here
using StatsBase
numChr,numLoci,chrLength,mutRate = 2,10,0.1,0.0
mapPos     = collect(0.005:0.01:0.1)
geneFreq   = fill(0.5,numLoci)
qtlMarker  = fill(false,numLoci)
qtlMarker[sample(1:numLoci)]= true
qtlEffects = randn(numLoci)
XSim.init(numChr,numLoci,chrLength,geneFreq,mapPos,qtlMarker,qtlEffects,mutRate);

popSizeFounder = 2
sires = sampleFounders(popSizeFounder);
dams  = sampleFounders(popSizeFounder);

#random mating
ngen,popSize = 5,10
sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);

@test typeof(sires1)==XSim.Cohort
@test typeof(dams1)==XSim.Cohort
