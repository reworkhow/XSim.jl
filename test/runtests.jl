using XSim
using Test

# write your own tests here
chrLength= 0.1  #length of each chromosome
numChr   = 2    #number of chromosomes
nmarkers = 10   #number of loci for each chromosome
nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
build_genome(numChr,chrLength,nmarkers,nQTL)

popSizeFounder = 2
sires = sampleFounders(popSizeFounder);
dams  = sampleFounders(popSizeFounder);

#random mating
ngen,popSize = 5,10
sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);

@test typeof(sires1)==XSim.Cohort
@test typeof(dams1)==XSim.Cohort
