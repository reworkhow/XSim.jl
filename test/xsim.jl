pwd()
# cd("src")
include("XSim.jl")
using XSim

chrLength = 1.0
numChr    = 1
numLoci   = 2000
mutRate   = 0.0
locusInt  = chrLength/numLoci
mapPos   = collect(0:locusInt:(chrLength-0.0001))
geneFreq = fill(0.5, numLoci)
build_genome(numChr,
             chrLength,
             numLoci,
             geneFreq,
             mapPos,
             mutRate)


popSizeFounder = 2
sires = sampleFounders(popSizeFounder);
dams  = sampleFounders(popSizeFounder);


ngen, popSize = 5,10
sires1, dams1, gen1 = sampleRan(popSize, ngen, sires, dams);





sires1.animalCohort[1].genomePat

nSires,nDams = 2,2
popSize,ngen = 10,5
varRes = 1.0
sire2,dam2,gen2 = sampleSel(popSize, nSires, nDams, ngen,sires, dams, varRes);


# two separate selct : select, and mating

sires1=select(num_selected_ind,sires)
dams1 =select(num_selected_ind,dams)
sires2,dams2=random_mating(noffspring,sires1,dams1)

function random_mating(noffspring,sires1)
    sires2,dams2=random_mating(noffspring,sires1,sires1)
    ind3=concat_cohort(sires2,dams2)
end

sires3,dams3=random_mating(noffspring,ind3,ind3)
ind4=concat_cohort(sires3,dams3)