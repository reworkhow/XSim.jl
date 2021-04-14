using Test
using XSim

function t1()
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)
	@test XSim.common.G.numChrom == numChr
	@test length(XSim.common.G.chr) == numChr
	@test XSim.common.countChromosome== 1
end
