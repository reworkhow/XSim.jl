using XSim
using LinearAlgebra
using Test

# write your own tests here

@testset "building_genomes_append" begin
	#
	#build genomes appends to the existing XSim.common
	# this may or may not pass depending how often build genomes has been called
	# hence we call clearGlobals at the start 
	#
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

@testset "building_genomes_with_founders" begin
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)

	popSizeFounder = 2
	sires = sampleFounders(popSizeFounder);
	dams  = sampleFounders(popSizeFounder);
	@test typeof(sires)==XSim.Cohort
	@test typeof(dams)==XSim.Cohort
	
	ploidy = 2
	@test length(XSim.common.founders) == 2+2 #sires & dams
	@test XSim.common.countChromosome == 1 + ploidy*(popSizeFounder+popSizeFounder)	
end

@testset "Random Mating" begin
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)

	popSizeFounder = 2
	sires = sampleFounders(popSizeFounder);
	dams  = sampleFounders(popSizeFounder);

	ngen,popSize = 5,10
	sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);
	@test typeof(sires1)==XSim.Cohort
	@test typeof(dams1)==XSim.Cohort
end
#
# This testset is designed to show how the selectSel function
# can be called with a correctly initilized set of globals for 1 trait
# Basically call in order
# 0. clearGlobals
# 1. build_genome
# 2. sampleFounders
# 3. XSim.setResidualVariance
# 4. sampleSel
#
@testset "selectSel_goldenpath_oneTrait" begin
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)

	popSize =10
	ngen =5
	popSizeFounder = 2
	sires = sampleFounders(popSizeFounder);
	dams  = sampleFounders(popSizeFounder);

	# Initilize LRes for 1 trait
	A = fill(1.0,1,1)
	XSim.setResidualVariance(cholesky(A).U')

	@test_throws BoundsError sire1,dam1,gen1=sampleSel(popSize, 4, 4, ngen,sires, dams);
end
#
# run a sampleSel with two traits
#
@testset "selectSel_goldenpath_twoTrait" begin
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)

	popSize = 10
	ngen = 5
	popSizeFounder = 2
	sires = sampleFounders(popSizeFounder);
	dams  = sampleFounders(popSizeFounder);

	# initialize for 2 trait
	A = [10.0 0 ; 0 3.0]
	XSim.setResidualVariance(cholesky(A).U')
	
	@test_throws DimensionMismatch sire1,dam1,gen1=sampleSel(popSize, 4, 4, ngen,sires, dams);
end

#
# This testset is designed to show how not initializing the globals
# causes various failure modes.
#
@testset "Mating with selection" begin
	clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	popSize =10
	ngen =5
	popSizeFounder = 2
	sires = sampleFounders(popSizeFounder);
	dams  = sampleFounders(popSizeFounder);

	# if LRes is not initialized 
	@test_throws ErrorException sire2,dam2,gen2=sampleSel(popSize, 4, 4, ngen,sires, dams);

	# Once LRes has been initialzed for 1 trait
	A = fill(1.0,1,1)
	XSim.setResidualVariance(cholesky(A).U')
	
	# this happens as the QTL have not been initialized...
	@test_throws DimensionMismatch sire2,dam2,gen2=sampleSel(popSize, 4, 4, ngen,sires, dams);
end

