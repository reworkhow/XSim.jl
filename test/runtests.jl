using XSim
using LinearAlgebra
using Test

# write your own tests here
@testset "building_genomes_append" begin
	#
	#build genomes appends to the existing XSim.GLOBAL
	# this may or may not pass depending how often build genomes has been called
	# hence we call clear_globals at the start 
	#
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)
	@test XSim.GLOBAL.G.numChrom == numChr
	@test length(XSim.GLOBAL.G.chr) == numChr
	@test XSim.GLOBAL.countChromosome== 1
end

@testset "building_genomes_with_founders" begin
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr, chrLength, nmarkers, nQTL)

	popSizeFounder = 2
	sires = Cohort(popSizeFounder);
	dams  = Cohort(popSizeFounder);
	@test typeof(sires)==XSim.Cohort
	@test typeof(dams)==XSim.Cohort

	ploidy = 2
	@test length(XSim.GLOBAL.founders) == 2 + 2 #sires & dams
	@test XSim.GLOBAL.countChromosome == 1 + ploidy*(popSizeFounder+popSizeFounder)	
end

@testset "Random Mating" begin
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	chrLength= 0.1  #length of each chromosome
	numChr   = 2    #number of chromosomes
	nmarkers = 10   #number of loci for each chromosome
	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
	build_genome(numChr,chrLength,nmarkers,nQTL)

	popSizeFounder = 2
	sires = Cohort(popSizeFounder);
	dams  = Cohort(popSizeFounder);

	ngen,popSize = 5,10
	sires1, dams1 = sample_random(sires, dams, popSize, ngen);
	@test typeof(sires1)==XSim.Cohort
	@test typeof(dams1)==XSim.Cohort
end
#
# This testset is designed to show how the selectSel function
# can be called with a correctly initilized set of globals for 1 trait
# Basically call in order
# 0. clear_globals
# 1. build_genome
# 2. Cohort
# 3. sampleSel
#
@testset "selectSel_goldenpath_oneTrait" begin
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	numChr=3
	chrLength = [1.0, 1.1, 0.9]
	numLoci=[10 ,15 ,20]
	geneFreq   = Array{Array{Float64,1},1}(undef,0)
	qtlIndex  = Array{Array{Int64,1},1}(undef,0)
	qtlEffects = Array{Array{Float64,2},1}(undef,0)
	numQTLOnChr =[2, 3, 1]
	numQTL=6
	qtlIndex = [[1; 4],[3;7;9],[5]]
	geneFreq = [[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5]]
	mPos = [[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5],
			[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5;.55;.65;.75;.85;.95],
			[.05;.1;.15;.2;.25;.3;.35;.4;.45;.5;.54;.58;.65;.75;.758;.759;.8;.85;.851;.859]]
	push!(qtlEffects,randn(2,1))
	push!(qtlEffects,randn(3,1))
	push!(qtlEffects,randn(1,1))
	G0 = Array{Float64,2}(undef,1,1)
	G0[1,1] = 1
	nTraits = 1
	XSim.build_genome(numChr,chrLength,numLoci,geneFreq,mPos,qtlIndex,qtlEffects,nTraits,G0)
	XSim.set_residual_var(G0)

	popSize =10
	ngen =5
	popSizeFounder = 20
	sires = Cohort(popSizeFounder);
	dams  = Cohort(popSizeFounder);

	sire1, dam1 = sample_select(sires, dams, popSize, 4, 4, ngen);
end
#
# run a sampleSel with two traits
# based on the wiki multitrait example
#
@testset "selectSel_goldenpath_twoTrait" begin
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	numChr = 3
	chrLength = [1.0, 1.1, 0.9]
	numLoci = [10 ,15 ,20]
	nTraits = 2
	geneFreq   = Array{Array{Float64,1},1}(undef,0)
	qtlIndex  = Array{Array{Int64,1},1}(undef,0)
	qtlEffects = Array{Array{Float64,2},1}(undef,0)
	numQTLOnChr =[2, 3, 1]
	numQTL=sum(numQTLOnChr)
	qtlIndex = [[1; 4],[3;7;9],[5]]
	geneFreq = [[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5]]
	mPos = [[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5],
			[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5;.55;.65;.75;.85;.95],
			[.05;.1;.15;.2;.25;.3;.35;.4;.45;.5;.54;.58;.65;.75;.758;.759;.8;.85;.851;.859]]
	push!(qtlEffects,randn(2,nTraits))
	push!(qtlEffects,randn(3,nTraits))
	push!(qtlEffects,randn(1,nTraits))
	G0 = [1. 0.5;0.5 2.]
	XSim.build_genome(numChr,chrLength,numLoci,geneFreq,mPos,qtlIndex,qtlEffects,nTraits,G0)
	XSim.set_residual_var(G0)

	popSize = 10
	ngen = 3
	popSizeFounder = 20
	sires = Cohort(popSizeFounder)
	dams  = Cohort(popSizeFounder)

	sire1, dam1 = sample_select(sires, dams, popSize, 4, 4, ngen);
end

#
# This testset is designed to show how not initializing the globals
# causes various failure modes.
#
@testset "Mating with selection" begin
	clear_globals()
	@test XSim.GLOBAL.G.numChrom == 0
	@test length(XSim.GLOBAL.G.chr) == 0
	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

	popSize =10
	ngen =5
	popSizeFounder = 20
	sires = Cohort(popSizeFounder);
	dams  = Cohort(popSizeFounder);

	# if LRes is not initialized
	@test_throws ErrorException sire2,dam2=sample_select(sires, dams, popSize, 4, 4, ngen);

	# Once LRes has been initialzed for 1 trait
	A = fill(1.0,1,1)
	XSim.set_residual_var(cholesky(A).U')

	# this happens as the QTL have not been initialized...
	@test_throws DimensionMismatch sire2,dam2=sample_select(sires, dams, popSize, 4, 4, ngen);
end

