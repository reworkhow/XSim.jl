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

function t2()
    clearGlobals()
	@test XSim.common.G.numChrom == 0
	@test length(XSim.common.G.chr) == 0
	@test XSim.common.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

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
	XSim.setResidualVariance(G0)

	popSize = 10
	ngen = 3
	popSizeFounder = 20
	sires = sampleFounders(popSizeFounder)
	dams  = sampleFounders(popSizeFounder)

	sire1,dam1,gen1=sampleSel(popSize, 4, 4, ngen,sires, dams);
end