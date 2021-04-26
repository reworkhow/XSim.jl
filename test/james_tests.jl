using Test
using XSim

function random_mate()
    clearGlobals()
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
	# sires1,dams1,gen1 = sampleRan(popSize, ngen, sires, dams);
    G0 = Array{Float64,2}(undef,1,1)
	G0[1,1] = 1
    XSim.setResidualVariance(G0)

    sires1,dams1,gen1 = selection_for_ngenerations(popSize, 1, 2, sires, dams;
                                    ngenerations=ngen, strategy="random")



	@test typeof(sires1)==XSim.Cohort
	@test typeof(dams1)==XSim.Cohort
end

# include("XSim/src/XSim.jl")
# using .XSim
# chrLength = 0.1  #length of each chromosome
# numChr    = 2    #number of chromosomes
# nmarkers  = 10   #number of loci for each chromosome
# nQTL      = 1    #number of QTL for each chromosomefects,mutRate);
# build_genome(numChr, chrLength, nmarkers, nQTL)

# popSizeFounder = 5
# sires = Cohort(popSizeFounder);
# dams  = Cohort(popSizeFounder);

# boys, girls = sample_select(sires, dams, 20, 3, 3, 5)
# b2, g2 = sample_random(sires, dams, 20, 3, 3, 5)



# G0 = Any[]
# Sampling 2 animals into base population.
# Sampling 2 animals into base population.
# Generation     2: sampling     5 males and     5 females
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 2
# Number of mothers used: 2
# Number of offspring generated: 5
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 2
# Number of mothers used: 2
# Number of offspring generated: 5
# Generation     3: sampling     5 males and     5 females
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 3
# Number of mothers used: 5
# Number of offspring generated: 5
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 2
# Number of mothers used: 3
# Number of offspring generated: 5
# Generation     4: sampling     5 males and     5 females
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 3
# Number of mothers used: 3
# Number of offspring generated: 5
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 3
# Number of mothers used: 3
# Number of offspring generated: 5
# Generation     5: sampling     5 males and     5 females
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 5
# Number of mothers used: 3
# Number of offspring generated: 5
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 2
# Number of mothers used: 4
# Number of offspring generated: 5
# Generation     6: sampling     5 males and     5 females
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 3
# Number of mothers used: 4
# Number of offspring generated: 5
# RUNNING OLD FUNCTIONNNNNNumber of fathers used: 2
# Number of mothers used: 4
# Number of offspring generated: 5


# # example to get trait in cohort
# mutable struct objA
#     X::Int64
#     Y::Int64
#     objA(x, y)=new(x, y)
# end

# mutable struct objAnimal
#     obj::Array{objA, 1}
#     objAnimal()=new()
# end

# mutable struct objCohort
#     obj::Array{objAnimal, 1}
# end

# ani = objAnimal(Array{objA}(undef, 5))
# for i in 1:5
#     ani.obj[i] = objA(1, 2)
# end

# out = (obj->[obj.X, obj.Y]).(ani.obj)
# hcat(out...)'


# anj = objAnimal(Array{objA}(undef, 5))
# for i in 1:5
#     anj.obj[i] = objA(3, 4)
# end

# function gx(obj::objAnimal)
#     return (x -> x.X).(obj.obj)
# end


# co = objCohort([ani, anj])

# out = (animal->gx(animal)).(co.obj)
# outt = hcat(out...)'

# # Test animal
# mutable struct Animal
#     ID          ::Int64
#     sire        ::Animal
#     dam         ::Animal
#     genome_sire ::Array{Float16, 1}
#     n_traits    ::Int
#     breedComp   ::Array{Float64   , 1}
#     is_founder  ::Bool

#     # Constructor: allow undef animal when it's a founder's parents
#     Animal() = new(0)
#     function Animal(sire      ::Animal,
#                     dam       ::Animal;
#                     is_founder::Bool=false)
#         animal = new(GLOBAL.countId, sire, dam,
#                      Array{Chromosome}(undef, 3),
#                      0,
#                      Array{Float64   }(undef, 0),
#                      is_founder)
#         GLOBAL.countId += 1
#         set_genome(animal)

#         if is_founder
#             push!(GLOBAL.founders, animal)
#         end

#         return animal
#     end

#     function Animal(is_founder::Bool)
#         # instantiate a founder
#         return Animal(Animal(), Animal(), is_founder=is_founder)
#     end

#     function set_genome(animal::Animal)
#         is_founder = animal.is_founder
#         for i in 1:GLOBAL.G.numChrom
#             animal.genome_sire[i] = Chromosome(i, is_founder=is_founder)
#             animal.genome_dam[i]  = Chromosome(i, is_founder=is_founder)
#         end
#     end
# end






# function t1()
# 	clearGlobals()
# 	@test XSim.GLOBAL.G.numChrom == 0
# 	@test length(XSim.GLOBAL.G.chr) == 0
# 	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

# 	chrLength= 0.1  #length of each chromosome
# 	numChr   = 2    #number of chromosomes
# 	nmarkers = 10   #number of loci for each chromosome
# 	nQTL     = 1    #number of QTL for each chromosomefects,mutRate);
# 	build_genome(numChr,chrLength,nmarkers,nQTL)
# 	@test XSim.GLOBAL.G.numChrom == numChr
# 	@test length(XSim.GLOBAL.G.chr) == numChr
# 	@test XSim.GLOBAL.countChromosome== 1
# end

# function t2()
#     clearGlobals()
# 	@test XSim.GLOBAL.G.numChrom == 0
# 	@test length(XSim.GLOBAL.G.chr) == 0
# 	@test XSim.GLOBAL.countChromosome == 0 #Number of chromosomes in founder population ; starts at 1

# 	numChr = 3
# 	chrLength = [1.0, 1.1, 0.9]
# 	numLoci = [10 ,15 ,20]
# 	nTraits = 2
# 	geneFreq   = Array{Array{Float64,1},1}(undef,0)
# 	qtlIndex  = Array{Array{Int64,1},1}(undef,0)
# 	qtlEffects = Array{Array{Float64,2},1}(undef,0)
# 	numQTLOnChr =[2, 3, 1]
# 	numQTL=sum(numQTLOnChr)
# 	qtlIndex = [[1; 4],[3;7;9],[5]]
# 	geneFreq = [[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
# 				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5],
# 				[.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5;.5]]
# 	mPos = [[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5],
# 			[.05;.1;.15;.2;.25;.3;.35;.40;.45;.5;.55;.65;.75;.85;.95],
# 			[.05;.1;.15;.2;.25;.3;.35;.4;.45;.5;.54;.58;.65;.75;.758;.759;.8;.85;.851;.859]]
# 	push!(qtlEffects,randn(2,nTraits))
# 	push!(qtlEffects,randn(3,nTraits))
# 	push!(qtlEffects,randn(1,nTraits))
# 	G0 = [1. 0.5;0.5 2.]
# 	XSim.build_genome(numChr,chrLength,numLoci,geneFreq,mPos,qtlIndex,qtlEffects,nTraits,G0)
# 	XSim.setResidualVariance(G0)

# 	popSize = 10
# 	ngen = 3
# 	popSizeFounder = 20
# 	sires = sampleFounders(popSizeFounder)
# 	dams  = sampleFounders(popSizeFounder)

# 	sire1,dam1,gen1=sampleSel(popSize, 4, 4, ngen,sires, dams);
# end