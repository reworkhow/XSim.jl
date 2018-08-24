function init(numChr,numLoci,chrLength,geneFreq,
        mapPos,qtl_marker,qtl_effect,mutRate,genotypeErrorRate=0.0,myCommon=common) #assume same chromosomes

    @warn "This function is deprecated. Please use build_genome()."
    #create genome
    locus_array = Array{LocusInfo}(undef,numLoci)
    QTL_index = Array{Int64}(undef,0)
    QTL_effect = Array{Float64}(undef,0)
    chr = Array{ChromosomeInfo}(undef,0)

    for j in 1:numChr
      for i in 1:numLoci
          if mapPos[i]>=chrLength
           error("Map position is not on the chromosome (map position >= chromosome length)")
          end

          locus_index    = numLoci*(j-1)+i
          locus_array[i] = LocusInfo(mapPos[i],[geneFreq[i],1-geneFreq[i]],qtl_marker[locus_index],qtl_effect[locus_index])
          if qtl_marker[locus_index]
            push!(QTL_index,locus_index) #make an array of QTL index for whole Genome
            push!(QTL_effect,locus_index)  #make an array of QTL effects for whole Genome
          end
      end
      chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
      push!(chr,chromosome)
    end
    chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
    chr = fill(chromosome,numChr)
    G = GenomeInfo(chr,numChr,mutRate,genotypeErrorRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array{Animal}(undef,0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
end

function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},mutRate::Float64,genotypeErrorRate=0.0,myCommon=common)
    qtl_marker = fill(false,numLoci)
    qtl_effect = fill(0.0,numLoci)
    init(numChr,numLoci,chrLength,geneFreq, mapPos,qtl_marker,qtl_effect,mutRate,genotypeErrorRate,myCommon)
end


mutable struct XSimMembers
    popSample::Function
    popNew::Function
    popAdd::Function
    popSel::Function
    getGenotypes::Function
    parents::Cohort
    children::Cohort
    gen::Int64
end

function startPop()
    parents  = Cohort(Array(Animal,0),Array(Int64,0,0))
    children = Cohort(Array(Animal,0),Array(Int64,0,0))

    function popSample(numGen::Int64,popSize::Int64)
        popSampleW(numGen,popSize,members)
    end
    function popNew(popSize::Int64)
        popNewW(popSize,members)
    end
    function popAdd(pop::XSimMembers)
    	popAddW(pop,members)
    end
    function popSel(popSize::Int64)
      popSelW(popSize,members)
    end
    function getGenotypes()
        getGenotypesW(members)
    end
    members = XSimMembers(popSample,popNew,popAdd,popSel,getGenotypes,parents,children,0)
    return(members)
end

function getGenotypesW(my)
    getOurGenotypes(my.children)
end

function popSelW(popSize::Int64, ancestors::XSimMembers)
    pheno=get_our_phenotypes(ancestors.children)
    myrank =(length(pheno)+1)-sortperm(pheno)##index from smallest to largest, no functions found
    sel = myrank.<=popSize

    newPop = XSim.startPop()
    newPop.children.animalCohort = ancestors.children.animalCohort[sel]
    newPop.parents  = newPop.children
    newPop.popSample(1,popSize) ##? allow these animals to mate randomly??
    return newPop
end

function popNewW(popSize::Int64,ancestors::XSimMembers)
    newPop = XSim.startPop()
    newPop.parents = ancestors.children
    newPop.popSample(1,popSize)
    return newPop
end

function popAddW(pop::XSimMembers, my::XSimMembers)
	my.parents.animalCohort = [pop.children.animalCohort, my.children.animalCohort]
	my.children.animalCohort = my.parents.animalCohort
end


function popSampleW(numGen::Int64,popSize::Int64,my::XSimMembers)
    if size(my.parents.animalCohort,1)==0
        my.children = sampleFounders(popSize)
        my.parents  = my.children
        numGen -=1
    end
    for i in 1:numGen
        my.gen += 1
        println("Sampling ",popSize," animals into generation: ",my.gen)
        my.children = sampleChildren(my.parents,my.parents,popSize)
        my.parents  = my.children
    end
end

function get_our_phenotypes(my::Cohort)
    QTL_index = common.G.qtl_index #QTL_pos is an array of index of QTLs
    QTL_effects = common.G.qtl_effects
    M = getOurGenotypes(my)
    Q = M[:,QTL_index]
    pheno = Q*QTL_effects #QTL_effects is an array of effects of QTLs
    return pheno
end

function printOurHaps(my::Cohort)
    for i in my.animalCohort
        printMyHaps(i)
    end
end

function printMyHaps(my)
    numberChromosomePair=get_num_chrom(common.G)
    println("I'm an animal with ID( ", my.myID," ) and sire( ", my.sireID, " ) and dam( ", my.damID," )")

    for i in 1:numberChromosomePair
        println(my.genomePat[i].haplotype)
        println(my.genomePat[i].pos)
        println(my.genomeMat[i].haplotype)
        println(my.genomeMat[i].pos)
    end
end

function popCross(popSize::Int64,breed1::XSimMembers,breed2::XSimMembers)
    newPop = XSim.startPop()
    newPop.children = sampleChildren(breed1.children,breed2.children,popSize)
    newPop.parents  = newPop.children
    return newPop
end
