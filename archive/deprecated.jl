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
    newPop.children.animals = ancestors.children.animals[sel]
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
	my.parents.animals = [pop.children.animals, my.children.animals]
	my.children.animals = my.parents.animals
end


function popSampleW(numGen::Int64,popSize::Int64,my::XSimMembers)
    if size(my.parents.animals,1)==0
        my.children = Cohort(popSize)
        my.parents  = my.children
        numGen -=1
    end
    for i in 1:numGen
        my.gen += 1
        println("Sampling ",popSize," animals into generation: ",my.gen)
        my.children = get_children(my.parents,my.parents,popSize)
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
    for i in my.animals
        printMyHaps(i)
    end
end

function printMyHaps(my)
    numberChromosomePair=get_num_chrom(common.G)
    println("I'm an animal with ID( ", my.myID," ) and sire( ", my.sireID, " ) and dam( ", my.damID," )")

    for i in 1:numberChromosomePair
        println(my.genome_sire[i].haplotype)
        println(my.genome_sire[i].pos)
        println(my.genome_dam[i].haplotype)
        println(my.genome_dam[i].pos)
    end
end

function popCross(popSize::Int64,breed1::XSimMembers,breed2::XSimMembers)
    newPop = XSim.startPop()
    newPop.children = get_children(breed1.children,breed2.children,popSize)
    newPop.parents  = newPop.children
    return newPop
end


function getOneHapsDeprecated(genome::Array{Chromosome,1})
    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair
        numLoci=common.G.chr[i].numLoci
        resize!(genome[i].haplotype,numLoci)

        numOri=length(genome[i].ori)
        push!(genome[i].pos,common.G.chr[i].chrLength) #this may make it not efficient

        iLoci   = 1
        position = common.G.chr[i].mapPos[iLoci]
        #println("getOneHaps(): number of segments in chromosome ",i,": ",numOri)
        chrom = common.G.chr[i]
        for segment in 1:numOri
            place = 0
            whichFounder=ceil(Integer,genome[i].ori[segment]/2)
            genome_sireorMatInThisFounder=(genome[i].ori[segment]%2==0) ? common.founders[whichFounder].genome_dam[i] : common.founders[whichFounder].genome_sire[i]

            startPos = genome[i].pos[segment]
            endPos   = genome[i].pos[segment+1]
            segLen   = 0
            endLocus = 0
            if segment < numOri
                if chrom.mapPos[numLoci] > endPos #there is at least a single locus in the next segment
                    while position >= startPos && position<endPos
                        iLoci = iLoci+1
                        segLen = segLen+1
                        position = chrom.mapPos[iLoci]
                    end
                    endLocus = iLoci-1
                else
                    break
                end

            elseif iLoci < numLoci
                segLen = numLoci - iLoci + 1
                endLocus = numLoci
            end

            genome[i].haplotype[(endLocus-segLen+1):endLocus]=genome_sireorMatInThisFounder.haplotype[(endLocus-segLen+1):endLocus]
        end

        for j in 1:length(genome[i].mut)
            whichlocus = findfirst(common.G.chr[i].mapPos .== genome[i].mut[j])
            genome[i].haplotype[whichlocus] = 1 - genome[i].haplotype[whichlocus]
        end

        pop!(genome[i].pos) #remove temporary added chrLength
    end
end

        
function getOneHapsBinSearchMap(genome::Array{Chromosome,1})
    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair
        numLoci=common.G.chr[i].numLoci
        resize!(genome[i].haplotype,numLoci)

        numOri=length(genome[i].ori)
        push!(genome[i].pos,common.G.chr[i].chrLength) #this may make it not efficient

        iLocus   = 1
        position = common.G.chr[i].mapPos[iLocus]
        println("getOneHaps(): number of segments in chromosome ",i,": ",numOri)
        chrom = common.G.chr[i]
        lociPerM = round(Int64, numLoci/genome[i].pos[numOri+1])
        segLen = 0
        prevSegLen = 0
        endLocus = 0

        for segment in 1:numOri
            println("getOneHaps(): segment=",segment)
            flush(stdout)
            whichFounder=ceil(Integer,genome[i].ori[segment]/2)
            genome_sireorMatInThisFounder=(genome[i].ori[segment]%2==0) ? common.founders[whichFounder].genome_dam[i] : common.founders[whichFounder].genome_sire[i]

            startPos = genome[i].pos[segment]
            endPos   = genome[i].pos[segment+1]
            prevSegLen += segLen       
            segLen   = 0
            if segment < numOri 
                println("getOneHaps(): startPos=",startPos, " endPos=",endPos," lociPerM=",lociPerM)
                flush(stdout)
                numLociUntilGuessedPos = round(Int64, endPos * lociPerM)
                if numLociUntilGuessedPos > numLoci
                    numLociUntilGuessedPos = numLoci
                end
                guessedPos = chrom.mapPos[numLociUntilGuessedPos]
                prevGuessedPos = startPos
                prevNumLociUntilGuessedPos = 0
                #println("getOneHaps(): guessedPos=",guessedPos," prevGuessedPos=",prevGuessedPos," numLociUntilGuessedPos=",numLociUntilGuessedPos)
                topDown = 0
                bottomUp = 0 
                stepSize = 100
                prevStepSize = 0
                while stepSize > 25 && (stepSize != prevStepSize) && (bottomUp<3 || topDown<3)
                    if guessedPos < endPos
                        bottomUp += 1       
                        prevGuessedPos = guessedPos        
                        guessedPos = (endPos-guessedPos)/2
                        prevNumLociUntilGuessedPos = numLociUntilGuessedPos
                        numLociUntilGuessedPos += ceil(Integer, guessedPos * lociPerM)
                        #println("getOneHaps1(): guessedPos=",guessedPos," prevGuessedPos=",prevGuessedPos," numLociUntilGuessedPos=",numLociUntilGuessedPos," prevNumLociUntilGuessedPos=",prevNumLociUntilGuessedPos," stepSize=",stepSize," prevStepSize=",prevStepSize," bottomUp=",bottomUp)
                        if numLociUntilGuessedPos > numLoci
                            numLociUntilGuessedPos = numLoci
                        end
                        prevStepSize = stepSize
                        stepSize = abs(numLociUntilGuessedPos - prevNumLociUntilGuessedPos)
                        guessedPos = chrom.mapPos[numLociUntilGuessedPos]
                        #println("getOneHaps1(): guessedPos=",guessedPos," prevGuessedPos=",prevGuessedPos," numLociUntilGuessedPos=",numLociUntilGuessedPos," prevNumLociUntilGuessedPos=",prevNumLociUntilGuessedPos," stepSize=",stepSize," prevStepSize=",prevStepSize," bottomUp=",bottomUp)
                    else
                        topDown += 1       
                        prevGuessedPos = guessedPos        
                        guessedPos = (guessedPos-endPos)/2
                        prevNumLociUntilGuessedPos = numLociUntilGuessedPos
                        numLociUntilGuessedPos -= floor(Integer, guessedPos * lociPerM)
                        #println("getOneHaps2(): guessedPos=",guessedPos," prevGuessedPos=",prevGuessedPos," numLociUntilGuessedPos=",numLociUntilGuessedPos," prevNumLociUntilGuessedPos=",prevNumLociUntilGuessedPos," stepSize=",stepSize," prevStepSize=",prevStepSize," topDown=",topDown)
                        prevStepSize = stepSize
                        stepSize = abs(numLociUntilGuessedPos - prevNumLociUntilGuessedPos)
                        guessedPos = chrom.mapPos[numLociUntilGuessedPos]
                        #println("getOneHaps2(): guessedPos=",guessedPos," prevGuessedPos=",prevGuessedPos," numLociUntilGuessedPos=",numLociUntilGuessedPos," prevNumLociUntilGuessedPos=",prevNumLociUntilGuessedPos," stepSize=",stepSize," prevStepSize=",prevStepSize," topDown=",topDown)
                    end
                end
                position = guessedPos
                iLocus = numLociUntilGuessedPos
                #println("getOneHaps(): position=",position)
                #println("getOneHaps(): iLocus=",iLocus)
                           
                println("Diff between end of segment and pos after binary search = ",round(endPos - position,digits=6), " current locus=",iLocus)
                println("Number of top down searches=",topDown)
                println("Number of bottom up searches=",bottomUp)

                            
                startLoc = iLocus
                if startLoc < numLoci
                    if position < endPos
                        while position < endPos
                            iLocus += 1
                            position = chrom.mapPos[iLocus]
                            #println("getOneHaps3(): iLocus=",iLocus," position=",position)
                        end 
                        endLocus = iLocus - 1
                        segLen = endLocus - prevSegLen
                        #println("getOneHaps3(): endLocus=",endLocus)
                    elseif position > endPos
                        while position > endPos
                            iLocus -= 1
                            position = chrom.mapPos[iLocus]
                            #println("getOneHaps4(): iLocus=",iLocus," position=",position)
                        end 
                        endLocus = iLocus
                        segLen = endLocus - prevSegLen
                        #println("getOneHaps4(): endLocus=",endLocus)
                    else
                        endLocus = iLocus
                        segLen = endLocus - prevSegLen
                        #println("getOneHaps5(): endLocus=",endLocus)
                    end                    
                    #println("getOneHaps6(): endLocus=",endLocus)
                    println("getOneHaps(): linear search finished after ",endLocus - startLoc," loci.")
                else
                    segLen = numLoci - endLocus
                    endLocus = numLoci      
                end
            elseif iLocus <= numLoci
                segLen = numLoci - endLocus
                endLocus = numLoci
                #println("getOneHaps(): endLocus=",endLocus)
            end  
            println("getOneHapsNew(): segment=",segment," endLocus=",endLocus, " seglen=", segLen, " prevSeglen=", prevSegLen)
            if segLen > 0 
               genome[i].haplotype[(endLocus-segLen+1):endLocus]=genome_sireorMatInThisFounder.haplotype[(endLocus-segLen+1):endLocus]
            end
        end
                        
                        
        for j in 1:length(genome[i].mut)
            whichlocus = findfirst(common.G.chr[i].mapPos .== genome[i].mut[j])
            genome[i].haplotype[whichlocus] = 1 - genome[i].haplotype[whichlocus]
        end

        pop!(genome[i].pos) #remove temporary added chrLength
    end
end

