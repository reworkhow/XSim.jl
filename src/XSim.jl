module XSim

using Distributions

tempPos=Array(Float64,100000)
tempOri=Array(Int64,  100000)

# Types and Methods for Information on the Genome
type LocusInfo
    map_pos::Float64
    allele_freq::Array
    QTL::Bool
    QTL_effect::Float64
end

function get_num_alleles(my)
    return length(my.allele_freq)
end

type ChromosomeInfo
    chrLength::Float64
    numLoci::Int64
    mapPos::Array{Float64,1}
    loci::Array{LocusInfo,1}
end

function set_num_loci(my,n::Int64)
    my.numLoci=n; nothing
end

function get_num_loci(my)
    return my.numLoci
end

function mk_mappos_from_locus_info(my)
    my.mapPos.resize(my.chrLength)
    for i in 1:my.chrLength
        my.mapPos[i]=my.locus[i].map_pos
    end
end

type GenomeInfo
    chr::Array{ChromosomeInfo,1}
    numChrom::Int64
    mutRate::Float64
    qtl_index::Array{Int64,1}
    qtl_effects::Array{Float64,1}
end

function set_num_chrom(my,n::Int64)
    my.numChrom=n; nothing
end

function get_num_chrom(my)
    return my.numChrom
end

function make_map_pos(my)
    for i in 1:my.numChrom
        mk_mappos_from_locus_info(my.chr[i])
    end
    mapPosDone = true
end

#Types and Methods for Simulating Genotypes of Animals

type Chromosome
    haplotype::Array{Int64,1}
    ori::Array{Int64,1}
    pos::Array{Float64,1}
end

type Animal
    genomePat::Array{Chromosome,1}
    genomeMat::Array{Chromosome,1}
    breedComp::Array{Float64,1}        ##################### 2015.6.4
    myID::Int64
    sireID::Int64
    damID::Int64
    phenVal::Float64
    genVal::Float64
end

function Animal(mySire::Int64, myDam::Int64)
    my=Animal(Array(Chromosome,common.G.numChrom),Array(Chromosome,common.G.numChrom), Array(Float64,0),0,0,0,
              -9999.0,-9999.0)
    my.sireID = mySire
    my.damID  = myDam
    my.myID   = common.countId
    common.countId += 1
    for i in 1:common.G.numChrom
        my.genomePat[i]=Chromosome(Array(Int64,0),Array(Int64,0),Array(Float64,0))
        my.genomeMat[i]=Chromosome(Array(Int64,0),Array(Int64,0),Array(Float64,0))
    end
    return my
end

type Cohort
    animalCohort::Array{Animal,1}
    npMatrix::Array{Int64,2}
end

# Type for Storing Globals

type CommonToAnimals
    founders::Array{Animal,1}
    G::GenomeInfo
    countChromosome::Int64
    countId::Int64
end

function sampleFounder()
    my = Animal(0,0)
    initFounderPosOri(my)
    initFounderHaps(my)
    return(my)
end

function initFounderHaps(my)

    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair

        numLoci=common.G.chr[i].numLoci

        Base.resize!(my.genomePat[i].haplotype,numLoci)
        Base.resize!(my.genomeMat[i].haplotype,numLoci)

        for j in 1:numLoci
            p=common.G.chr[i].loci[j].allele_freq[1]
            my.genomePat[i].haplotype[j]=rand(Bernoulli(p))
            my.genomeMat[i].haplotype[j]=rand(Bernoulli(p))
        end
    end
end

function initFounderPosOri(my)
        numberChromosomePair=get_num_chrom(common.G)
        for i in 1:numberChromosomePair
            my.genomePat[i].ori=[common.countChromosome]
            my.genomePat[i].pos=[0.0]

            my.genomeMat[i].ori=[common.countChromosome+1]
            my.genomeMat[i].pos=[0.0]
        end
        common.countChromosome += 2
end

function sampleFounders(numAnimals::Int64)
        my=Cohort(Array(Animal,0),Array(Int64,0,0))

        println("Sampling ",numAnimals," animals into base population.")
        resize!(my.animalCohort,numAnimals)
        for i in 1:numAnimals
            animal=sampleFounder()
            my.animalCohort[i] = animal
            push!(common.founders,animal)
        end
        return(my)
end

function getRandomInd(my::Cohort)
    cohortSize=length(my.animalCohort)
    thisOne=rand(1:cohortSize)
    return(my.animalCohort[thisOne])
end

function sampleChildren(fathers::Cohort,mothers::Cohort,numAnimals::Int64)
    my=Cohort(Array(Animal,0),Array(Int64,0,0))

    #println("Sampling ",numAnimals," animals into next generation.")
    resize!(my.animalCohort,numAnimals)

    for i in 1:numAnimals
        animal=sampleNonFounder(getRandomInd(fathers),getRandomInd(mothers))
        #println("Sampled animal: ",animal.myID)
        my.animalCohort[i]=animal
    end

    #print("The length of pos vector is ------------->")
    #printMyHaps(my.animalCohort[numAnimals])

    return(my)
end

function sampleNonFounder(father,mother)
        my=Animal(father.myID,mother.myID)

        numberChromosomePair=get_num_chrom(common.G)
        resize!(my.genomePat,numberChromosomePair)
        resize!(my.genomeMat,numberChromosomePair)

        sampleMyPosOri(my,father,mother)
        return(my)
end

function sampleMyPosOri(my,father,mother)
    sampleOnePosOri(my.genomePat,father)
    sampleOnePosOri(my.genomeMat,mother)
end

function sampleOnePosOri(genome::Array{Chromosome,1},parent::Animal)
    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair

        genome[i]=Chromosome(Array(Int64,0),Array(Int64,1),Array(Float64,1))

        currentChrom=(rand(Bernoulli(0.5))==1)?parent.genomePat[i]:parent.genomeMat[i]

        chrLength=common.G.chr[i].chrLength

        binomialN=convert(Int64,ceil(chrLength*3+1))
        numCrossover=rand(Binomial(binomialN,chrLength/binomialN))
        rec=[0.0]

        for irec in 1:numCrossover
            push!(rec,chrLength*rand(1)[1])
        end

        push!(rec,chrLength)
        sort!(rec)                #rec is like 0.00,0.23,0.45,0.76,1.00

        numTemp=1

        for j in 1:(length(rec)-1)

            for k in 1:length(currentChrom.pos)
              if currentChrom.pos[k] >= rec[j] && currentChrom.pos[k] < rec[j+1]
                tempPos[numTemp]=currentChrom.pos[k]
                tempOri[numTemp]=currentChrom.ori[k]
                numTemp=numTemp+1
              elseif currentChrom.pos[k]>=rec[j+1]
                break
              end
            end

            currentChrom=(currentChrom==parent.genomeMat[i])?parent.genomePat[i]:parent.genomeMat[i]

            findRecOri=0
            m=1
            lengthChrPos = length(currentChrom.pos)
            while m<=lengthChrPos && currentChrom.pos[m]<=rec[j+1] #pos[m] cannot be the length of chromosome
               m += 1
               findRecOri += 1
            end
            tempPos[numTemp]=rec[j+1]
            tempOri[numTemp]=currentChrom.ori[findRecOri]
            numTemp=numTemp+1 #**
        end

        numTemp=numTemp-1 #remove the last one got from **
        numTemp=numTemp-1 #remove the last one which is the length of chromsome

        resize!(genome[i].pos,numTemp)
        resize!(genome[i].ori,numTemp)

        ##merging of nearby same ori;not sure whether more effcient

        genome[i].pos[1]=tempPos[1]
        genome[i].ori[1]=tempOri[1]
        this=1
        for(m in 2:numTemp)
            if(tempOri[m]!=tempOri[m-1])
                this=this+1
                genome[i].pos[this]=tempPos[m]
                genome[i].ori[this]=tempOri[m]
            end
        end

        resize!(genome[i].pos,this) #resize after merging
        resize!(genome[i].ori,this)

        #println(genome[i].pos)
        #println(genome[i].ori)

    end
end

function getMyHaps(my)
    getOneHaps(my.genomePat)
    getOneHaps(my.genomeMat)
end

function getOneHaps(genome::Array{Chromosome,1})
    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair
        numLoci=common.G.chr[i].numLoci
        resize!(genome[i].haplotype,numLoci)

        numOri=length(genome[i].ori)
        push!(genome[i].pos,common.G.chr[i].chrLength) #this may make it not efficient

        iLoci   = 1
        position = common.G.chr[i].mapPos[iLoci]

        for segment in 1:numOri
            whichFounder=iceil(genome[i].ori[segment]/2)
            genomePatorMatInThisFounder=(genome[i].ori[segment]%2==0)?common.founders[whichFounder].genomeMat[i]:common.founders[whichFounder].genomePat[i]

            startPos = genome[i].pos[segment]
            endPos   = genome[i].pos[segment+1]
            segLen   = 0

            while position >= startPos && position<endPos
                iLoci = iLoci+1
                segLen = segLen+1

                if iLoci>numLoci
                    break
                end

                position = common.G.chr[i].mapPos[iLoci]
            end

            endLoci=iLoci-1
            genome[i].haplotype[(endLoci-segLen+1):endLoci]=genomePatorMatInThisFounder.haplotype[(endLoci-segLen+1):endLoci]
        end

        pop!(genome[i].pos)
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

function outputPedigree(my::Cohort, fileName::String)
    pedText  = fileName * ".ped"
    genText  = fileName * ".gen"
    brcText  = fileName * ".brc"
    pedStream = open(pedText,"a")
    genStream = open(genText,"a")
    brcStream = open(brcText,"a")
    for animal in my.animalCohort
        getMyHaps(animal)
        genotypes=getMyGenotype(animal)
        @printf(pedStream,  "%9d %9d %9d \n", animal.myID, animal.sireID, animal.damID)
        for j=1:length(animal.breedComp)
            @printf(brcStream, "%5.3f ", animal.breedComp[j])
        end
        @printf(brcStream, "\n")
        @printf(genStream, "%19d", animal.myID)
        for j=1:length(genotypes)
            @printf(genStream, "%3d", genotypes[j])
        end
        @printf(genStream, "\n")
    end
    close(brcStream)
    close(pedStream)
    close(genStream)

end

function outputPedigree(this::Cohort, fileName::String, sel::Array{Int64,1})
    cohort = cohortSubset(this,sel)
    outputPedigree(this, fileName)
end

function getOurGenVals(my::Cohort)
    n = size(my.animalCohort,1)
    genVals = Array(Float64,n)
    for (i,animal) = enumerate(my.animalCohort)
        #println(i)
        if animal.genVal==-9999
            getMyHaps(animal)
            myGenotypes = getMyGenotype(animal)
            animal.genVal = dot(myGenotypes[common.G.qtl_index],common.G.qtl_effects)
        end
        genVals[i] = animal.genVal
    end
    return genVals
end

function getOurPhenVals(my::Cohort, varRes)
    stdRes = sqrt(varRes)
    n = size(my.animalCohort,1)
    phenVals = Array(Float64,n)
    genVals  = getOurGenVals(my)
    for (i,animal) = enumerate(my.animalCohort)
        if animal.phenVal==-9999
            animal.phenVal =  animal.genVal + (randn(1)*stdRes)[1]
        end
        phenVals[i] = animal.phenVal
    end
    return phenVals
end



function getOurHaps(my::Cohort)
    for i in my.animalCohort
        getMyHaps(i)
    end
end

function printOurHaps(my::Cohort)
    for i in my.animalCohort
        printMyHaps(i)
    end
end

function printOurHaps(my::Cohort)
    for i in my.animalCohort
        printMyHaps(i)
    end
end

function getOurGenotypes(my::Cohort)
    getOurHaps(my)
    npMatrix=Array(Int64,length(my.animalCohort), common.G.numChrom*common.G.chr[1].numLoci)
    for (i,value) in enumerate(my.animalCohort)
        npMatrix[i,:]=getMyGenotype(value)
    end
    return npMatrix
end

function cohortSubset(my::Cohort,sel::Array{Int64,1})
    animals = Array(Animal,size(sel,1))
    for (i,j) = enumerate(sel)
        animals[i] = my.animalCohort[j]
    end
    return Cohort(animals,Array(Int64,0,0))
end

function getOurGenotypes(my::Cohort,sel::Array{Int64,1})
    cohort = cohortSubset(my,sel)
    return getOurGenotypes(cohort)
end

function getMyGenotype(my)
    myGenotype=Array(Int64,0)
    for i in 1:common.G.numChrom
        append!(myGenotype, my.genomePat[i].haplotype+my.genomeMat[i].haplotype)
    end
    #println(myGenotype)
    return myGenotype
end

#####
function get_our_phenotypes(my::Cohort)
    QTL_index = common.G.qtl_index #QTL_pos is an array of index of QTLs
    QTL_effects = common.G.qtl_effects
    M = getOurGenotypes(my)
    Q = M[:,QTL_index]
    pheno = Q*QTL_effects #QTL_effects is an array of effects of QTLs
    return pheno
end

# Make object for storing globals
G = GenomeInfo(Array(ChromosomeInfo,0),0,0.0,[],[])
common = CommonToAnimals(Array(Animal,0),G,0,0)


##encapsulation
type XSimMembers
    popSample::Function
    popNew::Function
    popAdd::Function
    popSel::Function
    getGenotypes::Function
    parents::Cohort
    children::Cohort
    gen::Int64
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

function popCross(popSize::Int64,breed1::XSimMembers,breed2::XSimMembers)
    newPop = XSim.startPop()
    newPop.children = sampleChildren(breed1.children,breed2.children,popSize)
    newPop.parents  = newPop.children
    return newPop
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

function getGenotypesW(my)
    getOurGenotypes(my.children)
end

function init(numChr::Int64,numLoci::Int64,chrLength::Float64,geneFreq::Array{Float64,1},
        mapPos::Array{Float64,1},qtl_marker::Array{Bool,1},qtl_effect::Array{Float64,1},mutRate::Float64,myCommon=common)

    #create genome
    locus_array = Array(LocusInfo,numLoci)
    QTL_index = Array(Int64,0)
    QTL_effect = Array(Float64,0)
    chr = Array(ChromosomeInfo,0)

    for j in 1:numChr
      for i in 1:numLoci
          locus_array[i] = LocusInfo(mapPos[i],[geneFreq[i],1-geneFreq[i]],qtl_marker[i],qtl_effect[i])
          if qtl_marker[i]
            push!(QTL_index,numLoci*(j-1)+i) #make an array of QTL index for whole Genome
            push!(QTL_effect,qtl_effect[i])  #make an array of QTL effects for whole Genome
          end
      end
      chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
      push!(chr,chromosome)
    end
    chromosome = ChromosomeInfo(chrLength,numLoci,mapPos,locus_array)
    chr = fill(chromosome,numChr)
    G = GenomeInfo(chr,numChr,mutRate,QTL_index,QTL_effect)

    # Init common
    myCommon.founders=Array(Animal,0)
    myCommon.G = G
    myCommon.countId = 1
    myCommon.countChromosome = 1; nothing
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

export startPop
export popCross

function concatCohorts(cohortLst...)
    # returns a cohort with concatenation of the animalCohorts from the arguments
    res = Cohort(Array{Animal,1},Array{Int64,2})
    for i in cohortLst
        res.animalCohort = [res.animalCohort, i.animalCohort]
    end
    return res
end

# Types and methods for sampling individuals from a given pedigree

type PedNode
    ind::Int64
    sire::Int64
    dam::Int64
end

function samplePed(ped::Array{PedNode,1})
    animals = Array(Animal,size(ped,1))
    for i in ped
        if i.ind <= i.sire || i.ind <= i.dam
            throw(Exception("ind < sire or dam \n"))
        end
        if i.sire == 0
            animal = sampleFounder()
            animals[i.ind] = animal
            push!(common.founders,animal)
        else
            animals[i.ind] = sampleNonFounder(animals[i.sire],animals[i.dam])
        end
    end
    res = Cohort(animals,Array(Int64,0,0))
end

function sampleSel(popSize, nSires, nDams, nGen, varRes)
    maleCandidates   = sampleFounders(int(popSize/2))
    femaleCandidates = sampleFounders(int(popSize/2))
    return sampleSel(popSize, nSires, nDams, nGen,maleCandidates,femaleCandidates, varRes)
end

function sampleSel(popSize, nSires, nDams, nGen,maleCandidates,femaleCandidates,varRes;gen=1,fileName="")
    sires = Cohort(Array(Animal,0),Array(Int64,0,0))
    dams  = Cohort(Array(Animal,0),Array(Int64,0,0))
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    for i=1:nGen
        @printf "Generation %5d: sampling %5d animals\n" gen+i popSize
        y = getOurPhenVals(maleCandidates,varRes)
        sires.animalCohort = maleCandidates.animalCohort[sortperm(y)][(end-nSires+1):end]
        y = getOurPhenVals(femaleCandidates,varRes)
        dams.animalCohort = femaleCandidates.animalCohort[sortperm(y)][(end-nDams+1):end]
        boys = sampleChildren(sires,dams,int(popSize/2))
        gals = sampleChildren(sires,dams,int(popSize/2))
        if fileName!=""
            outputPedigree(boys,fileName)
            outputPedigree(gals,fileName)
        end
        maleCandidates.animalCohort   = [sires.animalCohort, boys.animalCohort]
        femaleCandidates.animalCohort = [dams.animalCohort,  gals.animalCohort]
    end
    gen += nGen
    return boys,gals, gen
end

function sampleRan(popSize, nGen,sires,dams;gen=1)
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    for i=1:nGen
        @printf "Generation %5d: sampling %5d animals\n" gen+i popSize
        boys = sampleChildren(sires,dams,int(popSize/2))
        gals = sampleChildren(sires,dams,int(popSize/2))
        sires = boys
        dams  = gals
    end
    gen += nGen
    return boys,gals, gen
end

end # module
