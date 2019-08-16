include("genotype.jl")
include("phenotype.jl")
include("moreInfo.jl")

function outputPedigree(this::Cohort, fileName::AbstractString, sel::Array{Int64,1})
    cohortSel = cohortSubset(this,sel)
    outputPedigree(cohortSel, fileName)
end

function outputPedigree(my::Cohort, fileName::AbstractString)
    pedText  = fileName * ".ped"
    genText  = fileName * ".gen"
    brcText  = fileName * ".brc"
    pheText  = fileName * ".phe"
    pedStream = open(pedText,"a")
    genStream = open(genText,"a")
    brcStream = open(brcText,"a")
    pheStream = open(pheText,"a")
    getOurPhenVals(my)
    nTraits = size(common.LRes,2)
    for (k,animal) = enumerate(my.animalCohort)
        if k%1000 == 0
            println("outputPedigree(): ", k)
        end
        genotypes=getMyGenotype(animal)
        @printf(pedStream,  "%19d %19d %19d \n", animal.myID, animal.sireID, animal.damID)
        @printf(brcStream,  "%19d ", animal.myID)
        for j=1:length(animal.breedComp)
            @printf(brcStream, "%5.3f ", animal.breedComp[j])
        end
        @printf(brcStream, "\n")
        @printf(pheStream, "%19d ", animal.myID)
        for i=1:nTraits
            @printf(pheStream, "%11.3f ",animal.phenVal[i])
        end
        for i=1:nTraits
            @printf(pheStream, "%11.3f ",animal.genVal[i])
        end
        @printf(pheStream, "\n")
        @printf(genStream, "%19d", animal.myID)
        for j=1:length(genotypes)
            @printf(genStream, "%2d", genotypes[j])
        end
        @printf(genStream, "\n")
    end
    close(pheStream)
    close(brcStream)
    close(pedStream)
    close(genStream)
end

function outputHapData(this::Cohort, fileName::AbstractString)
    getOurHaps(this)
    hapStream = open(fileName,"w")
    @printf(hapStream, "%10s ","animalID")
    for chr in 1:common.G.numChrom
        for locus in 1:common.G.chr[chr].numLoci
            @printf(hapStream, "%2s", "p")
            @printf(hapStream, "%2s", "m")
        end
    end
    @printf(hapStream, "\n")
    for animal in this.animalCohort
        @printf(hapStream, "%10d ", animal.myID)
        for chr in 1:common.G.numChrom
            for locus in 1:common.G.chr[chr].numLoci
                @printf(hapStream, "%2d", animal.genomePat[chr].haplotype[locus])
                @printf(hapStream, "%2d", animal.genomeMat[chr].haplotype[locus])
            end
        end
        @printf(hapStream, "\n")
    end
    close(hapStream)
end

function outputGenData(this::Cohort, fileName::AbstractString )
    getOurHaps(this)
    genStream = open(fileName,"w")
    @printf(genStream, "%10s ","animalID")
    snpID = 1
    for chr in 1:common.G.numChrom
        for locus in 1:common.G.chr[chr].numLoci
            @printf(genStream, "%8d ", snpID)
            snpID += 1
        end
    end
    @printf(genStream, "\n")
    for animal in this.animalCohort
        @printf(genStream, "%10d ", animal.myID)
        for chr in 1:common.G.numChrom
            genotypes = animal.genomePat[chr].haplotype + animal.genomeMat[chr].haplotype
            if common.G.genotypeErrorRate > 0.0
                numberErrors = int(common.G.genotypeErrorRate * common.G.chr[chr].numLoci)
                errorPos = sample([1:common.G.chr[chr].numLoci],numberErrors)
                for i in errorPos
                    if genotypes[i]==2
                        genotypes[i] = 1
                    elseif genotypes[i]==1
                        genotypes[i] = bool(sample([0,1])) ? 0 : 2
                    else
                        genotypes[i] = 1
                    end
                end
            end
            for locus in 1:common.G.chr[chr].numLoci
                genScore = (genotypes[locus]-1)*10
                @printf(genStream, "%4d", genScore)
            end
        end
        @printf(genStream, "\n")
    end
    close(genStream)
end

function outputCatData(fileName::AbstractString )
    catStream = open(fileName,"w")
    @printf(catStream, "BuildV1 Chr Pos Freq \n")
    snpID = 1
    for chr in 1:common.G.numChrom
        for locus in 1:common.G.chr[chr].numLoci
            mapPos = round(Int,common.G.chr[chr].mapPos[locus]*100000000)
            @printf(catStream, "%8d %3d %10d 0.5 \n", snpID, chr, mapPos)
            snpID += 1
        end
    end
    close(catStream)
end
