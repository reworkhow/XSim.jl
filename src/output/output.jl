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
    nTraits = size(GLOBAL.LRes,2)
    for (k,animal) = enumerate(my.animals)
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
            @printf(pheStream, "%11.3f ",animal.val_p[i])
        end
        for i=1:nTraits
            @printf(pheStream, "%11.3f ",animal.val_g[i])
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
    set_genome(this)
    hapStream = open(fileName,"w")
    @printf(hapStream, "%10s ","animalID")
    for chr in 1:GLOBAL.G.n_chr
        for locus in 1:GLOBAL.G.chrs[chr].n_loci
            @printf(hapStream, "%2s", "p")
            @printf(hapStream, "%2s", "m")
        end
    end
    @printf(hapStream, "\n")
    for animal in this.animals
        @printf(hapStream, "%10d ", animal.myID)
        for chr in 1:GLOBAL.G.n_chr
            for locus in 1:GLOBAL.G.chrs[chr].n_loci
                @printf(hapStream, "%2d", animal.genome_sire[chr].haplotype[locus])
                @printf(hapStream, "%2d", animal.genome_dam[chr].haplotype[locus])
            end
        end
        @printf(hapStream, "\n")
    end
    close(hapStream)
end

function outputGenData(this::Cohort, fileName::AbstractString )
    set_genome(this)
    genStream = open(fileName,"w")
    @printf(genStream, "%10s ","animalID")
    snpID = 1
    for chr in 1:GLOBAL.G.n_chr
        for locus in 1:GLOBAL.G.chrs[chr].n_loci
            @printf(genStream, "%8d ", snpID)
            snpID += 1
        end
    end
    @printf(genStream, "\n")
    for animal in this.animals
        @printf(genStream, "%10d ", animal.myID)
        for chr in 1:GLOBAL.G.n_chr
            genotypes = animal.genome_sire[chr].haplotype + animal.genome_dam[chr].haplotype
            if GLOBAL.G.rate_error > 0.0
                numberErrors = int(GLOBAL.G.rate_error * GLOBAL.G.chrs[chr].n_loci)
                errorPos = sample([1:GLOBAL.G.chrs[chr].n_loci],numberErrors)
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
            for locus in 1:GLOBAL.G.chrs[chr].n_loci
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
    for chr in 1:GLOBAL.G.n_chr
        for locus in 1:GLOBAL.G.chrs[chr].n_loci
            mapPos = round(Int,GLOBAL.G.chrs[chr].mapPos[locus]*100000000)
            @printf(catStream, "%8d %3d %10d 0.5 \n", snpID, chr, mapPos)
            snpID += 1
        end
    end
    close(catStream)
end

function outputMarkerWithQTLEffects(fileName::AbstractString )
    nTraits = size(GLOBAL.LRes,2)

    totalLoci=0
    for i=1: GLOBAL.G.n_chr
        totalLoci=totalLoci+GLOBAL.G.chrs[i].n_loci
    end

    allMarkerQTLinfo = fill(0.0,totalLoci,nTraits)
    allMarkerQTLinfo[GLOBAL.G.qtl_index,:] = GLOBAL.G.qtl_effects

    catStream = open(fileName,"w")
    @printf(catStream, "CHROM\tPOS")
    for i in 1:nTraits
        @printf(catStream,"\tqtl_effect_%d",i)
    end
    @printf(catStream,"\n")
    snpID = 1
    for chr in 1:GLOBAL.G.n_chr
        for locus in 1:GLOBAL.G.chrs[chr].n_loci
            mapPos = round(Int,GLOBAL.G.chrs[chr].mapPos[locus]*100000000)
            @printf(catStream, "Chr%d\t%d",chr, mapPos)
            for i in 1:nTraits
                @printf(catStream,"\t%f",allMarkerQTLinfo[snpID,i])
            end
            @printf(catStream,"\n")
            snpID += 1
        end
    end
    close(catStream)
end