#genotyping a cohort
function getOurGenotypes(my::Cohort,sel::Array{Int64,1})
    cohort = cohortSubset(my,sel)
    return getOurGenotypes(cohort)
end

function getOurGenotypes(my::Cohort)
    getOurHaps(my)

    nLoci=0
    for i=1: common.G.numChrom
        nLoci=nLoci+common.G.chr[i].numLoci
    end

    npMatrix=Array{Int64}(undef,length(my.animalCohort), nLoci)
    for (i,value) in enumerate(my.animalCohort)
        npMatrix[i,:]=getMyGenotype(value)
    end
    return npMatrix
end


#detailed functions
function getOurHaps(my::Cohort)
    for i in my.animalCohort
        getMyHaps(i)
    end
end

function getMyGenotype(my)
    myGenotype=Array{Int64}(undef,0)
    for i in 1:common.G.numChrom
        append!(myGenotype, my.genomePat[i].haplotype+my.genomeMat[i].haplotype)
    end
    return myGenotype
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
            whichFounder=ceil(Integer,genome[i].ori[segment]/2)
            genomePatorMatInThisFounder=(genome[i].ori[segment]%2==0) ? common.founders[whichFounder].genomeMat[i] : common.founders[whichFounder].genomePat[i]

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

        for j in 1:length(genome[i].mut)
            whichlocus = findfirst(common.G.chr[i].mapPos .== genome[i].mut[j])
            genome[i].haplotype[whichlocus] = 1 - genome[i].haplotype[whichlocus]
        end

        pop!(genome[i].pos)
    end
end
