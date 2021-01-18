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

    npMatrix=Array{AlleleIndexType}(undef,length(my.animalCohort), nLoci)
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
    myGenotype=Array{AlleleIndexType}(undef,0)
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
        push!(genome[i].pos,common.G.chr[i].chrLength)

        iLocus   = 1
        position = common.G.chr[i].mapPos[iLocus]
        #println("getOneHaps(): number of segments in chromosome ",i,": ",numOri)
        chrom = common.G.chr[i]
        lociPerM = round(Int64, numLoci/genome[i].pos[numOri+1])
        segLen = 0
        prevSegLen = 0
        endLocus = 0

        for segment in 1:numOri
            #println("getOneHaps(): segment=",segment)
            flush(stdout)
            whichFounder=ceil(Integer,genome[i].ori[segment]/2)
            genomePatorMatInThisFounder=(genome[i].ori[segment]%2==0) ? common.founders[whichFounder].genomeMat[i] : common.founders[whichFounder].genomePat[i]

            startPos = genome[i].pos[segment]
            endPos   = genome[i].pos[segment+1]
            prevSegLen += segLen
            segLen   = 0
            if segment < numOri
               numLociUntilGuessedPos = clamp(round(Int64, endPos * lociPerM),1::Int64, numLoci)
               guessedPos = chrom.mapPos[numLociUntilGuessedPos]
               if guessedPos > endPos
                  ul = numLociUntilGuessedPos
                  ll = iLocus
                    if ll > endPos
                        ll -=1
                    end
               else
                  ll = numLociUntilGuessedPos
                  ul = numLoci
               end

               iter = 0
               while ul-ll > 1
                    #println("getOneHaps(): iter=",iter," endPos=",endPos," guessedPos=",guessedPos," ll=",ll," ul=",ul," segLen=",ul-ll)
                    prevNumLociUntilGuessedPos = numLociUntilGuessedPos
                    iter+=1
                    numLociUntilGuessedPos = (ul-ll)/2
                    numLociUntilGuessedPos = ll + round(Int64, numLociUntilGuessedPos)
                    guessedPos = chrom.mapPos[numLociUntilGuessedPos]
                    if prevNumLociUntilGuessedPos == numLociUntilGuessedPos
                        if guessedPos == ll
                            if chrom.mapPos[ll+1] < endPos
                                ll += 1
                            else
                                ul -= 1
                            end
                        else
                            if chrom.mapPos[ul-1] > endPos
                                ul -= 1
                            else
                                ll += 1
                            end
                        end
                    elseif guessedPos > endPos
                        ul = numLociUntilGuessedPos
                    else
                        ll = numLociUntilGuessedPos
                    end
                end
                endLocus = ll
                segLen = ll - prevSegLen
                if iLocus < numLoci
                    iLocus = ll+1
                end
                #println("getOneHaps(): iter=",iter," endPos=",endPos," guessedPos=",guessedPos," ll=",ll," ul=",ul," segLen=",segLen," iLocus=",iLocus)
             elseif iLocus <= numLoci
                segLen = numLoci - endLocus
                endLocus = numLoci
                #println("getOneHaps(): last segment, segLen=",segLen," iLocus=",iLocus," endLocus=",endLocus)
             else
                segLen = 0
                #println("getOneHaps(): iLocus = ",iLocus,", numLoci = ",numLoci,", startPos=",startPos," endPos=",endPos,", chrLength=",common.G.chr[i].chrLength)
             end
             #println("getOneHaps(): segment=",segment," endLocus=",endLocus, " seglen=", segLen, " prevSeglen=", prevSegLen)
             if segLen > 0
               genome[i].haplotype[(endLocus-segLen+1):endLocus]=genomePatorMatInThisFounder.haplotype[(endLocus-segLen+1):endLocus]
            end
        end


        for j in 1:length(genome[i].mut)
            whichlocus = findfirst(common.G.chr[i].mapPos .== genome[i].mut[j])
            genome[i].haplotype[whichlocus] = 1 - genome[i].haplotype[whichlocus]
        end

        pop!(genome[i].pos) #remove temporary added chrLength
    end
end
