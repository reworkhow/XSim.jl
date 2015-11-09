function sampleFounders(numAnimals::Int64)
    my=Cohort(Array{Animal,1}(),Array{Int64,2}())
    println("Sampling ",numAnimals," animals into base population.")
    resize!(my.animalCohort,numAnimals)
    for i in 1:numAnimals
        animal=sampleFounder()
        my.animalCohort[i] = animal
        push!(common.founders,animal)
    end
    return(my)
end

function sampleFounder()
    my = Animal(0,0) #function Animal
    initFounderPosOri(my)
    initFounderHaps(my)
    return(my)
end

function initFounderHaps(my::Animal)

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

function initFounderPosOri(my::Animal)
        numberChromosomePair=get_num_chrom(common.G)
        for i in 1:numberChromosomePair
            my.genomePat[i].ori=[common.countChromosome]
            my.genomePat[i].pos=[0.0]

            my.genomeMat[i].ori=[common.countChromosome+1]
            my.genomeMat[i].pos=[0.0]
        end
        common.countChromosome += 2
end
