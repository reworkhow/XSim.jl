#generate haplotypes base allel frequency OR read haplotypes from files

function sampleFounders(file::String;header=false)
  numAnimals = round(Int,countlines(file)/2)
    sampleFounders(numAnimals,file,header=header)
end

function sampleFounders(numAnimals::Int64,file::String = "";fileName="",header=false)
  hapFile = false
  if file!=""
      hapFile = open(file)
  end
    my=Cohort(Array{Animal}(0),Array{Int64}(0,0))
    println("Sampling ",numAnimals," animals into base population.")
    resize!(my.animalCohort,numAnimals)
    for i in 1:numAnimals
        animal=sampleFounder(hapFile)
        my.animalCohort[i] = animal
        push!(common.founders,animal)
    end
    if fileName!=""
        outputPedigree(my,fileName)
    end
    return(my)
end

function sampleFounder(hapFile)
    my = Animal(0,0) #function Animal
    initFounderPosOri(my)
    initFounderHaps(my,hapFile)
    return(my)
end

function initFounderHaps(my::Animal,hapFile)
    numberChromosomePair=get_num_chrom(common.G)
    if hapFile != false
      hap1  =float(split(readline(hapFile))[2:end])
      hap2  =float(split(readline(hapFile))[2:end])
      k=1
      for i in 1:numberChromosomePair
        numLoci=common.G.chr[i].numLoci
        Base.resize!(my.genomePat[i].haplotype,numLoci)
        Base.resize!(my.genomeMat[i].haplotype,numLoci)
        for j in 1:numLoci
            my.genomePat[i].haplotype[j]=hap1[k]
            my.genomeMat[i].haplotype[j]=hap2[k]
            k=k+1
        end
      end
      return
    end
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

# Function to reset founders to be the animals in newFounders

function resetFounders(newFounders::Cohort)
    getOurHaps(newFounders)
    common.founders = newFounders.animalCohort
    common.countChromosome = 1
    for i in common.founders
        initFounderPosOri(i)
    end
end
