function sampleChildren(fathers::Cohort,mothers::Cohort,numAnimals::Int64)
    my=Cohort(Array(Animal,0),Array(Int64,0,0))

    #println("Sampling ",numAnimals," animals into next generation.")
    resize!(my.animalCohort,numAnimals)

    for i in 1:numAnimals
        animal=sampleNonFounder(getRandomInd(fathers),getRandomInd(mothers))
        #println("Sampled animal: ",animal.myID)
        my.animalCohort[i]=animal
    end

    return(my)
end

function getRandomInd(my::Cohort)
    cohortSize=length(my.animalCohort)
    thisOne=rand(1:cohortSize)
    return(my.animalCohort[thisOne])
end

function sampleNonFounder(father,mother)
    my=Animal(father.myID,mother.myID)

    numberChromosomePair=get_num_chrom(common.G)
    resize!(my.genomePat,numberChromosomePair)
    resize!(my.genomeMat,numberChromosomePair)

    sampleMyPosOri(my,father,mother)
    my.breedComp = (father.breedComp + mother.breedComp)/2
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
