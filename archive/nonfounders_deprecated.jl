function sampleChildren(fathers::Cohort, mothers::Cohort, numAnimals::Int64;
                        fWoRepl::Bool=false, mWoRepl::Bool=false,
                        ET::Bool=false, numOffET::Int64=1)
    my=Cohort(Array{Animal}(undef,0),Array{Int64}(undef,0,0))
    fatherDict = Dict{Int64,Int64}()
    motherDict = Dict{Int64,Int64}()
    println("RUNNING OLD FUNCTIONNNNN")
    #println("Sampling ",numAnimals," offspring into next generation.")
    resize!(my.animals,numAnimals)
    if(ET == false)
        for i in 1:numAnimals
                if(fWoRepl == true)
                    if(mWoRepl == true)
                        mother = getRandomIndWoReplacement(mothers,i)
                        father = getRandomIndWoReplacement(fathers,i)
                    else
                        mother = getRandomInd(mothers)
                        father = getRandomIndWoReplacement(fathers,i)
                    end
                elseif(mWoRepl == true)
                    father = getRandomInd(fathers)
                    mother = getRandomIndWoReplacement(mothers,i)
                else
                    mother = getRandomInd(mothers)
                    father = getRandomInd(fathers)
                end
                motherDict[mother.myID] = get(motherDict, mother.myID, 0) + 1
                fatherDict[father.myID] = get(fatherDict, father.myID, 0) + 1
                animal=sampleNonFounder(father,mother)
                #println("Sampled animal: ",animal.myID)
                my.animals[i]=animal
        end
    else
        numMatings::Int64 = div(numAnimals, numOffET)
        if(rem(numAnimals, numOffET) != 0)
            println("get_children(): Number of offspring to be generated is ", numAnimals)
            println("get_children(): Number of offspring per ET dam is ", numOffET)
            println("get_children(): Number of matings is not an integer ", numAnimals/numOffET)
            numMatings+=1
        end
        for i in 1:numMatings
            mother = getRandomIndWoReplacement(mothers,i)
            if(i == numMatings)
                if(rem(numAnimals, numOffET) != 0)
                    numOffET = numAnimals - ((numMatings-1)*numOffET)
                    println("get_children(): Last mating (",numMatings,") is producing ", numOffET," ET-offspring.")
                end
            end
            for j in 1:numOffET
                if(fWoRepl == true)
                    father = getRandomIndWoReplacement(fathers,i)
                else
                    father = getRandomInd(fathers)
                end
                motherDict[mother.myID] = get(motherDict, mother.myID, 0) + 1
                fatherDict[father.myID] = get(fatherDict, father.myID, 0) + 1
                animal=sampleNonFounder(father,mother)
                #println("Sampled animal: ",animal.myID)
                my.animals[j+(numMatings * (i-1))]=animal
            end
        end
    end
    println("Number of fathers used: ", length(fatherDict))
    println("Number of mothers used: ", length(motherDict))
    println("Number of offspring generated: ", length(my.animals))
    return(my)
end




function sampleOffAllMatings(fathers::Cohort, mothers::Cohort, numAnimalsPerMating::Int64)
    my=Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))

    numAnimals = length(fathers.animals)*length(mothers.animals)*numAnimalsPerMating

    println("sampleOffAllMatings(): Sampling ",numAnimalsPerMating," offspring per mating into next generation.")
    println("There are a total of ",length(fathers.animals)*length(mothers.animals),"  matings and thus ", numAnimals, "  offspring in the next generation.")
    numMothers = length(mothers.animals)
    numFathers = length(fathers.animals)
    resize!(my.animals,numAnimals)
    for m in 1:numMothers
        for f in 1:numFathers
            animal=sampleNonFounder(fathers.animals[f],mothers.animals[m])
            #println("Sampled animal: ",animal.myID)
            my.animals[((m-1)*numFathers) + f]=animal
            if (((m-1)*numFathers) + f)%10000 == 0
                println(((m-1)*numFathers) + f," offspring processed.")
            end
        end
    end
    return(my)
end


function getRandomIndWoReplacement(my::Cohort,i::Int64)
   cohortSize=length(my.animals)
   if(i>=cohortSize)
     println("getRandomIndWoReplacement() No animal left in cohort to sample.")
     exit(1)
   end
   thisOne=rand(i:cohortSize)
   tempAnimal = my.animals[i]
   my.animals[i] = my.animals[thisOne]
   my.animals[thisOne] = tempAnimal
   #println("Sampled animal ID=",my.animals[i].myID," from a total of ",cohortSize-i+1," animals.")
   return(my.animals[i])
end

function getRandomSampleOfIndWoReplacement!(my::Cohort,num::Int64)
    println("getRandomSampleOfIndWoReplacement!()")
    cohortSize=length(my.animals)
    if num>=cohortSize
        #println("getRandomSampleOfIndWoReplacement!(): Number of animals to sample is larger than the cohort size.")
        #exit(1)
    end
    for i in 1:num
        thisOne=rand(i:cohortSize)
        tempAnimal = my.animals[i]
        my.animals[i] = my.animals[thisOne]
        my.animals[thisOne] = tempAnimal
        if i%5000 == 0
            println(i," individuals of total ", num, " individuals sampled.")
        end
    end
    println("getRandomSampleOfIndWoReplacement!(): Sampled ",num," animals randomly out of ", cohortSize)

    return my
end


function getSampleWOutRepl(fromCohort::Cohort, numAnimals::Int64)
    # sampling without replacement
    fromSize = length(fromCohort.animals)
    if numAnimals > fromSize
        println("In getSampleWOutRepl():  numAnimals > size of fromCohort")
        exit(1)
    end
    fromIndx = sample(1:fromSize,numAnimals,replace=false)
    toCohort=Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    resize!(toCohort.animals,numAnimals)
    toCohort.animals = fromCohort.animals[fromIndx]
    # for i=1:numAnimals
    #     toCohort.animals[i] = fromCohort.animals[fromIndx[i]]
    # end
    return(toCohort)
end

function getSampleWOutRepl!(fromCohort::Cohort, numAnimals::Int64)
    # sampling without replacement
    fromSize = length(fromCohort.animals)
    if numAnimals > fromSize
        println("In getSampleWOutRepl!(): numAnimals > size of cohort to sample from")
        exit(1)
    end
    println("getSampleWOutRepl!(): Sampling ", numAnimals," from cohort of orginal size ", length(fromCohort.animals))
    fromIndx = sort(sample(1:fromSize,numAnimals,replace=false))
    toCohort=Cohort(Array{Animal,1}(undef,0),Array{Int64,2}(undef,0,0))
    resize!(toCohort.animals,numAnimals)
    sizeCohort = length(fromCohort.animals)
    println("Ready to sample from 'fromCohort', size of fromCohort ",length(fromCohort.animals)," size of toCohort ", length(toCohort.animals)," idx = ", fromIndx)
    toCohort.animals = fromCohort.animals[fromIndx]
    #fromCohort.animals = deleteat!(vec(fromCohort.animals), fromIndx)
    println("getSampleWOutRepl!(): ", numAnimals," animals sampled from cohort of orginal size ", sizeCohort,"; new size is ", length(fromCohort.animals))
    return(toCohort)
end

function sampleOneDHOffspringFrom(parent::Animal)
    offspring = Animal(parent.myID,0)
    sampleOnePosOri(offspring.genome_sire,parent)
    offspring.genome_dam = deepcopy(offspring.genome_sire)
    return offspring
end

function sampleDHOffspringFrom(parents::Cohort, numDHOffs::Int64)
    println("Sampling $numDHOffs offspring from parents selected at random from a cohort of size ",size(parents.animals,1))
    offspring=Cohort(Array{Animal}(undef,0),Array{Int64}(undef,0,0))
    resize!(offspring.animals,numDHOffs)
    for i in 1:numDHOffs
        parent = getRandomInd(parents)
        offspring.animals[i] = sampleOneDHOffspringFrom(parent)
    end
    return offspring
end



function getRandomInd(my::Cohort)
    cohortSize=length(my.animals)
    thisOne=rand(1:cohortSize)
    return(my.animals[thisOne])
end

function sampleNonFounder(father,mother)
    my=Animal(father.myID,mother.myID)

    numberChromosomePair=get_num_chrom(common.G)
    resize!(my.genome_sire,numberChromosomePair)
    resize!(my.genome_dam,numberChromosomePair)

    sampleMyPosOri(my,father,mother)
    my.breedComp = (father.breedComp + mother.breedComp)/2
    return(my)
end

function sampleMyPosOri(my,father,mother)
    sampleOnePosOri(my.genome_sire,father)
    sampleOnePosOri(my.genome_dam,mother)
end

function sampleOnePosOri(genome::Array{Chromosome,1},parent::Animal)
    numberChromosomePair=get_num_chrom(common.G)

    for i in 1:numberChromosomePair

        genome[i]=Chromosome(Array{Int64}(undef,0),Array{Int64}(undef,1),Array{Float64}(undef,1),Array{Float64}(undef,1))

        currentChrom=(rand(Bernoulli(0.5))==1) ? parent.genome_sire[i] : parent.genome_dam[i]

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
        numTempMut=0
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
            for k in 1:length(currentChrom.mut)
              if currentChrom.mut[k] >= rec[j] && currentChrom.mut[k] < rec[j+1]
                  numTempMut = numTempMut+1
                  tempMut[numTempMut] = currentChrom.mut[k]
              elseif currentChrom.mut[k]>=rec[j+1]
                  break
              end
            end

            currentChrom=(currentChrom==parent.genome_dam[i]) ? parent.genome_sire[i] : parent.genome_dam[i]

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

        #add new mutation (no merge required as in pos/ori below)
        resize!(genome[i].mut,numTempMut)

        for muti in 1:numTempMut
            genome[i].mut[muti]=tempMut[muti]
        end

        mutation_rate = common.G.mutRate
        numLoci       = common.G.chr[i].numLoci
        nmut          = rand(Binomial(numLoci,mutation_rate))
        if nmut != 0
            muts = sample(common.G.chr[i].mapPos,nmut)
            genome[i].mut = vcat(genome[i].mut,muts)
            sort!(genome[i].mut)
        end

        ##merging of nearby same ori;not sure whether more effcient

        genome[i].pos[1]=tempPos[1]
        genome[i].ori[1]=tempOri[1]

        this=1
        for m in 2:numTemp
            if(tempOri[m]!=tempOri[m-1])
                this=this+1
                genome[i].pos[this]=tempPos[m]
                genome[i].ori[this]=tempOri[m]
            end
        end

        resize!(genome[i].pos,this) #resize after merging
        resize!(genome[i].ori,this)
    end
end
