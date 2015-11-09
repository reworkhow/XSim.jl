##random mating
function sampleRan(popSize, nGen,sires,dams;gen=1,fileName="")
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    mypopSize = round(Int,popSize/2)
    for i=1:nGen
        @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i mypopSize mypopSize
        boys = sampleChildren(sires,dams,mypopSize)
        gals = sampleChildren(sires,dams,mypopSize)
        sires = boys
        dams  = gals
        if fileName!=""
            outputPedigree(boys,fileName)
            outputPedigree(gals,fileName)
        end
    end
    gen += nGen
    return boys,gals, gen
end

## mating individuals for a given pedigree
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

function samplePed(ped::Array{PedNode,1},animalVec)
    atFounder = 1
    founders  = copy(animalVec)
    animals = Array(Animal,size(ped,1))
    for i in ped
        if i.ind <= i.sire || i.ind <= i.dam
            throw(Exception("ind < sire or dam \n"))
        end
        if i.sire == 0
            animal = founders.animalCohort[atFounder]
            atFounder += 1
            animals[i.ind] = animal
        else
            animals[i.ind] = sampleNonFounder(animals[i.sire],animals[i.dam])
        end
    end
    res = Cohort(animals,Array(Int64,0,0))
end

##mating with selections
function sampleSel(popSize, nSires, nDams, nGen, varRes)
    maleCandidates   = sampleFounders(round(Int,popSize/2))
    femaleCandidates = sampleFounders(round(Int,popSize/2))
    return sampleSel(popSize, nSires, nDams, nGen,maleCandidates,femaleCandidates, varRes)
end

function sampleSel(popSize, nSires, nDams, nGen,maleParents,femaleParents,varRes;gen=1,fileName="", direction=1)
    maleCandidates   = copy(maleParents)
    femaleCandidates = copy(femaleParents)
    sires = Cohort(Array(Animal,0),Array(Int64,0,0))
    dams  = Cohort(Array(Animal,0),Array(Int64,0,0))
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    for i=1:nGen
        @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i round(Int,popSize/2) round(Int,popSize/2)
        y = direction*getOurPhenVals(maleCandidates,common.varRes)
        sires.animalCohort = maleCandidates.animalCohort[sortperm(y)][(end-nSires+1):end]
        y = direction*getOurPhenVals(femaleCandidates,common.varRes)
        dams.animalCohort = femaleCandidates.animalCohort[sortperm(y)][(end-nDams+1):end]
        boys = sampleChildren(sires,dams,round(Int,popSize/2))
        gals = sampleChildren(sires,dams,round(Int,popSize/2))
        if fileName!=""
            outputPedigree(boys,fileName)
            outputPedigree(gals,fileName)
        end
        maleCandidates.animalCohort   = [sires.animalCohort; boys.animalCohort]
        femaleCandidates.animalCohort = [dams.animalCohort;  gals.animalCohort]
    end
    gen += nGen
    return boys,gals, gen
end

function copy(c::Cohort)
	return Cohort(Base.copy(c.animalCohort), Base.copy(c.npMatrix) )
end

##mating with breed components
function setBreedComp(c::Cohort,comp::Array{Float64,1})
    for animal in c.animalCohort
        animal.breedComp = comp
    end
end

##concat several cohorts
function concatCohorts(cohortLst...)
    # returns a cohort with concatenation of the animalCohorts from the arguments
    res = Cohort(Array{Animal,1}(),Array{Int64,2}())
    for i in cohortLst
        res.animalCohort = [res.animalCohort; i.animalCohort]
    end
    return res
end

#get a subset of a cohort
function cohortSubset(my::Cohort,sel::Array{Int64,1})
    animals = Array(Animal,size(sel,1))
    for (i,j) = enumerate(sel)
        animals[i] = my.animalCohort[j]
    end
    return Cohort(animals,Array(Int64,0,0))
end
