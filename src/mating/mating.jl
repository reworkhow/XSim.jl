##random mating
function sampleRan(popSize, nGen,sires,dams;gen=1,fileName="",printFlag=true)
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    mypopSize = round(Int,popSize/2)
    for i=1:nGen
        if printFlag==true
          @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i mypopSize mypopSize
        end
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

function mkPedArray(myPed::Array{Int64,2})
  pedArray = Array(XSim.PedNode,size(myPed,1));
  for i in 1:size(myPed,1)
    indi  = myPed[i,1]
    sirei = myPed[i,2]
    dami  = myPed[i,3]
    pedArray[myPed[i,1]] = XSim.PedNode(indi,sirei,dami)
  end
  return pedArray
end

function samplePed(myPed::Array{Int64,2})
  pedArray = mkPedArray(myPed)
  samplePed(pedArray)
end

function samplePed(ped::Array{PedNode,1})
    animals = Array(Animal,size(ped,1))
    hapFile = false
    for i in ped
        if i.ind <= i.sire || i.ind <= i.dam
            throw(Exception("ind < sire or dam \n"))
        end
        if i.sire == 0
            animal = sampleFounder(hapFile)
            animals[i.ind] = animal
            push!(common.founders,animal)
        else
            animals[i.ind] = sampleNonFounder(animals[i.sire],animals[i.dam])
        end
    end
    res = Cohort(animals,Array(Int64,0,0))
end

function samplePed(myPed::Array{Int64,2},animalVec::Cohort)
    pedArray = mkPedArray(myPed)
    samplePed(pedArray,animalVec)
end

function samplePed(ped::Array{PedNode,1},animalVec::Cohort)
    atFounder = 1
    founders  = XSim.copy(animalVec)

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


function sampleSel(popSize, nSires, nDams, nGen, varRes=common.varRes)
    maleCandidates   = sampleFounders(round(Int,popSize/2))
    femaleCandidates = sampleFounders(round(Int,popSize/2))
    return sampleSel(popSize, nSires, nDams, nGen,maleCandidates,femaleCandidates, varRes)
end

function sampleSel(popSize, nSires, nDams, nGen,maleParents,femaleParents,varRes=common.varRes;gen=1,fileName="", direction=1)
    common.varRes    = varRes # common.varRes is used in outputPedigree(cohort,fileName)
    maleCandidates   = copy(maleParents)
    femaleCandidates = copy(femaleParents)
    sires = Cohort(Array(Animal,0),Array(Int64,0,0))
    dams  = Cohort(Array(Animal,0),Array(Int64,0,0))
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    for i=1:nGen
        @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i round(Int,popSize/2) round(Int,popSize/2)
        y = direction*getOurPhenVals(maleCandidates,varRes)
        sires.animalCohort = maleCandidates.animalCohort[sortperm(y)][(end-nSires+1):end]
        y = direction*getOurPhenVals(femaleCandidates,varRes)
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

function sampleBLUPSel(popSize, nSires, nDams, nGen,maleParents,femaleParents,varRes=common.varRes,varGen=1;gen=1,fileName="XSim", direction=1)
    # initial BLUP evaluation--- parents
    run(`\rm -f $fileName.ped`)
    run(`\rm -f $fileName.phe`)
    run(`\rm -f $fileName.brc`)
    run(`\rm -f $fileName.gen`)
    outputPedigree(maleParents,fileName)
    outputPedigree(femaleParents,fileName)
    pedfile   = fileName*".ped"
    phenofile = fileName*".phe"
    colNames = [:Animal;:y;:bv]
    dfPhen = readtable(phenofile,separator = ' ',header=false,names=colNames)
    ped = get_pedigree(pedfile)
    mme = build_model("y = intercept + Animal",varRes)
    set_random(mme,"Animal",ped,varGen)
    out = solve(mme,dfPhen,solver="GaussSeidel",printout_frequency=40)
    # trasfer BLUP-EBV to animals
    putEBV(maleParents,ped,mme,out)
    putEBV(femaleParents,ped,mme,out)

    maleCandidates   = copy(maleParents)
    femaleCandidates = copy(femaleParents)
    sires = Cohort(Array(Animal,0),Array(Int64,0,0))
    dams  = Cohort(Array(Animal,0),Array(Int64,0,0))
    boys  = Cohort(Array(Animal,0),Array(Int64,0,0))
    gals  = Cohort(Array(Animal,0),Array(Int64,0,0))
    for i=1:nGen
        @printf "Generation %5d: sampling %5d males and %5d females\n" gen+i round(Int,popSize/2) round(Int,popSize/2)
        y = direction*[animal.ebv for animal in maleCandidates.animalCohort]
        sires.animalCohort = maleCandidates.animalCohort[sortperm(y)][(end-nSires+1):end]
        y = direction*[animal.ebv for animal in femaleCandidates.animalCohort]
        dams.animalCohort = femaleCandidates.animalCohort[sortperm(y)][(end-nDams+1):end]
        boys = sampleChildren(sires,dams,round(Int,popSize/2))
        gals = sampleChildren(sires,dams,round(Int,popSize/2))
        outputPedigree(boys,fileName)
        outputPedigree(gals,fileName)
        maleCandidates.animalCohort   = [sires.animalCohort; boys.animalCohort]
        femaleCandidates.animalCohort = [dams.animalCohort;  gals.animalCohort]

        # BLUP ebvs
        dfPhen = readtable(phenofile,separator = ' ',header=false,names=colNames)
        ped = get_pedigree(pedfile)
        mme = build_model("y = intercept + Animal",varRes)
        set_random(mme,"Animal",ped,varGen)
        out = solve(mme,dfPhen,solver="GaussSeidel",printout_frequency=40)
        # trasfer BLUP-EBV to animals
        putEBV(maleCandidates,ped,mme,out)
        putEBV(femaleCandidates,ped,mme,out)
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
