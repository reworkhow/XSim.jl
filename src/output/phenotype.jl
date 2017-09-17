function getOurPhenVals(my::Cohort, varRes)
    stdRes = sqrt(varRes)
    n = size(my.animalCohort,1)
    phenVals = Array{Float64}(n)
    genVals  = getOurGenVals(my)
    for (i,animal) = enumerate(my.animalCohort)
        if animal.phenVal==-9999
            animal.phenVal =  animal.genVal + (randn(1)*stdRes)[1]
        end
        phenVals[i] = animal.phenVal
    end
    return phenVals
end

function getOurGenVals(my::Cohort)
    n = size(my.animalCohort,1)
    genVals = Array{Float64}(n)
    for (i,animal) = enumerate(my.animalCohort)
        #println(i)
        if animal.genVal==-9999
            getMyHaps(animal)
            myGenotypes = getMyGenotype(animal)
            animal.genVal = dot(myGenotypes[common.G.qtl_index],common.G.qtl_effects)
        end
        genVals[i] = animal.genVal
    end
    return genVals
end
