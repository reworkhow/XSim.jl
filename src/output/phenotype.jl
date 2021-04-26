function getOurPhenVals(my::Cohort, varRes=0)
    if varRes == 0
       if Base.length(GLOBAL.LRes) != 0
           LRes = GLOBAL.LRes
        else
           error("getOurPhenVals(): GLOBAL.LRes is not initialized")
        end
    else
       if typeof(varRes)==Float64
          LRes = cholesky(fill(varRes,1,1)).U'
       else
          LRes = cholesky(varRes).U'
        end
    end

    n = size(my.animals, 1)
    nTraits = size(GLOBAL.LRes, 2)
    phenVals = Array{Float64,2}(undef, n, nTraits)
    genVals  = getOurGenVals(my, nTraits)
    for (i, animal) = enumerate(my.animals)
        if length(animal.val_p) == 0
            animal.phenVal =  animal.val_g + LRes * randn(nTraits)
        end
        phenVals[i,:] = animal.val_p
    end
    return phenVals
end

function getOurGenVals(my::Cohort,nTraits=0)
    if nTraits == 0
        if  size(GLOBAL.LRes, 2) != 0
            nTraits = size(GLOBAL.LRes, 2)
        else
            error("getOurGenVals(): GLOBAL.LRes is not initialized,
                   cannot get the number of traits from the 2nd dimension of
                   the residual variance matrix.")
        end
    end

    n = size(my.animals, 1)
    genVals = Array{Float64}(undef, n, nTraits)
    for (i, animal) = enumerate(my.animals)
        if i % 1000 == 0
            println("getOurGenVals(): ", i)
        end
        if animal.nTraits == 0
            set_genomes(animal)
            myGenotypes = getMyGenotype(animal)
            resize!(animal.traits, nTraits)

            animal.val_g = (myGenotypes[GLOBAL.G.qtl_index]'GLOBAL.G.qtl_effects)'
	    end
        genVals[i,:] = animal.val_g
    end
    return genVals
end
