function getIDs(animals::Cohort)
    animalIDs=Array{Int64,1}()
    for i in animals.animalCohort
        push!(animalIDs,i.myID)
    end
    return animalIDs
end

function getPedigree(animals::Cohort)
    pedMat=zeros(Int64,length(animals.animalCohort),3)
    for (i,j) in enumerate(animals.animalCohort)
        pedMat[i,1]=j.myID
        pedMat[i,2]=j.sireID
        pedMat[i,3]=j.damID
    end
    return pedMat
end

function putEBV(cohort,ped,mme,out)
    # transfer ebv from mme to XSim
    trmAnimal = mme.modelTermDict["1:Animal"]
    for animal in cohort.animalCohort
        id = animal.myID
        strID = string(id)
        mmePos = ped.idMap[strID].seqID + trmAnimal.startPos - 1
        animal.ebv = out[mmePos,2]
    end
end
