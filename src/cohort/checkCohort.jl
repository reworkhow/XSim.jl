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