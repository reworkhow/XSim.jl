

function find_haplotype_founder(animals::Cohort)
    haps=[]
    for i in animals.animalCohort
        for j in [i.genomePat;i.genomeMat]
            for k in j.ori
                push!(haps,k)
            end
        end
    end
    println("Identifying founders for haplotype segments ")
    haps
end
