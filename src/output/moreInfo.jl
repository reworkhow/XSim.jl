

function find_haplotype_founder(animals::Cohort)
    haps=[]
    for i in animals.animals
        for j in [i.genome_sire;i.genome_dam]
            for k in j.ori
                push!(haps,k)
            end
        end
    end
    println("Identifying founders for haplotype segments ")
    haps
end
