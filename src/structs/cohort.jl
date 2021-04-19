mutable struct Cohort
    n           ::Int16
    animals     ::Array{Animal, 1}
    is_founders ::Bool
    npMatrix    ::Array{Int64 , 2}

    # Null Constructor
    Cohort() = new(0)

    # Constructor for founders
    function Cohort(n         ::Int;
                    name_file ::String="")
        cohort = new(n,
                     Array{Animal}(undef, n),
                     true,
                     Array{Int64 }(undef, 0, 0))

        println("Sampling ", n, " animals into base population.")
        for i in 1:n
            cohort.animals[i] = Animal(true)
        end

        # if have file
        if name_file != ""
            outputPedigree(cohort, name_file)
        end

        return cohort
    end

    # Constructor for progeny
    function Cohort(animals   ::Array{Animal, 1};
                    is_founder::Bool=false)
        n = length(animals)
        return new(n,
                   animals,
                   is_founder,
                   Array{Int64 }(undef, 0, 0))
    end
end

# available types: phenotypic, genotypic, estimated
function get_traits(cohort::Cohort, type::String)
    traits_2d = (animal->get_traits(animal, type)).(cohort.animals)
    # return a n by p matrix
    return hcat(traits_2d...)'
end

function get_IDs(cohort::Cohort)
    return (animal->animal.ID).(cohort.animals)
end

function get_pedigree(cohort::Cohort)
    # return a 3-column matrix: ID, SireID, DamID
    ped = (animal->[animal.ID, animal.sire.ID, animal.dam.ID]).(cohort.animal)
    return hcat(ped...)'
end

"""
    DHs = get_double_haploids(parents::Cohort, nDHs::Int64)

* Produce *nDHs* double-haploids from a population **parents**.
"""
function get_DHs(parents::Cohort, n::Int64)
    println("Producing $n double-haploids from parents randomly
             selected from a population of size ", parents.n)

    animals = Array{Animal}(undef, n)
    select_idx = sample(parents.n, n, replace=true)
    for i in 1:nDHs
        parent = parents.animals[select_idx[i]]
        animals[i] = getDH(parent)
    end
    return cohort(animals)
end

function set_genome(cohort::Cohort)
    for animal in cohort.animals
        set_genomes(animal)
    end
end

# ped, mme, out is from JWAS get_pedigree(), build_model(), and solve()
function putEBV(cohort::Cohort, ped, mme, out)
    # transfer ebv from mme to XSim
    trmAnimal = mme.modelTermDict["1:Animal"]
    for animal in cohort.animals
        id = animal.ID
        strID = string(id)
        mmePos = ped.idMap[strID].seqID + trmAnimal.startPos - 1
        animal.traits[1].estimated = out[mmePos, 2]
    end
end


# # Public
# function concat_cohort(cohort::Cohort, concated_animal::Animal)
# end

# function concat_cohort(cohort::Cohort, concated_cohort::Cohort)
# end

# function get_g_incidence(cohort::Cohort)
# end

#include("nonfounders.jl")
include("get_children_from_manyparents.jl")
include("get1childfrom2parents.jl")
include("nonfounders_deprecated.jl")
include("checkCohort.jl")

# put mating function as constructor cohort

