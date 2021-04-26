mutable struct Cohort
    animals     ::Array{Animal, 1}
    n           ::Int64

    # Null Constructor
    Cohort() = new()

    # Constructor for founders
    function Cohort(n::Int64)
        cohort = new(Array{Animal}(undef, n), n)
        for i in 1:n
            cohort.animals[i] = Animal(true)
        end

        return cohort
    end

    # Constructor for progeny
    function Cohort(animals::Array{Animal, 1})
        n = length(animals)
        return new(animals, n)
    end

    function Cohort(animal::Animal)
        return new([animal], 1)
    end

end

function length(cohort::Cohort)
    return length(cohort.animals)
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
    ped = (animal->[animal.ID, animal.sire.ID, animal.dam.ID]).(cohort.animals)
    return hcat(ped...)'
end


function get_DHs(parents::Cohort, n::Int64)
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


##mating with breed components
function setBreedComp(c::Cohort, comp::Array{Float64,1})
    for animal in c.animals
        animal.breedComp = comp
    end
end


function Base.:+(x::Cohort, y::Cohort)
    return Cohort(vcat(x.animals, y.animals))
end

function Base.:+(x::Cohort, y::Animal)
    return Cohort(vcat(x.animals, y))
end

function Base.:+(x::Animal, y::Cohort)
    return Cohort(vcat(x, y.animals))
end

function sample(cohort::Cohort, n::Int64; replace::Bool=true)
    select = sample(1:cohort.n, n, replace=replace)
    return cohort[select]
end

function getindex(cohort::Cohort, I...)
    if length(cohort) == 1
        return getindex(cohort.animals, I...)
    else
        return Cohort(getindex(cohort.animals, I...))
    end
end

function print(cohort::Cohort, type::String="ID")
    if type == "ID"
        for animal in cohort.animals
            println("Individual: ", animal.ID)
        end
    elseif type == "Pedigree"
        return get_pedigree(cohort)
    end
end
Base.show(io::IO, z::Cohort) = print(z)


# function sample(pool::Int64, n::Int64; replace::Bool=false)
#     # select n samples from 1:pool numbers, return a size of n 1-d array
#     if replace
#         samples = rand(1:pool, n)
#     else
#         n = n > pool ? pool : n
#         samples = shuffle(1:pool)[1:n]
#     end

#     return n == 1 ? samples[1] : samples
# end



# # Public
# function concat_cohort(cohort::Cohort, concated_animal::Animal)
# end

# function concat_cohort(cohort::Cohort, concated_cohort::Cohort)
# end

# function get_g_incidence(cohort::Cohort)
# end


# put mating function as constructor cohort

