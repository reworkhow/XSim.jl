mutable struct Animal
    ID          ::Int64
    sire        ::Animal
    dam         ::Animal
    genome_sire ::Array{Chromosome, 1}
    genome_dam  ::Array{Chromosome, 1}
    traits      ::Array{Trait     , 1}
    n_traits    ::Int
    breedComp   ::Array{Float64   , 1}
    is_founder  ::Bool

    # Null constructor
    Animal() = new(0)

    # Constructor for founders
    function Animal(is_founder::Bool)
        # instantiate a founder
        return Animal(Animal(), Animal(), is_founder=is_founder)
    end

    # Constructor for progenies
    function Animal(sire      ::Animal,
                    dam       ::Animal;
                    is_founder::Bool=false)

        animal = new(GLOBAL.countId, sire, dam,
                     Array{Chromosome}(undef, GLOBAL.G.numChrom),
                     Array{Chromosome}(undef, GLOBAL.G.numChrom),
                     Array{Trait     }(undef, 0),
                     0,
                     Array{Float64   }(undef, 0),
                     is_founder)
        GLOBAL.countId += 1

        if is_founder
            set_genome!(animal)
            push!(GLOBAL.founders, animal)
        else
            sampleOnePosOri(animal.genome_sire, sire)
            sampleOnePosOri(animal.genome_dam,  dam)
            animal.breedComp = (sire.breedComp + dam.breedComp)/2
        end

        return animal
    end

    function set_genome!(animal::Animal)
        is_founder = animal.is_founder
        for i in 1:GLOBAL.G.numChrom
            animal.genome_sire[i] = Chromosome(i,
                                               ori=GLOBAL.countChromosome,
                                               is_founder=is_founder)
            animal.genome_dam[i]  = Chromosome(i,
                                               ori=GLOBAL.countChromosome + 1,
                                               is_founder=is_founder)
        end
        GLOBAL.countChromosome += 2
    end
end


function Base.:+(x::Animal, y::Animal)
    return Cohort([x, y])
end


# available types: phenotypic, genotypic, estimated
function get_traits(animal::Animal, type::String)
    return (traits -> getfield(traits, Symbol(type))).(animal.traits)
end


"""
    DH = get_DH(individual::Animal)

* Produce one double-haploid from the **individual**.
"""
function get_DH(individual::Animal)
    progeny = Animal(individual, Animal())
    sampleOnePosOri(progeny.genome_sire, individual)
    progeny.genome_dam = deepcopy(progeny.genome_sire)
    return progeny
end



