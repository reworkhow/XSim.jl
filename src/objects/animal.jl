mutable struct Animal
    ID          ::Int64
    sire        ::Animal
    dam         ::Animal
    genome_sire ::Array{Chromosome, 1}
    genome_dam  ::Array{Chromosome, 1}
    val_p       ::Array{Float64,    1}
    val_g       ::Array{Float64,    1}
    is_founder  ::Bool

    # Null constructor
    Animal() = new(0)

    # Constructor for founders
    function Animal(is_founder::Bool)
        # instantiate a founder
        return Animal(Animal(), Animal(), is_founder=is_founder)
    end

    # Constructor
    function Animal(sire      ::Animal,
                    dam       ::Animal;
                    is_founder::Bool=false)

        animal = new(GLOBAL("count_id"), sire, dam,
                     Array{Chromosome}(undef, GLOBAL("n_chr")),
                     Array{Chromosome}(undef, GLOBAL("n_chr")),
                     Array{Float64   }(undef, 0),
                     Array{Float64   }(undef, 0),
                     is_founder)
        add_count_ID!(by=1)

        # Setup genome
        if is_founder
            set_genome!(animal)
            add_founder!(animal)
        else
            set_genome!(animal, dam, sire)
        end

        # Compute breeding values
        set_BV!(animal)

        return animal
    end
end

function set_BV!(animal::Animal)
    set_haplotypes!(animal)
    genotypes = get_genotypes(animal)
    animal.val_g = (genotypes'GLOBAL("effects"))'
end


function set_genome!(animal::Animal)
    # Founders
    for i in 1:GLOBAL("n_chr")
        animal.genome_sire[i] = Chromosome(i, GLOBAL("count_hap"))
        animal.genome_dam[i]  = Chromosome(i, GLOBAL("count_hap") + 1)
    end
    add_count_haplotype!(by=2)
end

function set_genome!(animal::Animal, dam::Animal, sire::Animal)
    # Progenies
    for i in 1:GLOBAL("n_chr")
        animal.genome_sire[i] = Chromosome(i, sire)
        animal.genome_dam[i]  = Chromosome(i, dam)
    end
end

function set_haplotypes!(animal::Animal)
    set_haplotypes!(animal.genome_sire)
    set_haplotypes!(animal.genome_dam)
end

function get_traits(animal::Animal,
                     option::String="Ve",
                     values::Union{Array{Float64}, Float64})
    n_traits = GLOBAL("n_traits")
    Vg       = GLOBAL("Vg")
    if option == "h2"
        if n_traits > 1 & !isa(values, Array)
            values = fill(values, n_traits)
        end
        Ve = ((ones(n_traits) .- values) .* diag(Vg)) ./ values
        Ve = n_traits == 1 ? Ve[1] : Ve

    elseif option == "Ve"
        Ve = handle_diagonal(values, n_traits)
    end

    animal.val_p = animal.val_g + Ve * randn(n_traits)

    return animal.val_p
end


function get_DH(individual::Animal)
    return Animal(individual, individual)
end

function get_genotypes(animal::Animal)
    myGenotype = Array{AlleleIndexType}(undef, 0)
    for i in 1:GLOBAL("n_chr")
        append!(myGenotype,
                animal.genome_sire[i].haplotype + animal.genome_dam[i].haplotype)
    end

    return myGenotype
end


function Base.:+(x::Animal, y::Animal)
    return Cohort([x, y])
end


