mutable struct Animal <: AbstractAnimal
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

function get_BVs(animal::Animal)
    return animal.val_g
end

function get_phenotypes(animal::Animal;
                        h2::Union{Array{Float64}, Float64}=.5,
                        Ve::Union{Array{Float64}, Float64}=-999.99)

    n_traits = GLOBAL("n_traits")
    Vg       = GLOBAL("Vg")

    if Ve == -999.99
        if n_traits > 1 && !isa(h2, Array)
            h2 = fill(h2, n_traits)
        end
        # Handle inf variance when h2 = 0
        is_zeros = h2 .== 0
        h2[is_zeros] .= 1e-5

        Ve = ((ones(n_traits) .- h2) .* diag(Vg)) ./ h2
        Ve = n_traits == 1 ? Ve[1] : Ve
    end

    Ve = handle_diagonal(Ve, n_traits)
    animal.val_p = animal.val_g .+ cholesky(Ve).U * randn(n_traits)

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


function print(animal::Animal)
    println("ID   : ", animal.ID)
    println("Sire : ", animal.sire.ID)
    println("Dam  : ", animal.dam.ID)
    println("BV   : ", animal.val_g)
end

function Base.:+(x::Animal, y::Animal)
    return Cohort([x, y])
end
Base.show(io::IO, animal::Animal) = print(animal)
Base.iterate(animal::Animal, i...) = Base.iterate([animal], i...)


