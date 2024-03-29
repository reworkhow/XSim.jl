mutable struct Animal <: AbstractAnimal
    ID::Int64
    sire::Animal
    dam::Animal
    genome_sire::Array{Chromosome,1}
    genome_dam::Array{Chromosome,1}
    val_p::Array{Float64,1}
    val_g::Array{Float64,1}
    is_founder::Bool
    has_hap::Bool # has haplotype initialized (for offspring)

    # Null constructor
    Animal() = new(0)

    # Base constructor
    function Animal(sire::Animal,
        dam::Animal;
        haplotypes::Array{AlleleIndexType} = [0],
        is_DH::Bool = false,
        is_haplotype::Bool = false)

        if sire.ID == 0 || dam.ID == 0
            is_founder = true
            has_hap = true
        else
            is_founder = false
            has_hap = false
        end

        animal = new(GLOBAL("count_id"), sire, dam,
            Array{Chromosome}(undef, GLOBAL("n_chr")),
            Array{Chromosome}(undef, GLOBAL("n_chr")),
            Array{Float64}(undef, 0),
            Array{Float64}(undef, 0),
            is_founder, has_hap)
        add_count_ID!(by = 1)

        # Setup genome
        if is_founder
            set_genome!(animal, haplotypes, is_haplotype=is_haplotype)
            add_founders!(animal)
        else
            if is_DH
                set_genome!(animal, sire)
            else
                set_genome!(animal, dam, sire)
            end
        end
        add_animal!(animal)

        return animal
    end

    # DH constructor
    function Animal(parent::Animal;
        args...)
        return Animal(parent, parent; is_DH = true)
    end
end

"""
check if the haplotypes are set
"""
function handle_init_haplotypes!(animal::Animal)
    if !animal.has_hap
        set_haplotypes!(animal)
        animal.has_hap = true
    end
end

function get_genotypes(animal::Animal)
    handle_init_haplotypes!(animal)

    myGenotype = Array{AlleleIndexType}(undef, 0)
    for i in 1:GLOBAL("n_chr")
        append!(myGenotype,
            animal.genome_sire[i].haplotype + animal.genome_dam[i].haplotype)
    end

    return myGenotype
end

function get_BVs(animal::Animal)
    handle_init_haplotypes!(animal)

    genotypes = get_genotypes(animal)
    animal.val_g = (genotypes'GLOBAL("effects")) |> transpose

    return animal.val_g
end


"""
Set genome function for founders
"""
function set_genome!(animal       ::Animal,
                     haplotypes   ::Array{AlleleIndexType,1} = [0];
                     is_haplotype ::Bool = false)
    # Founders
    is_file = haplotypes != [0]

    if is_file
        if is_haplotype
            # the phase is given, assign parental haplotypes accordingly
            n_loci = length(haplotypes)
            idx_even = [2:2:n_loci;]
            idx_odd = [1:2:n_loci;]
            hap_sire = convert(Array{AlleleIndexType}, haplotypes[idx_even])
            hap_dam = convert(Array{AlleleIndexType}, haplotypes[idx_odd])
        else
            # the phase is missing, arbitrarily assign haplotypes
            # Sire haplotype: 0->0, 1->1, 2->1
            hap_sire = convert(Array{AlleleIndexType}, haplotypes .> 0)
            # Dam haplotype:  0->0, 1->0, 2->1
            hap_dam = convert(Array{AlleleIndexType}, haplotypes .== 2)
        end
    else
        hap_sire = [0]
        hap_dam = [0]
    end

    idx_chr = GLOBAL("idx_chr")
    ori_sire = GLOBAL("count_hap")
    ori_dam = GLOBAL("count_hap") + 1
    for i in 1:GLOBAL("n_chr")
        if is_file
            from, to = idx_chr[i, :]
            animal.genome_sire[i] = Chromosome(i, ori_sire, hap_sire[from:to])
            animal.genome_dam[i] = Chromosome(i, ori_dam, hap_dam[from:to])
        else
            animal.genome_sire[i] = Chromosome(i, ori_sire, hap_sire)
            animal.genome_dam[i] = Chromosome(i, ori_dam, hap_dam)
        end
    end
    add_count_haplotype!(by = 2)
end

"""
Set genome function for progenies
"""
function set_genome!(animal::Animal, dam::Animal, sire::Animal)
    for i in 1:GLOBAL("n_chr")
        animal.genome_sire[i] = Chromosome(i, sire)
        animal.genome_dam[i] = Chromosome(i, dam)
    end
end

# For DH
function set_genome!(animal::Animal, parent::Animal)
    for i in 1:GLOBAL("n_chr")
        animal.genome_sire[i] = Chromosome(i, parent)
        animal.genome_dam[i] = deepcopy(animal.genome_sire[i])
    end
end

"""
Assign haplotypes to the animal
"""
function set_haplotypes!(animal::Animal)
    set_haplotypes!(animal.genome_sire)
    set_haplotypes!(animal.genome_dam)
end

function get_phenotypes(animal::Animal)
    nothing
end

function get_DH(animal::Animal, n::Int = 1)
    animals = Array{Animal}(undef, n)
    for i in 1:n
        animals[i] = Animal(animal)
    end
    return Cohort(animals)
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
