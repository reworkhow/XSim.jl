"""
# Mating function
    mate(cohort_A         ::Cohort,
         cohort_B         ::Cohort;
         nA               ::Int64=cohort_A.n,
         nB_per_A         ::Int64=1,
         n_per_mate       ::Int64=1,
         replace_A        ::Bool =false,
         replace_B        ::Bool =false,
         ratio_malefemale ::Union{Float64, Int64}=0,
         scheme           ::String ="none",
         args...)

    mate(cohort::Cohort; args...) =  mate(cohort, cohort; args...)

## Arguments
Positional arguments
- `cohort_A` : A `cohort` object that is treated as common mating parents.
- `cohort_B` : A `cohort` object that is a mating pool from which individuals are sampled to mate with `cohort_A`.

Keyword arguments
- `nA` : `nA` individuals will be sampled from `cohort_A` and treated as common parents.
- `nB_per_A` : `nB_per_A` individuals sampled from `cohort_B` will mate with each individual from `cohort_A`.
- `n_per_mate` : `n_per_mate` progenies will be reproduced from each pair of mating parent.
- `replace_A` : Whether the sampling is replacable in `cohort_A`.
- `replace_B` : Whether the sampling is replacable in `cohort_B`.
- `ratio_malefemale` : By default, two cohorts which are male and female progenies will be returned. `ratio_malefemale` defined the progenies ratio of males over females. If `ratio_malefemale=0`, only one cohort will be returned.
- `scheme` : Available options are ["random", "diallel cross", "selfing", "DH"]. See the examples for more details.

## Outputs
By default, two `cohort` objects will be returned. The first `cohort` is assumed to be male progenies and the other `cohort` are female progenies. The size of two cohorts will folow the ratio `raiot_malefemale`. When `ratio_malefemale` is set to `0`, only one `cohort` will be returned.


## Example
### Random mating (Default)
Initialize cohorts
```jldoctest
julia> cohort_A = Cohort(5)
julia> cohort_B = Cohort(10)
```
Define mating events
```jldoctest
julia> args = Dict(:nA               => cohort_A.n,
                   :nB_per_A         => 1,
                   :replace_A        => false,
                   :replace_B        => false,
                   :n_per_mate       => 1)
julia> progenies = mate(cohort_A, cohort_B; args...)

# Equivalent
julia> progenies = mate(cohort_A, cohort_B)

# Equivalent
julia> progenies = mate(cohort_A, cohort_B; scheme="random")

# Equivalent
julia> progenies = cohort_A * cohort_B
```

Check the pedigree to see if the mating goes as desired.
```jldoctest
julia> get_pedigree(progenies)
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 19  1   8
 16  2   6
 17  3  10
 20  4  15
 18  5  14
```

### Diallel cross
Initialize cohorts
```jldoctest
julia> cohort_A = Cohort(2)
julia> cohort_B = Cohort(5)
```
Define mating events
```jldoctest
julia> args = Dict(:nA              => sires.n,
                   :nB_per_A        => dams.n,
                   :replace_A       => false,
                   :replace_B       => false,
                   :n_per_mate      => 1,
                   :ratio_malefemale=> 1)
julia> male, female = mate(cohort_A, cohort_B; args...)

# Equivalent
julia> male, female = mate(cohort_A, cohort_B; scheme="diallel cross")
```

Check the pedigree to see if the mating goes as desired.
```jldoctest
julia> get_pedigree(male + female)
10×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 12  2  7
 10  2  6
 11  2  4
 14  1  3
 15  1  5
  9  2  5
 13  1  6
 17  1  4
  8  2  3
 16  1  7
```

### Selfing
Initialize cohorts
```jldoctest
julia> parents = Cohort(5)
```

In the selfing scheme, only one `cohort` is required.
```jldoctest
julia> args = Dict(:nA          => 3,
                   :replace_A   => false,
                   :n_per_mate  => 5,
                   :scheme      => "selfing")
julia> progenies = mate(parents; args...)
```
Inspect the pedigree to verify the mating behavior
```jldoctest
julia> get_pedigree(progenies)
15×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
  6  4  4
  7  4  4
  8  4  4
  9  4  4
 10  4  4
 11  1  1
 12  1  1
 13  1  1
 14  1  1
 15  1  1
 16  5  5
 17  5  5
 18  5  5
 19  5  5
 20  5  5
```
"""
function mate(cohort_A::Union{Cohort,Animal},
    cohort_B::Union{Cohort,Animal};
    nA::Int64 = cohort_A.n,
    nB_per_A::Int64 = 1,
    n_per_mate::Int64 = 1,
    n_offspring::Int64 = -1,
    replace_A::Bool = false,
    replace_B::Bool = false,
    ratio_malefemale::Union{Float64,Int64} = 0,
    silent::Bool = GLOBAL("silent"),
    scheme::String = "none",
    args...)

    # determine scheme
    if scheme != "none"
        if scheme == "random"
            cohort = mate(cohort_A, cohort_B, silent = true)
            n_offspring = cohort_A.n

        elseif scheme == "diallel cross"
            cohort = mate(cohort_A, cohort_B; nA = cohort_A.n, nB_per_A = cohort_B.n, silent = true)
            n_offspring = cohort_A.n * cohort_B.n
            nB_per_A = cohort_B.n

        elseif scheme == "selfing"
            # Preallocate animals
            n_offspring = nA * n_per_mate
            animals = Array{Animal}(undef, n_offspring)
            # Sample parent from cohort_A
            select_A = sample(cohort_A, nA, replace = replace_A)
            # Mating
            for i = 1:nA
                parent = select_A[i]
                for j = 1:n_per_mate
                    idx = (i - 1) * n_per_mate + j
                    animals[idx] = Animal(parent, parent)
                end
            end
            cohort = Cohort(animals)
        else
            LOG("The available scheme options are: [random, diallel cross, selfing]", "error")
        end

    else
        # Preallocate animals, infer size
        animals, nA, nB_per_A, n_per_mate, n_offspring =
            preallocate_animals(nA, nB_per_A, n_per_mate, n_offspring)
        # Sample A and B breeds
        select_A = sample(cohort_A, nA, replace = replace_A)
        select_B = sample_B(cohort_B, nA, nB_per_A, replace_B)
        # Mating
        for i = 1:nA
            animal_A = select_A[i]

            for j = 1:nB_per_A
                animal_B = select_B[(i-1)*nB_per_A+j]

                for k = 1:n_per_mate
                    idx = (i - 1) * nB_per_A * n_per_mate + (j - 1) * n_per_mate + k
                    animals[idx] = Animal(animal_A, animal_B)
                end
            end
        end
        cohort = split_by_ratio(animals, ratio_malefemale, n_offspring)

    end

    log_mate(silent, n_offspring, nA, nB_per_A, n_per_mate, scheme, replace_A, replace_B)
    return cohort
end

# mate pedigree
function mate(pedigree::DataFrame)
    # cast pedigree to matrix
    ped = Matrix(pedigree)

    # collect ped info
    bol_fd = sum(ped[:, 2:end], dims = 2) .== 0
    idx_fd = ped[bol_fd[:], 1]
    bol_nfd = [~b for b in bol_fd][:]
    ped_nfd = ped[bol_nfd, :]

    # allocate cohort
    cohort_new = Array{Animal}(undef, maximum(ped))

    # fill in founders
    for i in idx_fd
        cohort_new[i] = GET_LINES([i])[1]
    end

    # fill in non-founders
    for row in eachrow(ped_nfd)
        id, sire_id, dam_id = row
        sire, dam = GET_LINES([sire_id, dam_id])
        tmp = (sire*dam)[1]
        tmp.ID = id
        cohort_new[id] = tmp
    end

    # return
    Cohort(cohort_new)
end

mate(path_ped::String) = mate(CSV.read(path_ped, DataFrame, header = false))


# intra-group mating
mate(cohort::Cohort; args...) = mate(cohort, cohort; args...)
mate(animal::Animal; args...) = mate(Cohort(animal); args...)

# handle cohort with one animal
mate(cohort_A::Animal, cohort_B::Animal; args...) =
    mate(Cohort(cohort_A), Cohort(cohort_B); args...)
mate(cohort_A::Cohort, cohort_B::Animal; args...) =
    mate(cohort_A, Cohort(cohort_B); args...)
mate(cohort_A::Animal, cohort_B::Cohort; args...) =
    mate(Cohort(cohort_A), cohort_B; args...)

Base.:*(x::Union{Cohort,Animal},
    y::Union{Cohort,Animal}) = mate(x, y)


function preallocate_animals(nA::Int64,
    nB_per_A::Int64,
    n_per_mate::Int64,
    n_offspring::Int64)

    # If no population size assigned, the rest three are assumed to be correct
    if n_offspring == -1
        n_offspring = nA * nB_per_A * n_per_mate

        # If assigned population size, but doens't match the rest
    elseif n_offspring != nA * nB_per_A * n_per_mate
        # If nB_per_A is not provided, and the rest are provided
        if nB_per_A == 1
            nB_per_A = trunc(Int, n_offspring / nA / n_per_mate) # ground to int

        # If n_per_mate is not provided, and the rest are provided
        elseif n_per_mate == 1
            n_per_mate = trunc(Int, n_offspring / nA / nB_per_A)
        end
    end

    return Array{Animal}(undef, n_offspring), nA, nB_per_A, n_per_mate, n_offspring
end

function sample_B(cohort_B,
    nA,
    nB_per_A,
    replace_B)

    n_crosses = nA * nB_per_A

    if n_crosses > cohort_B.n
        if replace_B
            # > 10x10, 5x3, T, T
            select_B = sample(cohort_B, n_crosses, replace = replace_B)

        else
            # > 10x10, 5x3, T, F
            select_B = Cohort()
            for _ = 1:nA
                select_B += sample(cohort_B, nB_per_A, replace = replace_B)
            end

        end
    else
        # > 10x10, 5x2, T, T    # > 10x10, 5x2, T, F
        select_B = sample(cohort_B, n_crosses, replace = replace_B)
    end

    return select_B

    # NOTE
    # > 10x10, 5x2, T, T
    # (10, 10, T)

    # > 10x10, 5x2, T, F
    # (10, 10, F)

    # > 10x10, 5x3, T, T
    # all crosses = 15
    # (10, 15, T)

    # > 10x10, 5x3, T, F
    # all crosses = 15
    # (10, 3, F)* 5
end

function split_by_ratio(animals,
    ratio_malefemale,
    n_offspring)

    cohort = Cohort(animals)
    # shuffle the order
    cohort = sample(cohort, cohort.n, replace = false)

    if ratio_malefemale != 0
        float_males = round(n_offspring * ratio_malefemale / (ratio_malefemale + 1))
        n_males = convert(Int64, float_males)
        return cohort[1:n_males], cohort[(n_males+1):end]

    else
        return cohort
    end
end

function log_mate(silent, n_offspring, nA, nB_per_A, n_per_mate, scheme, replace_A, replace_B)
    if !silent
        str_A = replace_A ? "with" : "without"
        str_B = replace_B ? "with" : "without"
        if scheme == "selfing"
            LOG("--------- Mating Summary ---------")
            LOG("Generated $n_offspring individuals")
            LOG("from $nA cohort_A individuals $str_A replacement")
            LOG("Each selfing parent reproduces $n_per_mate progenies")
            LOG("")
            LOG("--------- Offsprings Summary ---------")

        else
            LOG("--------- Mating Summary ---------")
            LOG("Generated $n_offspring individuals")
            LOG("from $nA cohort_A individuals $str_A replacement")
            LOG()
            LOG("Every cohort_A individual mates with $nB_per_A individuals")
            LOG("from cohort_B $str_B replacement")
            LOG()
            LOG("And each mating reproduces $n_per_mate progenies")
            LOG("")
            LOG("--------- Offsprings Summary ---------")
        end
        # print(cohort)
    end
end
