"""
# Initialize a cohort by population size
    Cohort(n::Int64=0)

### Arguments
- `n` : An integer to assign the population size.

### Examples
```jldoctest
julia> cohort = Cohort(5)
[ Info: Cohort (5 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.265 1.697]
[ Info:
[ Info: Variance of breeding values:
[ Info: [1.6 1.4]
```
──────────────────────────────────────────────────────────────
# Initialize a cohort by genotypes/haplotypes files
    Cohort(genetic_data ::Union{DataFrame, Array{Int64}}; args...)
    Cohort(filename     ::String; args...)

### Arguments
- `genetic_data` : A `dataframe`/`2D-array` that stores genotypes/haplotypes in the dimension of individuals by markers.
- `filename` : A `filepath` to a file storing genotypes/haplotypes data.
- `n` : Number of lines to be loaded from the file. The default value is `-1` and the entire file will be loaded.
- `random` : By default it's set to `true` to randomly select `n` lines (individuals) from the file to generate the cohort.
- `alter_maf` : It will update MAF based on the provided genotypes if it's set to `true` (default).

### Example of the `demo_genotypes.csv` and `demo_haplotypes.csv`
Both demo files store marker information for 5 individuals and 4 markers.
Use `DATA("demo_genotypes.csv")` to interact with demo files.
```
# demo_genotypes.csv
# rows: individuals, columns: markers
# Homozygote is coded as 0 and 2, otherwise is coded as 1
2,0,0,1
0,0,1,0
0,1,0,2
1,1,0,2
2,0,2,0

# demo_haplotypes.csv
# rows: individuals, columns: alleles
# Reference allele is coded as 0, otherwise is coded as 1
1,1,0,0,0,0,1,0
0,0,0,0,1,0,0,0
0,0,0,1,0,0,1,1
1,0,1,0,0,0,1,1
1,1,0,0,1,1,0,0
```

### Example
```jldoctest
# Load entire file
julia> cohort = Cohort("demo_haplotypes.csv")
julia> get_genotypes(cohort)
5×4 Array{Int64,2}:
 2  0  0  1
 2  0  2  0
 0  0  1  0
 1  1  0  2
 0  1  0  2

# Randomly load 3 individuals with a dataframe.
julia> data = DATA("demo_haplotypes.csv", header=false)
julia> cohort = Cohort(data, random=true, n=3)
julia> get_genotypes(cohort)
3×4 Array{Int64,2}:
 2  0  2  0
 0  1  0  2
 1  1  0  2

# Replace marker MAF by the provided file
julia> cohort = Cohort("demo_haplotypes.csv", alter_maf=true)
[ Info: MAF has been updated based on provided haplotypes/genotypes
[ Info: Cohort (5 individuals)
[ Info:
[ Info: Mean of breeding values:
[ Info: [1.418]
[ Info:
[ Info: Variance of breeding values:
[ Info: [2.012]
```
──────────────────────────────────────────────────────────────
# Functions that insepct `Cohort` properties:
All the listed functions can take a keyword argument `ID=true` to insert individuals' IDs as the first column.

### Genotypes
Genotype matirx in the dimension of `individuals` by `markers`
```jldoctest
julia> get_genotypes(cohort)
5×4 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
 0  0  1  0
 2  0  2  0
 2  0  0  1
 0  1  0  2
 1  1  0  2
```
### QTLs
QTLs matirx in the dimension of `individuals` by `markers`
```jldoctest
julia> get_QTLs(cohort)
5×3 Array{Int64,2}:
 2  2  0
 0  0  2
 0  1  0
 1  0  2
 2  0  1
```
### Breeding values
Breeding values in the dimenstion `individuals` by `traits`
```jldoctest
julia> get_BVs(cohort)
5×2 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
 1.26491   0.0
 3.79473   0.0
 1.26491   1.21268
 0.0       1.69775
 0.632456  1.69775
```
### Pedigree
Pedigree matrix, listed columns are in the order of individuals' ID, sire ID, and dam ID.
```jldoctest
julia> get_pedigree(cohort)
5×3 LinearAlgebra.Adjoint{Int64,Array{Int64,2}}:
1  0  0
2  0  0
3  0  0
4  0  0
5  0  0
```

### Minor Allele Frequencies (MAF)
In the case where we have 3 QTLs out of 4 markers, we want to compare their allel frequencies.

```jldoctest
julia> get_MAF(cohort)
4-element Array{Float64,1}:
 0.5
 0.2
 0.3
 0.5
```

### Phenotypes
Simulate cohort phenotypes based on the defined `phenome`. `h2` and `ve` can be assigned specifically for this one-time simulation.
```jldoctest
julia> get_phenotypes(cohort)
5×1 Array{Float64,2}:
  1.1126064336992942
 -0.8337021175232547
 -0.363621019381922
  4.042256656472762
  1.7828800511223049
```
"""
mutable struct Cohort
    animals      ::Array{Animal, 1}
    n            ::Int64

    ```Initiate cohorts with sampled genotypes```
    function Cohort(n          ::Int64=0)

        cohort = Array{Animal}(undef, n)
        for i in 1:n
            cohort[i] = Animal(Animal(), Animal())
        end

        center_BV!(cohort)
        return new(cohort, n)
    end

    ```Initiate cohorts with given genotypes```
    function Cohort(genetic_data::Union{DataFrame, Array{Int64}};
                    n           ::Int64=-1,
                    random      ::Bool=true,
                    alter_maf   ::Bool=false)

        # Extract genotypes meta
        if isa(genetic_data, DataFrame)
            n_founders    = nrow(genetic_data)
            n_loci        = ncol(genetic_data)
            genetic_data  = Array(genetic_data)

        elseif isa(genetic_data, Array)
            n_founders = size(genetic_data, 1)
            n_loci     = size(genetic_data, 2)

        else
            LOG("Input genotypes not supported", "error")
        end

        # Examine the column size
        if n_loci == GLOBAL("n_loci") * 2
            # It's haplotypes, convert it to genotypes
            for i in 1:GLOBAL("n_loci")
                i1 = (i * 2) - 1
                i2 = (i * 2)
                sub = genetic_data[:, i1:i2]
                genetic_data[:, i] = sum(sub, dims=2)
            end
            genetic_data = genetic_data[:, 1:GLOBAL("n_loci")]

        elseif n_loci != GLOBAL("n_loci")
            LOG("Number of loci mismatches the given genome", "error")
        end

        # Whether to alter MAF
        if alter_maf
            LOG("MAF has been updated based on provided haplotypes/genotypes")
            SET("maf", get_MAF(genetic_data))
        end

        if (n == -1) || (n > n_founders)
            n = n_founders
        end

        # n_founders
        cohort        = Array{Animal}(undef, n)
        pool          = random ? sample(1:n_founders, n, replace=false) : (1:n)
        for i in 1:n
            hap       = vcat(genetic_data[pool[i], :]...)
            cohort[i] = Animal(Animal(), Animal(), haplotypes=hap)
        end

        center_BV!(cohort)
        return new(cohort, n)
    end

    function Cohort(filename    ::String; args...)

        dt = CSV.read(filename, DataFrame,
                      header=false,
                      missingstrings=["9"])

        return Cohort(dt; args...)
    end

    # Constructor for non-founder
    function Cohort(animals     ::Array{Animal, 1};
                    is_center   ::Bool=false)

        n = length(animals)
        if is_center
            center_BV!(animals)
        end
        return new(animals, n)
    end

    function Cohort(animal      ::Animal;
                    is_center   ::Bool=false)

        if is_center
            center_BV!(animals)
        end
        return new([animal], 1)
    end
end

Founders(genetic_data::Union{DataFrame, Array{Int64}}; args...) = Cohort(genetic_data; args...)
Founders(filename ::String; args...)                         = Cohort(filename; args...)
Founders(n        ::Int64)                                   = Cohort(n)

function Base.summary(cohort::Cohort; is_return=true)
    bvs        = get_BVs(cohort)
    mu_g       = round.(XSim.mean(bvs, dims=1), digits=3)
    var_g      = round.(XSim.var(bvs,  dims=1), digits=3)

    if is_return
        return Dict("n"     => cohort.n,
                    "mu_g"  => mu_g,
                    "var_g" => var_g)
    else
        LOG("Cohort ($(cohort.n) individuals)")
        LOG()
        LOG("Mean of breeding values: ")
        LOG("$mu_g")
        LOG()
        LOG("Variance of breeding values: ")
        LOG("$var_g")
    end
end

get_QTLs(cohort::Cohort) = get_genotypes(cohort)[:, GLOBAL("is_QTLs")]
get_IDs(cohort::Cohort)  = (animal->animal.ID).(cohort)

function get_MAF(cohort::Cohort; ID::Bool=false)
    genotypes = get_genotypes(cohort)
    maf = get_MAF(genotypes)
    if ID == true
        return hcat(get_IDs(cohort), maf)
    else
        return maf
    end
end


function get_genotypes(cohort::Cohort, option::String="XSim"; ID::Bool=false)
    genotypes_2d_tmp = (animal->get_genotypes(animal)).(cohort)
    genotypes        = hcat(genotypes_2d_tmp...)'

    if option == "XSim"
        if ID == true
            return hcat(get_IDs(cohort), genotypes[:, :])
        else
            return genotypes[:, :]
        end

    elseif option == "JWAS"
        # dt_G = hcat(get_IDs(cohort), genotypes) |> XSim.DataFrame
        # dt_G = hcat("a".* string.(get_IDs(cohort)), genotypes) |> XSim.DataFrame
        # CSV.write("jwas_g.csv", dt_G)
        return genotypes|>DataFrame|>JWAS.get_genotypes
    end
end

function get_BVs(cohort::Cohort; ID::Bool=false)
    bv_2d = (animal->get_BVs(animal)).(cohort)
    bv_out = hcat(bv_2d...)'
    if ID == true
        return hcat(get_IDs(cohort), bv_out)
    else
        return bv_out
    end
end


function get_phenotypes(cohort   ::Cohort,
                        option   ::String="XSim";
                        h2       ::Union{Array{Float64}, Float64}=GLOBAL("h2"),
                        ve       ::Union{Array{Float64}, Float64}=GLOBAL("Ve"),
                        return_ve::Bool=false,
                        ID       ::Bool=false)

    if ve != GLOBAL("Ve")
        nothing
    elseif h2 != GLOBAL("h2")
        ve = get_Ve(GLOBAL("n_traits"), GLOBAL("Vg"), h2)
    else
        ve = handle_diagonal(ve, GLOBAL("n_traits"))
    end


    n_traits   = GLOBAL("n_traits")
    ve_u       = cholesky(ve).U
    eff_nonG   = hcat([ve_u * randn(n_traits) for _ in 1:cohort.n]...)' |> Array
    eff_G      = get_BVs(cohort)
    phenotypes = eff_G + eff_nonG
    if ID == true
        phenotypes = hcat(get_IDs(cohort), phenotypes)
    end

    if option == "XSim"
        return (return_ve) ? (phenotypes, ve) : (phenotypes)

    elseif option == "JWAS"
        jwas_P = hcat(get_IDs(cohort), phenotypes) |> XSim.DataFrame
        rename!(jwas_P, vcat(["ID"], ["y$i" for i in 1:GLOBAL("n_traits")]))
        # if nrow(cofactors) != 0
        #     jwas_P = hcat(jwas_P, cofactors)
        # end
        jwas_P[!, "ID"] = string.(Int.(jwas_P[!, "ID"]))
        return (return_ve) ? (jwas_P, ve) : (jwas_P)
    end
end

function get_pedigree(cohort::Cohort, option::String="XSim")
# Note replace 0 with Missing
    # return a 3-column matrix: ID, SireID, DamID
    ped_tmp  = (animal->[animal.ID, animal.sire.ID, animal.dam.ID]).(cohort)
    ped_array = hcat(ped_tmp...)'
    ped_array = ped_array[sortperm(ped_array[:, 1]), :]

    if option == "XSim"
        ped_array

    elseif option == "JWAS"
        # Get pedigree in dataframe
        df = ped_array |> Array |> XSim.DataFrame

        # Cast columns to strings
        df[!,1] = strip.(string.(df[!,1]))
        df[!,2] = strip.(string.(df[!,2]))
        df[!,3] = strip.(string.(df[!,3]))

        # Instantiate Pedigree
        ped = JWAS.PedModule.Pedigree(
                1,
                Dict{AbstractString, JWAS.PedModule.PedNode}(),
                Dict{Int64, Float64}(),
                Set(), Set(), Set(), Set(),
                Array{String, 1}())

        # JWAS things
        JWAS.PedModule.fillMap!(ped,df)
        for id in keys(ped.idMap)
         JWAS.PedModule.code!(ped, id)
        end

        for id in keys(ped.idMap)
          JWAS.PedModule.calcInbreeding!(ped,id)
        end

        ped.IDs = JWAS.PedModule.getIDs(ped)

        JWAS.PedModule.get_info(ped)
        JWAS.PedModule.writedlm("IDs_for_individuals_with_pedigree.txt",ped.IDs)

        return ped
    end
end


function get_DH(parents::Cohort, n::Int64=parents.n; replace=false)
    animals = Array{Animal}(undef, n)
    parents_DH = sample(parents, n, replace=replace)
    for i in 1:n
        animals[i] = get_DH(parents_DH[i])
    end
    return Cohort(animals)
end

function get_haplotype_founder(cohort::Cohort)
    haps = []
    for animal in cohort
        for chromosome in [animal.genome_sire; animal.genome_dam]
            for ori in chromosome.ori
                push!(haps, ori)
            end
        end
    end
    haps
end

function center_BV!(cohort::Array{Animal})
    if length(cohort) != 0
        bvs = [animal.val_g for animal in cohort]
        adj = XSim.mean(bvs, dims=1)
        bvs_0 = bvs .- adj

        for i in 1:length(cohort)
            cohort[i].val_g = bvs_0[i]
        end
    end
end

# ped, mme, out is from JWAS get_pedigree(), build_model(), and solve()
function putEBV(cohort::Cohort, ped, mme, out)
    # transfer ebv from mme to XSim
    trmAnimal = mme.modelTermDict["1:Animal"]
    for animal in cohort
        id = animal.ID
        strID = string(id)
        mmePos = ped.idMap[strID].seqID + trmAnimal.startPos - 1
        animal.traits[1].estimated = out[mmePos, 2]
    end
end

function sample(cohort ::Cohort,
                n      ::Int64;
                replace::Bool=true)

    select = sample(1:cohort.n, n, replace=replace)
    return Cohort(cohort.animals[select])
end


function sort(cohort::Cohort; by::String="BV")
    if by == "BV"
        bvs = get_BVs(cohort)
        bvs_all = XSim.sum(bvs, dims=2)[:, 1]
        return cohort[sortperm(bvs_all, rev=true)]
    elseif by == "pedigree"
        ids = get_IDs(cohort)
        return cohort[sortperm(ids)]
    else
        LOG("Available options are ['BV', 'pedigree']")
        return cohort
    end
end



function print(cohort::Cohort, option::String="None")
    if option == "None"
        if cohort.n != 0
            Base.summary(cohort, is_return=false)
        end

    elseif option == "ID"
        print("Individual: [ ")
        for animal in cohort
             print(animal.ID, " ")
        end
        println("]")

    elseif option == "Pedigree"
        return get_pedigree(cohort)
    end
end

function getindex(cohort::Cohort, I...)
    if length(I...) == 1
        return cohort.animals[I[1]]
    else
        return Cohort(getindex(cohort.animals, I...))
    end
end

Base.:+(x::Cohort, y::Cohort)       = Cohort(vcat(x.animals, y.animals))
Base.:+(x::Cohort, y::Animal)       = Cohort(vcat(x.animals, y))
Base.:+(x::Animal, y::Cohort)       = Cohort(vcat(x, y.animals))
Base.length(cohort::Cohort)         = length(cohort.animals)
Base.show(io::IO, cohort::Cohort)   = GLOBAL("silent") ? nothing : print(cohort)
Base.iterate(cohort::Cohort, i...)  = Base.iterate(cohort.animals, i...)
Base.lastindex(cohort::Cohort)      = length(cohort)
Base.setindex!(cohort::Cohort, animal::Animal, i::Int64) =
    Base.setindex!(cohort.animals, animal, i)




# function get(cohort::Cohort,
#              item  ::String,
#              option::Any)

#     if item == "Traits"
#         return get_phenotypes(cohort, option)

#     elseif item == "ID"
#         return get_IDs(cohort)

#     elseif item == "Pedigree"
#         return get_pedigree(cohort)

#     elseif item == "DH"
#         return get_DH(cohort, option)

#     else
#         println("""
#             The available options are: 'Triats', 'ID', 'Pedigree', and 'DH'
#         """)
#     end
# end/